# The shortcuts that make large models faster must return exactly what the code they replace
# returned: the same floats, bit for bit, so that no result moves. Each test compares a shortcut
# with the general computation it skips, written out here as it stood before the shortcut.

from numpy import array, insert, matmul, subtract, zeros
from numpy.linalg import inv
import pytest as pt

from Pynite import FEModel3D
from Pynite import Analysis

COMBO = 'Combo 1'


def _frame():
    """A loaded, solved frame out of plane, one of its members pinned at both ends."""

    model = FEModel3D()
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', 0, 12, 0)
    model.add_node('N3', 16, 13.5, 0.7)
    model.add_node('N4', 16, 0, 0)
    model.add_material('Steel', 29000*144, 11200*144, 0.3, 0.490, 50)
    model.add_section('W12x26', 7.65/144, 204/12**4, 17.3/12**4, 0.3/12**4)
    model.add_member('M1', 'N1', 'N2', 'Steel', 'W12x26')
    model.add_member('M2', 'N2', 'N3', 'Steel', 'W12x26', rotation=30)
    model.add_member('M3', 'N4', 'N3', 'Steel', 'W12x26')
    model.add_member('M4', 'N1', 'N3', 'Steel', 'W12x26')
    model.def_releases('M4', Ryi=True, Rzi=True, Ryj=True, Rzj=True)
    model.def_support('N1', True, True, True, True, True, True)
    model.def_support('N4', True, True, True, False, False, False)
    model.add_member_dist_load('M2', 'Fy', -1.3, -0.9, 2, 14)
    model.add_member_pt_load('M2', 'Fz', 0.7, 5)
    model.add_member_pt_load('M1', 'Mx', 0.2, 4)
    model.add_node_load('N2', 'FX', 3.1)
    model.add_node_load('N3', 'FZ', -1.7)
    model.analyze(check_statics=False)
    return model


def _sub_members(model):
    return [sub for phys in model.members.values() for sub in phys.sub_members.values()]


def _bits_equal(a, b):
    """Same shape, dtype and bytes: -0.0 and 0.0 differ, and so does any last bit."""

    a, b = array(a), array(b)
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


def _expand(matrix, releases, vector=False):
    """Re-inserts zero rows (and columns) at the released DOFs, as the condensation does."""

    for i, released in enumerate(releases):
        if released:
            matrix = insert(matrix, i, 0, axis=0)
            if not vector:
                matrix = insert(matrix, i, 0, axis=1)
    return matrix


def test_ke_is_the_condensed_matrix():

    for member in _sub_members(_frame()):
        k11, k12, k21, k22 = member._partition(member._ke_unc())
        expected = _expand(subtract(k11, matmul(matmul(k12, inv(k22)), k21)), member.Releases)
        assert _bits_equal(member.ke(), expected), member.name


def test_fer_is_the_condensed_vector():

    for member in _sub_members(_frame()):
        k11, k12, k21, k22 = member._partition(member._ke_unc())
        fer1, fer2 = member._partition(member._fer_unc(COMBO))
        expected = _expand(subtract(fer1, matmul(matmul(k12, inv(k22)), fer2)), member.Releases,
                           vector=True)
        assert _bits_equal(member.fer(COMBO), expected), member.name


@pt.mark.parametrize('elastic_condensation', [False, True])
@pt.mark.parametrize('P', [0.0, -12.5, 7.25])
def test_kg_without_releases_is_the_condensed_matrix(P, elastic_condensation):

    moments = dict(Myi=1.5, Mzi=-0.25, Myj=-2.0, Mzj=0.5, Mxj=0.125)
    for member in _sub_members(_frame()):
        if any(member.Releases):
            continue

        # The uncondensed matrix `kg` partitions, read off a call that has to condense
        captured = []
        original = member._partition
        member._partition = lambda matrix: captured.append(matrix) or original(matrix)
        member.Releases = [False]*11 + [True]
        try:
            member.kg(P, combo_name=COMBO, **moments)
        finally:
            member.Releases = [False]*12
            del member._partition
        kg11, kg12, kg21, kg22 = original(captured[0])

        if elastic_condensation:
            expected = kg11
        elif P == 0.0:
            expected = zeros(kg11.shape)
        else:
            expected = subtract(kg11, matmul(matmul(kg12, inv(kg22)), kg21))
        result = member.kg(P, elastic_condensation=elastic_condensation, combo_name=COMBO,
                           **moments)
        assert _bits_equal(result, expected), member.name


def test_cross_product_matches_numpy():

    from numpy import cross
    from numpy.random import default_rng
    from Pynite.Member3D import _cross

    # `_build_T` always crosses a vector of floats, with floats or with whole numbers
    rng = default_rng(3)
    floats = [list(v) for v in rng.normal(size=(200, 3)) * 10.0**rng.integers(-6, 6, (200, 1))]
    floats += [[0.0, -0.0, 1.0], [-0.0, 0.0, -0.0]]
    others = floats[::7] + [[0, 1, 0], [-1, 0, 0], [3, 0, -4]]
    for a in floats:
        for b in others:
            assert _bits_equal(_cross(a, b), cross(a, b))
            assert _bits_equal(_cross(b, a), cross(b, a))
            assert _bits_equal(_cross(array(a), b), cross(array(a), b))


def test_inverse_transformation_matrix_is_cached_with_the_matrix():

    member = _sub_members(_frame())[1]
    assert member._T_inv() is member._T_inv()
    assert _bits_equal(member._T_inv(), inv(member.T()))
    with pt.raises(ValueError):
        member._T_inv()[0, 0] = 1.0

    # Moving a node rebuilds the inverse along with the matrix
    before = member._T_inv()
    member.j_node.Y += 1.0
    assert member._T_inv() is not before
    assert _bits_equal(member._T_inv(), inv(member._build_T()))


def _unpartition_by_search(model, V1, V2, V1_indices, V2_indices):
    """`_unpartition` as it stood: a search of both lists for every DOF."""

    V = zeros((len(model.nodes)*6, 1))
    for node in model.nodes.values():
        for i in range(6):
            if node.ID*6 + i in V2_indices:
                V[(node.ID*6 + i, 0)] = V2[V2_indices.index(node.ID*6 + i), 0]
            else:
                V[(node.ID*6 + i, 0)] = V1[V1_indices.index(node.ID*6 + i), 0]
    return V


def test_unpartition_matches_the_search():

    model = _frame()
    D1_indices, D2_indices, _ = Analysis._partition_D(model)
    D1 = array([[0.5*k - 3.25] for k in range(len(D1_indices))])
    D2 = array([[-1.0 - k] for k in range(len(D2_indices))])
    assert _bits_equal(Analysis._unpartition(model, D1, D2, D1_indices, D2_indices),
                       _unpartition_by_search(model, D1, D2, D1_indices, D2_indices))


def test_unpartition_still_raises_for_a_dof_in_neither_list():

    model = _frame()
    D1_indices, D2_indices, D2 = Analysis._partition_D(model)
    D1 = zeros((len(D1_indices) - 1, 1))
    with pt.raises(ValueError):
        Analysis._unpartition(model, D1, D2, D1_indices[1:], D2_indices)


def test_geometric_stiffness_assembly_keeps_the_term_order():

    model = _frame()
    Kg = model.Kg(COMBO, first_step=False, sparse=True, include_moments=True)

    # The term-by-term assembly as it stood, on the same member matrices
    row, col, data = [], [], []
    for member in _sub_members(model):
        f = (member.ke() @ member.d(COMBO) + member.fer(COMBO)).ravel()
        member_Kg = member.Kg(0.5*(f[6] - f[0]), elastic_condensation=True, combo_name=COMBO,
                              Myi=f[4], Mzi=f[5], Myj=f[10], Mzj=f[11], Mxj=f[9])
        for a in range(12):
            m = member.i_node.ID*6 + a if a < 6 else member.j_node.ID*6 + (a - 6)
            for b in range(12):
                n = member.i_node.ID*6 + b if b < 6 else member.j_node.ID*6 + (b - 6)
                row.append(m)
                col.append(n)
                data.append(member_Kg[(a, b)])

    # scipy picks the index type of the COO matrix itself, so the indices are compared as numbers
    assert Kg.row.tolist() == row
    assert Kg.col.tolist() == col
    assert _bits_equal(Kg.data, array(data))


def test_nodal_instability_is_still_reported():

    # Nothing holds N1's rotation about X: the member's end is released and the support is free
    model = FEModel3D()
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', 10, 0, 0)
    model.add_material('Steel', 29000, 11200, 0.3, 0.490)
    model.add_section('Section', 10, 100, 150, 250)
    model.add_member('M1', 'N1', 'N2', 'Steel', 'Section')
    model.def_releases('M1', Rxi=True)
    model.def_support('N1', True, True, True, False, True, True)
    model.def_support('N2', True, True, True, True, True, True)
    with pt.raises(Exception, match='Unstable node'):
        model.analyze(check_statics=False)
