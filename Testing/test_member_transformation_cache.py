from numpy import array_equal
import pytest as pt

from Pynite import FEModel3D


def _model():

    model = FEModel3D()
    model.add_node('N1', 0, 0, 0)
    model.add_node('N2', 10, 2, 3)
    model.add_material('Steel', 29000/12**2, 11200/12**2, 0.3, 0.490, 60)
    model.add_section('W12x26', 7.65/12**2, 17.3/12**4, 204/12**4, 0.3/12**4)
    model.add_member('M1', 'N1', 'N2', 'Steel', 'W12x26', 30)
    return model


def test_transformation_matrix_is_built_once():

    member = _model().members['M1']
    assert member.T() is member.T()


def test_transformation_matrix_is_read_only():

    member = _model().members['M1']
    with pt.raises(ValueError):
        member.T()[0, 0] = 1.0


def test_moving_a_node_rebuilds_the_matrix():

    model = _model()
    member = model.members['M1']
    before = member.T()
    model.nodes['N2'].Y = 5
    after = member.T()
    assert not array_equal(before, after)
    assert array_equal(after, member._build_T())


def test_changing_the_rotation_rebuilds_the_matrix():

    member = _model().members['M1']
    before = member.T()
    member.rotation = 60
    after = member.T()
    assert not array_equal(before, after)
    assert array_equal(after, member._build_T())
