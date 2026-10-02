# Linear buckling analysis of a planar portal frame: in-plane and out-of-plane
#
# A 2D portal frame is modelled in the global XY plane. The same model is analysed three ways:
#
#   1. In-plane buckling (plane='XY'): the out-of-plane DOFs are restrained at every node. This is
#      the classical sway/member buckling of the frame and gives α_cr for EN 1993-1-1 5.2.1.
#   2. Out-of-plane buckling (plane='XY', out_of_plane=True): lateral and lateral-torsional
#      buckling of the members under the in-plane forces and moments. This is α_cr,op for the
#      general method of EN 1993-1-1 6.3.4. The model must include its out-of-plane supports.
#   3. Full 3D analysis: both at once. Each mode is labelled in-plane or out-of-plane.
#
# Units: N and m.

from Pynite import FEModel3D

# IPE200. Iy is the strong axis in Eurocode notation. The warping constant Iw is only needed for
# the equivalent torsion constant approximation of warping shown below.
E, G = 210e9, 80.77e9
A, Iy, Iz, J, Iw = 28.48e-4, 1943e-8, 142e-8, 7.0e-8, 12990e-12

h = 3.0     # Column height
b = 5.0     # Beam span
w = -10e3   # Uniform load on the beam (N/m, downward)

def build_frame(brace_eaves: bool) -> FEModel3D:

    model = FEModel3D()
    model.add_material('Steel', E, G, 0.3, 0.0)
    model.add_section('IPE200', A, Iy, Iz, J, Iw=Iw)

    model.add_node('A', 0, 0, 0)
    model.add_node('B', 0, h, 0)
    model.add_node('C', b, h, 0)
    model.add_node('D', b, 0, 0)

    # Pynite's default member orientation puts the local z-axis (and so Iz) in the plane of an XY
    # frame. Rotating the members by 90° puts the strong axis Iy in the plane instead.
    model.add_member('Col1', 'A', 'B', 'Steel', 'IPE200', rotation=90)
    model.add_member('Beam', 'B', 'C', 'Steel', 'IPE200', rotation=90)
    model.add_member('Col2', 'D', 'C', 'Steel', 'IPE200', rotation=90)

    # Fixed bases, in and out of plane
    model.def_support('A', True, True, True, True, True, True)
    model.def_support('D', True, True, True, True, True, True)

    # Optional lateral (DZ) and torsional (RX) restraint at the eaves, e.g. from roof bracing
    if brace_eaves:
        model.def_support('B', False, False, True, True, False, False)
        model.def_support('C', False, False, True, True, False, False)

    model.add_member_dist_load('Beam', 'FY', w, w)

    return model

for brace_eaves in (False, True):

    model = build_frame(brace_eaves)
    print(f"\n=== Portal frame, eaves {'braced' if brace_eaves else 'unbraced'} out of plane ===")

    # 1. In-plane buckling. Members are subdivided into 8 elements each for the analysis.
    in_plane = model.analyze_buckling(num_modes=3, plane='XY')
    print('In-plane α_cr        :', ', '.join(f'{lam:8.2f}' for lam in in_plane.load_multipliers))

    # Buckling length of the columns from the first in-plane mode. With rotation=90 the in-plane
    # bending of these members is about their local y-axis.
    L_cr = in_plane.effective_length('Col1', mode=0, plane='y')
    N_Ed = abs(model.members['Col1'].axial(0.0, 'Combo 1'))
    print(f'Column: N_Ed = {N_Ed/1e3:.1f} kN, L_cr = {L_cr:.2f} m = {L_cr/h:.2f} h')

    # 2. Out-of-plane buckling: lateral-torsional buckling of the beam and weak-axis buckling of the
    #    columns under the in-plane forces. Warping is neglected (conservative).
    out_of_plane = model.analyze_buckling(num_modes=3, plane='XY', out_of_plane=True)
    print('Out-of-plane α_cr,op :', ', '.join(f'{lam:8.2f}' for lam in out_of_plane.load_multipliers))

    # The same with warping approximated by an equivalent torsion constant J + π²EIw/(GL²)
    with_warping = model.analyze_buckling(num_modes=3, plane='XY', out_of_plane=True,
                                          warping='equivalent_torsion')
    print('  with warping approx:', ', '.join(f'{lam:8.2f}' for lam in with_warping.load_multipliers))

    # 3. Full 3D analysis: the modes of 1 and 2 together, each labelled
    full = model.analyze_buckling(num_modes=4)
    for i, (lam, label) in enumerate(zip(full.load_multipliers, full.classify_modes('XY'))):
        print(f'3D mode {i + 1}: α_cr = {lam:8.2f}  ({label})')

# The mode shapes are stored on the model as load combinations 'Buckling Mode 1', 'Buckling Mode 2',
# ... and can be rendered like any deformed shape:
#
# from Pynite.Visualization import Renderer
# rndr = Renderer(model)
# rndr.deformed_shape = True
# rndr.deformed_scale = 0.5
# rndr.combo_name = 'Buckling Mode 1'
# rndr.render_model()
