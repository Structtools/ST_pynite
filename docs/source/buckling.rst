=================
Buckling Analysis
=================

Linear buckling analysis finds the load factor α_cr at which a structure in its first-order
state loses stability, together with the buckling mode shapes. Pynite provides
``FEModel3D.analyze_buckling()``, which solves the eigenvalue problem

.. math::

   [K_e]\{\phi\} = \lambda\,[-K_g]\{\phi\}

where :math:`[K_e]` is the elastic stiffness matrix and :math:`[K_g]` the geometric stiffness
matrix assembled from the member forces of a first-order static solve. The smallest positive
:math:`\lambda` is the critical load factor: the loads of the analysed combination multiplied by
:math:`\lambda` reach the elastic critical load.

Frame elements only. Plates and quads contribute no geometric stiffness.

What is captured
================

The member geometric stiffness matrix contains two groups of terms (see ``Member3D.kg``):

- **Axial force terms**: the classical consistent geometric stiffness matrix of a 3D beam. These
  give flexural (Euler) buckling about both axes and, through the :math:`I_p/A` term, torsional
  buckling of doubly symmetric sections.
- **Bending moment, shear force and torque terms**: these couple the lateral displacement, the
  twist and the bending rotations, and give lateral-torsional buckling of beams and
  flexural-torsional buckling of beam-columns. They are what makes an out-of-plane analysis of a
  frame loaded in its plane meaningful. They are included by default and can be switched off with
  ``include_moments=False``, which reproduces the axial-force-only behaviour of earlier versions.
- **Load height terms**: a transverse member load applied at a distance :math:`z_g` from the
  shear centre moves by :math:`z_g\varphi^2/2` along its line of action when the section twists
  by :math:`\varphi`, and so adds :math:`-|p|\,z_g` to the twist stiffness. Give the height with
  ``add_member_pt_load(..., load_height=z_g)`` or ``add_member_dist_load(..., load_height=z_g)``.
  Following EN 1993-1-1 Annex F, :math:`z_g` is positive when the load acts towards the shear
  centre (a gravity load on the top flange, destabilising) and negative when it acts away from it
  (a load hanging from the bottom flange, or uplift on the top flange, stabilising). The
  convention is tied to the load's line of action, so it does not depend on the member's
  rotation or on whether the load is given in local or global coordinates. Load height only
  affects the buckling analysis; static results are unchanged. Nodal loads act at the node.

The moment terms are derived from the second-order work of the initial stresses with the
Green-Lagrange strain of a doubly symmetric thin-walled beam, with the moment varying linearly
along each element. The derivation is in
``Derivations/Geometric Stiffness Matrix - Moments and Torque.py``.

Planar frames
=============

A frame modelled in a plane (say XY) can be analysed three ways. The examples below use a
portal frame; see ``Examples/portal_frame_buckling.py`` for the complete model.

.. code-block:: python

   # In-plane buckling: out-of-plane DOFs (DZ, RX, RY) are restrained at every node,
   # including the interior nodes created by the subdivision. α_cr for EN 1993-1-1 5.2.1.
   in_plane = model.analyze_buckling(num_modes=3, plane='XY')

   # Out-of-plane buckling: the static pre-solve uses the supports as modelled, so the
   # model must include its out-of-plane supports. The eigenproblem then restrains the
   # in-plane DOFs (DX, DY, RZ), so only lateral and lateral-torsional modes come back.
   # This is α_cr,op for the general method of EN 1993-1-1 6.3.4.
   out_of_plane = model.analyze_buckling(num_modes=3, plane='XY', out_of_plane=True)

   # Full 3D analysis with the supports as modelled. Each mode can be labelled.
   full = model.analyze_buckling(num_modes=6)
   labels = full.classify_modes('XY')   # 'in-plane', 'out-of-plane' or 'coupled'

For a planar frame loaded in its plane the in-plane and out-of-plane degrees of freedom are
uncoupled, so the modes of the 3D analysis are exactly the union of the modes of the two
restricted analyses. The restricted analyses are the convenient way to get a given number of
modes of one kind; for an unbraced steel frame the out-of-plane modes usually come first and
would otherwise crowd out the in-plane ones.

.. note::

   Pynite's default member orientation puts the local z-axis, and so ``Iz``, in the plane of an
   XY frame. If the section properties follow the Eurocode convention (``Iy`` strong axis), add
   the members with ``rotation=90`` so that the strong axis bends in the plane of the frame.
   ``effective_length(..., plane='y')`` then refers to in-plane bending of such a member.

Member subdivision
==================

Buckling modes are curved, and one cubic element per member overestimates the Euler load of a
pin-ended member by 22 %. Every physical member is therefore subdivided into
``elements_per_member`` elements (default 8) for the duration of the analysis. The subdivision is
made on an internal copy of the model, so the caller's model is untouched apart from the results
written back to it. With 8 elements per member the benchmark problems below are within 2 % of
the closed-form values, with 16 within 0.5 %.

Results
=======

``analyze_buckling`` returns a ``BucklingResults`` object with:

- ``load_multipliers`` (alias ``critical_load_factors``): the load factors α_cr, ascending.
- ``mode_shapes``, ``free_dof_indices`` and ``dof_map``: the mode shapes in the free-DOF basis
  of the eigenproblem, scaled to a largest component of 1, and the node and DOF of each row.
- ``effective_length(member, mode, plane)``: the flexural buckling length
  :math:`L_{cr} = \pi\sqrt{EI/(\alpha_{cr} N_{Ed})}` of a member, for EN 1993-1-1 6.3.1 checks.
  Only meaningful for a compressed member that takes part in the mode.
- ``out_of_plane_share(plane)`` and ``classify_modes(plane)``: how much of each mode lives on
  the out-of-plane DOFs, and a label for each mode.
- ``analysis_model``: the subdivided copy the analysis ran on, holding the first-order results.

The mode shapes are also stored on the model as the displacements of load combinations
``Buckling Mode 1``, ``Buckling Mode 2``, ... (tagged ``'buckling'``), so they render like any
deformed shape, and the first-order displacements of the analysed combination are stored on
the model too, so member forces can be queried.

Verification
============

``Testing/test_buckling.py`` checks the implementation against closed-form solutions. Those for
lateral-torsional buckling are from Timoshenko & Gere, *Theory of Elastic Stability*, Chapter 6,
all for :math:`I_w = 0`:

===================================================  ======================================  ==========
Case                                                 Closed form                             Error (16 elements)
===================================================  ======================================  ==========
Pin-ended column                                     :math:`\pi^2 EI/L^2`                    < 0.1 %
Cantilever column                                    :math:`\pi^2 EI/(2L)^2`                 < 0.1 %
Rigid-beam sway frame, fixed / pinned bases          :math:`K = 1.0` / :math:`K = 2.0`       < 1 %
Fork-supported beam, uniform moment                  :math:`(\pi/L)\sqrt{EI_y GJ}`           0.16 %
Fork-supported beam, central point load              :math:`16.94\sqrt{EI_y GJ}/L^2`         0.18 %
Fork-supported beam, uniform load                    :math:`28.3\sqrt{EI_y GJ}/L^3`          0.5 %
Cantilever beam, tip load                            :math:`4.013\sqrt{EI_y GJ}/L^2`         0.09 %
Laterally clamped beam, uniform moment               :math:`(2\pi/L)\sqrt{EI_y GJ}`          0.64 %
Fork-supported beam, uniform moment, with warping    :math:`(\pi/L)\sqrt{EI_y(GJ+\pi^2EI_w/L^2)}`  0.16 %
Rigid cantilever, tip load at height :math:`z_g`     :math:`GJ/(L z_g)` (pure twist)           exact
Rigid cantilever, uniform load at height :math:`z_g` :math:`(\pi/2)^2 GJ/(L^2 z_g)`            0.08 %
===================================================  ======================================  ==========

The load height terms are also checked against Timoshenko & Gere's first-order correction for a
cantilever, :math:`P_{cr}(a) \approx P_{cr}(0)\,(1 - a\sqrt{EI/GJ}/L)`, and against the
EN 1993-1-1 Annex F three-factor formula (:math:`C_1 = 1.365`, :math:`C_2 = 0.553`) for a
fork-supported beam with a central point load on the top and bottom flange, which agree to
within the accuracy of those formulas (0.1 % and 3 %).

The in-plane results for several portal frames have also been compared with FEM Design; see
``Testing/comparison_summary.md``.

Limitations
===========

Warping
   The 6-DOF beam element has no warping degree of freedom, so warping stiffness is neglected by
   default. This is conservative for I-sections and the difference is significant for short
   spans. ``warping='equivalent_torsion'`` replaces the torsion constant of every member in the
   eigenproblem by

   .. math::

      J_{eff} = J + \frac{\pi^2 E I_w}{G L^2}

   with :math:`L` the length of the physical member. This is exact for a member in uniform
   bending with fork supports and free warping at its ends (:math:`k = k_w = 1` over the member
   length) and approximate otherwise. It requires ``Iw`` on every section (``add_section(...,
   Iw=...)``). A member whose twist is restrained at interior nodes should be split into separate
   physical members; otherwise the longer length understates the warping contribution and the
   approximation stays on the safe side.

   Without warping, the torsional buckling load of a column is :math:`G J A / I_p` for every
   twist mode regardless of its length, so a column braced closely about its weak axis may show
   a cluster of identical torsional eigenvalues. That is the correct consequence of
   :math:`I_w = 0`, not a numerical artifact.

Load height
   Member loads act at the shear centre unless a ``load_height`` is given. Nodal loads always act
   at the node; apply a load with a height as a member point load at the member end instead.

Sections
   Sections are taken as doubly symmetric with the shear centre at the centroid.

Moment variation
   The bending moment varies linearly over each element. Members carrying distributed loads rely
   on the subdivision for accuracy.

Torque
   The torque terms follow the same derivation with the torque split equally between the two
   shear stress first moments (the semitangential assumption). They are zero for a planar frame
   loaded in its plane. External nodal moments are treated as having no second-order work.

Non-linear members
   Tension-only and compression-only members are treated as active in the first-order pre-solve.
