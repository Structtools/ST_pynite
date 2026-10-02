"""
Derivation of the 12-DOF beam geometric stiffness matrix used by `Member3D.kg`, including the
bending moment, shear force and torque terms that produce lateral-torsional and flexural-torsional
buckling.

Run with Sympy installed (`uv run --extra derivations python "Derivations/Geometric Stiffness
Matrix - Moments and Torque.py"`). The script prints the matrix split into its axial force part,
which is the classical consistent geometric stiffness matrix (McGuire, Gallagher & Ziemian,
"Matrix Structural Analysis", 2nd Ed., Section 9.3), and the moment/torque part.

Kinematics and sign conventions (Pynite local axes)
---------------------------------------------------
x runs along the member from the i-end to the j-end; y and z complete a right-handed system.

DOF order: [u_i, v_i, w_i, θx_i, θy_i, θz_i, u_j, v_j, w_j, θx_j, θy_j, θz_j]

The displacement of a point (y, z) of the cross-section of a thin-walled beam, with warping
neglected and the shear centre at the centroid, is

    u_p = u - y v' - z w',     v_p = v - z φ,     w_p = w + y φ

with φ = θx the twist, v' = +θz and w' = -θy (matching the sign pattern of Pynite's elastic
stiffness matrix). The initial stresses of the pre-buckling state are

    σ_x  = N/A - M_z y/I_z + M_y z/I_y
    τ_xy, τ_xz : transverse shear stresses with resultants V_y, V_z and torque T

where N (tension positive), M_y(x), M_z(x) and T are the stress resultants on the +x face.

Second-order work of the initial stresses
------------------------------------------
The nonlinear parts of the Green-Lagrange strains are

    ε_xx^(2) = ½ (v_p,x² + w_p,x²)
    γ_xy^(2) = u_p,x u_p,y + v_p,x v_p,y + w_p,x w_p,y
    γ_xz^(2) = u_p,x u_p,z + v_p,x v_p,z + w_p,x w_p,z

Integrating σ_x ε^(2) + τ_xy γ_xy^(2) + τ_xz γ_xz^(2) over a doubly symmetric section (∫y dA =
∫z dA = ∫yz dA = 0, ∫y² dA = I_z, ∫z² dA = I_y) and dropping the products of u' with the lateral
slopes gives

    U_g = ½ ∫ [ N (u'² + v'² + w'²) + N (I_p/A) φ'²
                - 2 M_z φ' w' - 2 M_y φ' v'
                + 2 V_y w' φ  - 2 V_z v' φ
                + T (v'' w' - w'' v') ] dx

The torque term uses the equal split of the torque between ∫y τ_xz dA and -∫z τ_xy dA, which
is the semitangential (Argyris) assumption. Nothing here is integrated by parts, so no boundary
terms are dropped and the element needs no additional joint moment terms at non-collinear
joints: summing the elements gives the complete second-order strain energy of the frame.

Stress resultants in terms of the end actions on the element
------------------------------------------------------------
With the member's local end force vector f = [fx_i, fy_i, fz_i, mx_i, my_i, mz_i, fx_j, ...]
(the actions on the member's ends), and the moments varying linearly along the element,

    N      = fx_j                                   (= -fx_i)
    M_z(x) = -mz_i (1 - ξ) + mz_j ξ,  ξ = x/L
    M_y(x) = -my_i (1 - ξ) + my_j ξ
    V_y    = -dM_z/dx = -(mz_i + mz_j)/L
    V_z    = +dM_y/dx = (my_i + my_j)/L
    T      = mx_j

Interpolation: linear for u and φ, cubic Hermite for v and w (as in the elastic matrix).

Verification
------------
The resulting matrix reproduces, with 16 elements and to within 0.2 %, the classical closed-form
solutions for lateral-torsional buckling without warping (Timoshenko & Gere, "Theory of Elastic
Stability", Chapter 6): simply supported beam under uniform moment M_cr = (π/L)√(EI_y GJ);
simply supported beam with a central point load P_cr = 16.94 √(EI_y GJ)/L²; simply supported beam
with a uniform load q_cr = 28.3 √(EI_y GJ)/L³; cantilever with a tip load P_cr = 4.013 √(EI_y GJ)/L².
See `Testing/test_buckling.py`.
"""
import sympy as sp

x, L = sp.symbols('x L', positive=True)
N, Myi, Mzi, Myj, Mzj, T, IpA = sp.symbols('N M_yi M_zi M_yj M_zj T Ip/A', real=True)
xi = x/L

# Hermite cubics in xi (the rotation shape functions carry a factor L)
H1 = 1 - 3*xi**2 + 2*xi**3
H2 = (xi - 2*xi**2 + xi**3)*L
H3 = 3*xi**2 - 2*xi**3
H4 = (-xi**2 + xi**3)*L
N1 = 1 - xi
N2 = xi

# Shape function rows for u, v, w and phi over the 12 DOFs
Nu = [N1, 0, 0, 0, 0, 0, N2, 0, 0, 0, 0, 0]
Nv = [0, H1, 0, 0, 0, H2, 0, H3, 0, 0, 0, H4]       # v' = +theta_z
Nw = [0, 0, H1, 0, -H2, 0, 0, 0, H3, 0, -H4, 0]     # w' = -theta_y
Nphi = [0, 0, 0, N1, 0, 0, 0, 0, 0, N2, 0, 0]


def d(row, n=1):
    """Differentiates a row of shape functions n times with respect to x."""
    return [sp.diff(f, x, n) for f in row]


def outer(a, b):
    """Returns the 12x12 matrix a_i b_j."""
    return sp.Matrix(12, 12, lambda i, j: a[i]*b[j])


def sym(a, b):
    """Returns the symmetric matrix of the bilinear form (a.d)(b.d)."""
    return (outer(a, b) + outer(b, a))/2


# Stress resultants along the element
Mz = -Mzi*(1 - xi) + Mzj*xi
My = -Myi*(1 - xi) + Myj*xi
Vy = -sp.diff(Mz, x)
Vz = sp.diff(My, x)

# Integrand of the second-order work (the leading 1/2 is dropped to give the stiffness matrix)
integrand = (
    N*(outer(d(Nu), d(Nu)) + outer(d(Nv), d(Nv)) + outer(d(Nw), d(Nw)))
    + N*IpA*outer(d(Nphi), d(Nphi))
    - 2*Mz*sym(d(Nphi), d(Nw))
    - 2*My*sym(d(Nphi), d(Nv))
    + 2*Vy*sym(d(Nw), Nphi)
    - 2*Vz*sym(d(Nv), Nphi)
    + T*(sym(d(Nv, 2), d(Nw)) - sym(d(Nw, 2), d(Nv)))
)

kg = sp.Matrix(12, 12, lambda i, j: sp.simplify(sp.integrate(integrand[i, j], (x, 0, L))))

assert kg == kg.T, 'The geometric stiffness matrix must be symmetric'

# Split into the axial force part and the moment/torque part
kg_N = kg.subs({Myi: 0, Mzi: 0, Myj: 0, Mzj: 0, T: 0})
kg_M = sp.simplify(kg - kg_N)

print('Axial force part, multiplied by L/N (the classical consistent geometric stiffness matrix):')
sp.pprint(sp.simplify(kg_N*L/N))
print()
print('Bending moment, shear force and torque part (upper triangle, 1-based indices):')
for i in range(12):
    for j in range(i, 12):
        if kg_M[i, j] != 0:
            print(f'  ({i + 1:2d},{j + 1:2d}): {sp.factor(kg_M[i, j])}')
