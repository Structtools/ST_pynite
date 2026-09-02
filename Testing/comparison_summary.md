# Pynite reference-comparison run

- Run at: 2026-08-11 12:31:18
- Repository: `C:\Users\AsgerAndersen\Desktop\ST_pynite`
- Revision: 8cad571 (fix/timoshenko-fer)
- Python: 3.12.12
- Result: 4 passed, 0 failed of 4 script(s) in 3.7s

These scripts print analysis results intended for comparison against an external reference (PolyFrame, FEM Design, hand calculations). The full console output of each run is reproduced below.

## Status

| Script | Status | Time | What it models |
| --- | --- | --- | --- |
| `braced_pitched_frame_comparison.py` | passed | 1.0s | Pitched-roof portal frame with a diagonal brace (A -> C) |
| `buckling_comparison.py` | passed | 1.0s | Two-bay single-storey portal frame: column and beam buckling scenarios |
| `hinged_frame_comparison.py` | passed | 0.9s | Single-bay portal frame with a moment release at the top-left corner |
| `pitched_frame_comparison.py` | passed | 0.8s | Single-bay pitched-roof (gable) portal frame |

## Output

### braced_pitched_frame_comparison.py

Pitched-roof portal frame with a diagonal brace (A -> C). Status: **passed** (1.0s).

```text
================================================================
  PyNite -- Braced Pitched-Roof Portal Frame Buckling Analysis
================================================================

  Scenario 4: Braced pitched frame, wind-like: eave FX point loads + asymmetric rafter FY UDLs

  Frame geometry : H=3.00 m, B_SPAN=5.00 m, RIDGE_H=1.00 m
  Ridge at       : (2.50, 4.00, 0)
  Rafter length  : 2.693 m   slope = 21.8 deg
  Brace  (A->C)  : 5.831 m   angle = 31.0 deg
  Base condition : Fixed
  Sub-elements   : 8 per member

  Column section : IPE200
    A  = 28.48 cm^2
    Iy = 1943 cm^4  (strong axis)
    Iz = 142 cm^4  (weak axis)
    J  = 7.00 cm^4
    Asy= 17.00 cm^2  (shear area, weak-axis bending)
    Asz= 14.02 cm^2  (shear area, strong-axis bending)
  Rafter section : IPE200
    A  = 28.48 cm^2
    Iy = 1943 cm^4  (strong axis)
    Iz = 142 cm^4  (weak axis)
    J  = 7.00 cm^4
    Asy= 17.00 cm^2  (shear area, weak-axis bending)
    Asz= 14.02 cm^2  (shear area, strong-axis bending)
  Brace section : IPE200
    A  = 28.48 cm^2
    Iy = 1943 cm^4  (strong axis)
    Iz = 142 cm^4  (weak axis)
    J  = 7.00 cm^4
    Asy= 17.00 cm^2  (shear area, weak-axis bending)
    Asz= 14.02 cm^2  (shear area, strong-axis bending)
  Beam formulation : Timoshenko (shear deformation included)

  Applied loads
  Raft1 (B->R) : w = -10.00 kN/m  (distributed FY)
  Raft2 (R->C) : w = +10.00 kN/m  (distributed FY)
  Node B       : FX = +10.0 kN  (point load)
  Node C       : FX = -10.0 kN  (point load)
  Lateral bracing : Yes — 2D constraint (DZ, RX, RY restrained at all nodes)

  Running buckling analysis ...

----------------------------------------------------------------
  Mode    Lambda_cr   Notes
----------------------------------------------------------------
     1     417.1316  <-- critical (lowest)
     2     665.1076
     3    1076.9499
     4    1455.9565
     5    1919.9956
----------------------------------------------------------------

----------------------------------------------------------------
  Effective buckling lengths  (H = 3.00 m, L_raft = 2.693 m, L_brace = 5.831 m)
----------------------------------------------------------------
  Member    Mode       Plane   L_cr[m]   L_cr/L   N_cr[kN]
----------------------------------------------------------------
  Col1         1  y (strong)     2.548    0.849     6200.8
  Col1         1    z (weak)     0.689    0.230     6200.8
  Col1         2  y (strong)     2.018    0.673     9887.0
  Col1         2    z (weak)     0.546    0.182     9887.0
  Col1         3  y (strong)     1.586    0.529    16009.2
  Col1         3    z (weak)     0.429    0.143    16009.2
  Col1         4  y (strong)     1.364    0.455    21643.3
  Col1         4    z (weak)     0.369    0.123    21643.3
  Col1         5  y (strong)     1.188    0.396    28541.4
  Col1         5    z (weak)     0.321    0.107    28541.4

  Col2         1  y (strong)     2.705    0.902     5502.7
  Col2         1    z (weak)     0.731    0.244     5502.7
  Col2         2  y (strong)     2.142    0.714     8773.9
  Col2         2    z (weak)     0.579    0.193     8773.9
  Col2         3  y (strong)     1.684    0.561    14206.8
  Col2         3    z (weak)     0.455    0.152    14206.8
  Col2         4  y (strong)     1.448    0.483    19206.6
  Col2         4    z (weak)     0.391    0.130    19206.6
  Col2         5  y (strong)     1.261    0.420    25328.1
  Col2         5    z (weak)     0.341    0.114    25328.1

  Raft1        1  y (strong)     2.629    0.976     5825.4
  Raft1        1    z (weak)     0.711    0.264     5825.4
  Raft1        2  y (strong)     2.082    0.773     9288.5
  Raft1        2    z (weak)     0.563    0.209     9288.5
  Raft1        3  y (strong)     1.636    0.608    15040.0
  Raft1        3    z (weak)     0.442    0.164    15040.0
  Raft1        4  y (strong)     1.407    0.523    20333.0
  Raft1        4    z (weak)     0.380    0.141    20333.0
  Raft1        5  y (strong)     1.226    0.455    26813.4
  Raft1        5    z (weak)     0.331    0.123    26813.4

  Raft2        1  y (strong)     2.733    1.015     5390.9
  Raft2        1    z (weak)     0.739    0.274     5390.9
  Raft2        2  y (strong)     2.164    0.804     8595.6
  Raft2        2    z (weak)     0.585    0.217     8595.6
  Raft2        3  y (strong)     1.701    0.632    13918.1
  Raft2        3    z (weak)     0.460    0.171    13918.1
  Raft2        4  y (strong)     1.463    0.543    18816.3
  Raft2        4    z (weak)     0.395    0.147    18816.3
  Raft2        5  y (strong)     1.274    0.473    24813.4
  Raft2        5    z (weak)     0.344    0.128    24813.4

  Brc          1  y (strong)     5.448    0.934     1356.9
  Brc          1    z (weak)     1.473    0.253     1356.9
  Brc          2  y (strong)     4.314    0.740     2163.5
  Brc          2    z (weak)     1.166    0.200     2163.5
  Brc          3  y (strong)     3.391    0.581     3503.2
  Brc          3    z (weak)     0.917    0.157     3503.2
  Brc          4  y (strong)     2.916    0.500     4736.0
  Brc          4    z (weak)     0.788    0.135     4736.0
  Brc          5  y (strong)     2.539    0.435     6245.5
  Brc          5    z (weak)     0.686    0.118     6245.5

----------------------------------------------------------------
  Euler reference  (pinned-pinned, N_cr = pi^2*EI/L^2)
----------------------------------------------------------------
  Column  (L=H=3.00 m):
    Strong axis (Iy):  N_cr = 4474.5 kN   L_eff/H = 1.00
    Weak axis  (Iz):  N_cr = 327.0 kN   L_eff/H = 1.00
  Rafter  (L=L_raft=2.693 m):
    Strong axis (Iy):  N_cr = 5554.6 kN   L_eff/L = 1.00
    Weak axis  (Iz):  N_cr = 405.9 kN   L_eff/L = 1.00
  Brace   (L=L_brace=5.831 m):
    Strong axis (Iy):  N_cr = 1184.4 kN   L_eff/L = 1.00
    Weak axis  (Iz):  N_cr = 86.6 kN   L_eff/L = 1.00

----------------------------------------------------------------
  Column axial forces from static pre-solve
----------------------------------------------------------------
  Member       N_Ed [kN]  Sign
----------------------------------------------------------------
  Col1             14.87  compression
  Col2            -13.19  tension
----------------------------------------------------------------

----------------------------------------------------------------
  Rafter & brace internal forces  (compare with PolyFrame / FEM Design)
----------------------------------------------------------------
  Member       N_Ed [kN]   M_max [kNm]  Sign_N
----------------------------------------------------------------
  Raft1            13.97          6.30  compr.
  Raft2            12.92          8.41  compr.
  Brc              -3.25          0.00  tension
----------------------------------------------------------------

  Done.  Compare Lambda_cr and L_cr/L with PolyFrame / FEM Design output.

----------------------------------------------------------------
  FEM Design comparison
----------------------------------------------------------------
  Mode      PyNite  FEM Design      Diff
----------------------------------------------------------------
     1     417.132     414.113     +0.7%
     2     665.108     659.500     +0.9%
     3    1076.950    1064.552     +1.2%
     4    1455.957    1455.810     +0.0%
     5    1919.996         n/a
----------------------------------------------------------------
```

### buckling_comparison.py

Two-bay single-storey portal frame: column and beam buckling scenarios. Status: **passed** (1.0s).

```text
================================================================
  PyNite -- Linear Buckling Analysis
================================================================

  Scenario 2: 2-bay, fixed bases, midspan point loads on beams (global FY)

  Frame geometry : H=3.00 m,  B1=5.00 m, B2=5.00 m
  Base condition : Fixed
  Sub-elements   : 8 per member

  Column section : IPE200
    A  = 28.48 cm^2
    Iy = 1943 cm^4  (strong axis)
    Iz = 142 cm^4  (weak axis)
    J  = 7.00 cm^4
    Asy= 17.00 cm^2  (shear area, weak-axis bending)
    Asz= 14.02 cm^2  (shear area, strong-axis bending)
  Beam section : IPE200
    A  = 28.48 cm^2
    Iy = 1943 cm^4  (strong axis)
    Iz = 142 cm^4  (weak axis)
    J  = 7.00 cm^4
    Asy= 17.00 cm^2  (shear area, weak-axis bending)
    Asz= 14.02 cm^2  (shear area, strong-axis bending)
  Beam formulation : Timoshenko (shear deformation included)

  Applied loads
  Beam1 (B->C) : P = -10.0 kN @ midspan (global FY)
  Beam2 (C->D) : P = -10.0 kN @ midspan (global FY)
  Lateral brace at beam nodes : No

  Running buckling analysis ...

----------------------------------------------------------------
  Mode    Lambda_cr   Notes
----------------------------------------------------------------
     1     428.9898  <-- critical (lowest)
     2     884.6188
     3    1687.5880
     4    1828.8785
     5    2104.9027
----------------------------------------------------------------

----------------------------------------------------------------
  Effective buckling lengths  (column height H = 3.00 m)
----------------------------------------------------------------
  Member    Mode       Plane   L_cr[m]   L_cr/H   N_cr[kN]
----------------------------------------------------------------
  Col1         1  y (strong)     4.675    1.558     1842.6
  Col1         1    z (weak)     1.264    0.421     1842.6
  Col1         2  y (strong)     3.256    1.085     3799.7
  Col1         2    z (weak)     0.880    0.293     3799.7
  Col1         3  y (strong)     2.357    0.786     7248.7
  Col1         3    z (weak)     0.637    0.212     7248.7
  Col1         4  y (strong)     2.264    0.755     7855.6
  Col1         4    z (weak)     0.612    0.204     7855.6
  Col1         5  y (strong)     2.110    0.703     9041.2
  Col1         5    z (weak)     0.571    0.190     9041.2

  Col2         1  y (strong)     2.868    0.956     4894.5
  Col2         1    z (weak)     0.775    0.258     4894.5
  Col2         2  y (strong)     1.998    0.666    10092.9
  Col2         2    z (weak)     0.540    0.180    10092.9
  Col2         3  y (strong)     1.446    0.482    19254.3
  Col2         3    z (weak)     0.391    0.130    19254.3
  Col2         4  y (strong)     1.389    0.463    20866.3
  Col2         4    z (weak)     0.376    0.125    20866.3
  Col2         5  y (strong)     1.295    0.432    24015.6
  Col2         5    z (weak)     0.350    0.117    24015.6

  Col3         1  y (strong)     4.675    1.558     1842.6
  Col3         1    z (weak)     1.264    0.421     1842.6
  Col3         2  y (strong)     3.256    1.085     3799.7
  Col3         2    z (weak)     0.880    0.293     3799.7
  Col3         3  y (strong)     2.357    0.786     7248.7
  Col3         3    z (weak)     0.637    0.212     7248.7
  Col3         4  y (strong)     2.264    0.755     7855.6
  Col3         4    z (weak)     0.612    0.204     7855.6
  Col3         5  y (strong)     2.110    0.703     9041.2
  Col3         5    z (weak)     0.571    0.190     9041.2

----------------------------------------------------------------
  Single-column Euler reference  (pinned-pinned, N_cr = pi^2*EI/H^2)
----------------------------------------------------------------
  Strong axis (Iy):  N_cr = 4474.5 kN   L_eff/H = 1.00
  Weak axis  (Iz):  N_cr = 327.0 kN   L_eff/H = 1.00

----------------------------------------------------------------
  Column axial forces from static pre-solve
----------------------------------------------------------------
  Member       N_Ed [kN]  Sign
----------------------------------------------------------------
  Col1              4.30  compression
  Col2             11.41  compression
  Col3              4.30  compression
----------------------------------------------------------------

----------------------------------------------------------------
  Beam internal forces  (compare with PolyFrame N_Ed, M_Ed)
----------------------------------------------------------------
  Member       N_Ed [kN]   M_max [kNm]  Sign_N
----------------------------------------------------------------
  Beam1             1.91          7.40  compr.
  Beam2             1.91          7.40  compr.
----------------------------------------------------------------

  Done.  Compare Lambda_cr and L_cr/H with PolyFrame / FEM Design output.

----------------------------------------------------------------
  FEM Design comparison
----------------------------------------------------------------
  Mode      PyNite  FEM Design      Diff
----------------------------------------------------------------
     1     428.990     426.220     +0.6%
     2     884.619     875.321     +1.1%
     3    1687.588    1691.256     -0.2%
     4    1828.879    1818.937     +0.5%
     5    2104.903         n/a
----------------------------------------------------------------
```

### hinged_frame_comparison.py

Single-bay portal frame with a moment release at the top-left corner. Status: **passed** (0.9s).

```text
================================================================
  PyNite -- Hinged Portal Frame Buckling Analysis
================================================================

  Scenario 1: Vertical UDL on beam: 10 kN/m downward

  Frame geometry : H = 3.00 m, B_SPAN = 4.00 m
  Base condition : Fixed
  Hinge          : top-left corner (beam i-end at node B)
  Sub-elements   : 8 per member

  Section (all members): IPE300
    A  = 53.80 cm^2
    Iy = 8356 cm^4  (strong axis)
    Iz = 604 cm^4  (weak axis)
    J  = 20.10 cm^4
    Asy= 32.10 cm^2  (shear area, weak-axis bending)
    Asz= 25.67 cm^2  (shear area, strong-axis bending)
  Beam formulation : Timoshenko

  Applied loads
  Beam1 (B->C) : w = -10.00 kN/m  (distributed FY)
  Lateral bracing : Yes -- 2D constraint (DZ, RX, RY restrained at all non-base nodes)

  Running buckling analysis ...

----------------------------------------------------------------
  Mode    Lambda_cr   Notes
----------------------------------------------------------------
     1     361.5503  <-- critical (lowest)
     2    1755.2618
     3    1908.2979
     4    3393.9275
     5    4090.2281
----------------------------------------------------------------

----------------------------------------------------------------
  Effective buckling lengths  (H = 3.00 m, B_SPAN = 4.00 m)
----------------------------------------------------------------
  Member    Mode       Plane   L_cr[m]   L_cr/L   N_cr[kN]
----------------------------------------------------------------
  Col1         1  y (strong)     5.164    1.721     6493.9
  Col1         1    z (weak)     1.388    0.463     6493.9
  Col1         2  y (strong)     2.344    0.781    31526.5
  Col1         2    z (weak)     0.630    0.210    31526.5
  Col1         3  y (strong)     2.248    0.749    34275.3
  Col1         3    z (weak)     0.604    0.201    34275.3
  Col1         4  y (strong)     1.686    0.562    60958.9
  Col1         4    z (weak)     0.453    0.151    60958.9
  Col1         5  y (strong)     1.535    0.512    73465.3
  Col1         5    z (weak)     0.413    0.138    73465.3

  Col2         1  y (strong)     4.662    1.554     7968.1
  Col2         1    z (weak)     1.253    0.418     7968.1
  Col2         2  y (strong)     2.116    0.705    38683.9
  Col2         2    z (weak)     0.569    0.190    38683.9
  Col2         3  y (strong)     2.029    0.676    42056.7
  Col2         3    z (weak)     0.546    0.182    42056.7
  Col2         4  y (strong)     1.522    0.507    74798.2
  Col2         4    z (weak)     0.409    0.136    74798.2
  Col2         5  y (strong)     1.386    0.462    90143.9
  Col2         5    z (weak)     0.373    0.124    90143.9

  Beam1        1  y (strong)    15.569    3.892      714.5
  Beam1        1    z (weak)     4.186    1.046      714.5
  Beam1        2  y (strong)     7.066    1.766     3468.9
  Beam1        2    z (weak)     1.900    0.475     3468.9
  Beam1        3  y (strong)     6.777    1.694     3771.3
  Beam1        3    z (weak)     1.822    0.455     3771.3
  Beam1        4  y (strong)     5.081    1.270     6707.3
  Beam1        4    z (weak)     1.366    0.342     6707.3
  Beam1        5  y (strong)     4.629    1.157     8083.4
  Beam1        5    z (weak)     1.244    0.311     8083.4

----------------------------------------------------------------
  Euler reference  (pinned-pinned, N_cr = pi^2*EI/L^2)
----------------------------------------------------------------
  Column  (L = H = 3.00 m):
    Strong axis (Iy):  N_cr = 19243.1 kN
    Weak axis  (Iz):  N_cr = 1391.0 kN
  Beam    (L = B_SPAN = 4.00 m):
    Strong axis (Iy):  N_cr = 10824.2 kN
    Weak axis  (Iz):  N_cr = 782.4 kN

----------------------------------------------------------------
  Member axial forces from static pre-solve
----------------------------------------------------------------
  Member       N_Ed [kN]  Sign
----------------------------------------------------------------
  Col1             17.96  compression
  Col2             22.04  compression
  Beam1             1.98  compression
----------------------------------------------------------------

----------------------------------------------------------------
  Member internal forces  (N_Ed at i-end, M_max along member)
----------------------------------------------------------------
  Member       N_Ed [kN]   M_max [kNm]  Sign_N
----------------------------------------------------------------
  Col1             17.96          5.93  compr.
  Col2             22.04          8.16  compr.
  Beam1             1.98         16.13  compr.
----------------------------------------------------------------

  Done.  Compare Lambda_cr and L_cr/L with reference output.

----------------------------------------------------------------
  FEM Design comparison
----------------------------------------------------------------
  Mode       PyNite   FEM Design   Diff [%]   Status
----------------------------------------------------------------
     1      361.550      358.320      +0.90   OK
     2     1755.262     1721.181      +1.98   OK
     3     1908.298     1874.093      +1.83   OK
     4     3393.928     3408.881      -0.44   OK
----------------------------------------------------------------
  RESULT: ALL MODES WITHIN 5% TOLERANCE
```

### pitched_frame_comparison.py

Single-bay pitched-roof (gable) portal frame. Status: **passed** (0.8s).

```text
================================================================
  PyNite -- Pitched-Roof Portal Frame Buckling Analysis
================================================================

  Scenario 4: Pitched frame, wind-like: column FX UDL + asymmetric rafter FY UDLs

  Frame geometry : H=3.00 m, B_SPAN=5.00 m, RIDGE_H=1.00 m
  Ridge at       : (2.50, 4.00, 0)
  Rafter length  : 2.693 m   slope = 21.8 deg
  Base condition : Fixed
  Sub-elements   : 8 per member

  Column section : IPE200
    A  = 28.48 cm^2
    Iy = 1943 cm^4  (strong axis)
    Iz = 142 cm^4  (weak axis)
    J  = 7.00 cm^4
    Asy= 17.00 cm^2  (shear area, weak-axis bending)
    Asz= 14.02 cm^2  (shear area, strong-axis bending)
  Rafter section : IPE200
    A  = 28.48 cm^2
    Iy = 1943 cm^4  (strong axis)
    Iz = 142 cm^4  (weak axis)
    J  = 7.00 cm^4
    Asy= 17.00 cm^2  (shear area, weak-axis bending)
    Asz= 14.02 cm^2  (shear area, strong-axis bending)
  Beam formulation : Timoshenko (shear deformation included)

  Applied loads
  Raft1 (B->R) : w = -10.00 kN/m  (distributed FY)
  Raft2 (R->C) : w = +10.00 kN/m  (distributed FY)
  Col1          : w = +10.00 kN/m  (distributed FX)
  Col2          : w = +10.00 kN/m  (distributed FX)
  Lateral bracing : Yes — 2D constraint (DZ, RX, RY restrained at all nodes)

  Running buckling analysis ...

----------------------------------------------------------------
  Mode    Lambda_cr   Notes
----------------------------------------------------------------
     1     986.7696  <-- critical (lowest)
     2    1187.9311
     3    2955.1747
     4    3673.1517
     5    4746.9796
----------------------------------------------------------------

----------------------------------------------------------------
  Effective buckling lengths  (H = 3.00 m, L_raft = 2.693 m)
----------------------------------------------------------------
  Member    Mode       Plane   L_cr[m]   L_cr/L   N_cr[kN]
----------------------------------------------------------------
  Col1         1  y (strong)     2.059    0.686     9497.8
  Col1         1    z (weak)     0.557    0.186     9497.8
  Col1         2  y (strong)     1.877    0.626    11434.0
  Col1         2    z (weak)     0.507    0.169    11434.0
  Col1         3  y (strong)     1.190    0.397    28444.0
  Col1         3    z (weak)     0.322    0.107    28444.0
  Col1         4  y (strong)     1.067    0.356    35354.6
  Col1         4    z (weak)     0.289    0.096    35354.6
  Col1         5  y (strong)     0.939    0.313    45690.4
  Col1         5    z (weak)     0.254    0.085    45690.4

  Col2         1  y (strong)     2.059    0.686     9497.8
  Col2         1    z (weak)     0.557    0.186     9497.8
  Col2         2  y (strong)     1.877    0.626    11434.0
  Col2         2    z (weak)     0.507    0.169    11434.0
  Col2         3  y (strong)     1.190    0.397    28444.0
  Col2         3    z (weak)     0.322    0.107    28444.0
  Col2         4  y (strong)     1.067    0.356    35354.6
  Col2         4    z (weak)     0.289    0.096    35354.6
  Col2         5  y (strong)     0.939    0.313    45690.4
  Col2         5    z (weak)     0.254    0.085    45690.4

  Raft1        1  y (strong)     3.379    1.255     3527.4
  Raft1        1    z (weak)     0.913    0.339     3527.4
  Raft1        2  y (strong)     3.080    1.144     4246.5
  Raft1        2    z (weak)     0.833    0.309     4246.5
  Raft1        3  y (strong)     1.952    0.725    10563.8
  Raft1        3    z (weak)     0.528    0.196    10563.8
  Raft1        4  y (strong)     1.751    0.650    13130.4
  Raft1        4    z (weak)     0.473    0.176    13130.4
  Raft1        5  y (strong)     1.541    0.572    16969.0
  Raft1        5    z (weak)     0.416    0.155    16969.0

  Raft2        1  y (strong)     2.520    0.936     6340.3
  Raft2        1    z (weak)     0.681    0.253     6340.3
  Raft2        2  y (strong)     2.297    0.853     7632.8
  Raft2        2    z (weak)     0.621    0.231     7632.8
  Raft2        3  y (strong)     1.456    0.541    18987.9
  Raft2        3    z (weak)     0.394    0.146    18987.9
  Raft2        4  y (strong)     1.306    0.485    23601.1
  Raft2        4    z (weak)     0.353    0.131    23601.1
  Raft2        5  y (strong)     1.149    0.427    30500.8
  Raft2        5    z (weak)     0.311    0.115    30500.8

----------------------------------------------------------------
  Euler reference  (pinned-pinned, N_cr = pi^2*EI/L^2)
----------------------------------------------------------------
  Column  (L=H=3.00 m):
    Strong axis (Iy):  N_cr = 4474.5 kN   L_eff/H = 1.00
    Weak axis  (Iz):  N_cr = 327.0 kN   L_eff/H = 1.00
  Rafter  (L=L_raft=2.693 m):
    Strong axis (Iy):  N_cr = 5554.6 kN   L_eff/L = 1.00
    Weak axis  (Iz):  N_cr = 405.9 kN   L_eff/L = 1.00

----------------------------------------------------------------
  Column axial forces from static pre-solve
----------------------------------------------------------------
  Member       N_Ed [kN]  Sign
----------------------------------------------------------------
  Col1              9.63  compression
  Col2             -9.63  tension
----------------------------------------------------------------

----------------------------------------------------------------
  Rafter internal forces  (compare with PolyFrame / FEM Design)
----------------------------------------------------------------
  Member       N_Ed [kN]   M_max [kNm]  Sign_N
----------------------------------------------------------------
  Raft1             3.57         13.89  compr.
  Raft2             6.43         13.89  compr.
----------------------------------------------------------------

  Done.  Compare Lambda_cr and L_cr/L with PolyFrame / FEM Design output.

----------------------------------------------------------------
  FEM Design comparison
----------------------------------------------------------------
  Mode      PyNite  FEM Design      Diff
----------------------------------------------------------------
     1     986.770     980.953     +0.6%
     2    1187.931    1172.103     +1.4%
     3    2955.175         n/a
     4    3673.152         n/a
     5    4746.980         n/a
----------------------------------------------------------------
```
