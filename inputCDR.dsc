! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
! 2d Convection-Diffusion-Reaction simulation                 !
! Input data file                                             !
!                                   MAOG       Bcn, Dic. 2021 !
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

# > > > > > > > Geometry
ElemType = Quadrilateral !Triangle, Quadrilateral
ProbType = trans         !Problem type transient or static
DimPr    = 2             !Dimension del problema
nelem    = 9           !Number of lnods
nnodes   = 16           !Total number of nodal points
nne      = 4             !Number of nodes in the element
ndofn    = 2             !Degrees of freedom
totGp    = 4             !1,4,9 for Q, 1,3,4,6 for P
maxband  = 78             !Maximo ancho de banda

# > > > > > > > Time Discretization
time_ini = 1.0             !Starting time simulation
time_fin = 10.0             !Ending time simulation
max_time = 10              !Max time of simulation
u0_cond  = 0.0             !Value of initial condition

# > > > > > > > Stabilization
kstab    = 3             !Type of stabilization: 1(SUPG), 2(GLS), 3/5(SGS/TG), 4(CG)
ktaum    = 1             !Type of tau matrix: 0, 1, 2 
patau    = 1.0           !Parameter to obtain tau
hnatu    = 2.0           !Reference element length

# > > > > > > > Physical Properties
#DIFMA_11                        !Diffusion tensor
   1000    , 0.0E-0  , 0.0E-0
0.0E-0  ,   1000     , 0.0E-0
0.0E-0  , 0.0E-0  , 0
#DIFMA_12
0 ,  0.0e-0 ,  0.0e-0
0.0e-0  , 0  , 0.0e-0
0.0e-0  , 0.0e-0  , 0
#DIFMA_22
   1000   , 0.0E-0 ,  0.0E-0
0.0E-0 ,    1000  ,  0.0E-0
0.0E-0 , 0.0E-0 , 0
#COMAT_1                        !Convection tensor
100 , 0 , 0.0e-0
0.0e-0 , 0.0e-0 , 0.0e-0
0.0e-0 , 0.0e-0 , 0.0e-0
#COMAT_2
0 , 0.0e-0 , 0
0.0e-0 , 0.0e-0 , 0.0e-0
0.0e-0 , 0.0e-0 , 0.0e-0
#REAMA                          !Reaction tensor
0.0e-0 , 0.0e-0 , 0
0.0e-0 , 0 , 0.0e-0
0.0e-0 , 0.0e-0 , 0
#FORCE                          !Force tensor
0 , 0.5 , 0.0e+0

# > > > > > > > Mesh -> Element's Nodes -> Coordinate's Nodes
1     1     2     6     5
2     2     3     7     6
3     3     4     8     7
4     5     6    10     9
5     6     7    11    10
6     7     8    12    11
7     9    10    14    13
8    10    11    15    14
9    11    12    16    15
 1    0.0000    0.0000
 2    0.3333    0.0000
 3    0.6667    0.0000
 4    1.0000    0.0000
 5    0.0000    0.3333
 6    0.3333    0.3333
 7    0.6667    0.3333
 8    1.0000    0.3333
 9    0.0000    0.6667
10    0.3333    0.6667
11    0.6667    0.6667
12    1.0000    0.6667
13    0.0000    1.0000
14    0.3333    1.0000
15    0.6667    1.0000
16    1.0000    1.0000
