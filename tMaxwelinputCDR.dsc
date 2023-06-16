! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ∂tE + λ∇x∇xE + ∇P + β∇∇⋅E = -∂J/∂t  in  Ω 
! 2d Convection-Diffusion-Reaction simulation                 !                ∇⋅E - ɣΔP  = 0       in  Ω 
! Input data file                                             !                       nxE = 0       on  ∂Ω
!                                   MAOG       Bcn, Dic. 2021 ! with:
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       λ, β, ɣ coefficients.    
 
# > > > > > > > Model Parameters
ElemType = QUAD             !Could be QUAD or TRIA
ProbType = TIrE             !Problem type TIME=transient, other=static
DimPr    = 2                !Dimension del problema
ndofn    = 3                !Degrees of freedom
totGp    = 4                !1,4,9 for Q, 1,3,7 for P
simul    = 2                !1=LdomT2; 2=FullSpace; 3=PolyMaxwell; 4=PolyStokes; 5=Cavity-Driven Flow 
pospro   = 2                !Execution of post-processing routine 1=yes, 2=no

# > > > > > > > Geometry
meshfile = FixEMxyz.msh     !File .msh that contains the mesh
nnodes   = 1681             !Total nodal points
nelem    = 1600             !Total elements
nne      = 4                !Nodes per element Q:4-9; P:3-6
i_exp    = 0                !Exponent of characteristic mesh size 3,4,5 or 6 2^(-i)
hnatu    = 2.0              !Reference element length
refiType = NO               !NONE; PS=Powell-Sabin; CB=Crossed-Box

# > > > > > > > Stabilization
kstab    = 6                !Stabilization: 0(NONE), 1(SUPG), 2(GLS), 3/5(SGS/TG), 4(CG), 6(MVAF)
ktaum    = 0                !Tau matrix: 0, 1, 2 
patau    = 0.0              !Parameter to obtain tau
n_val    = 0.0              !n parameter in exact solution, for simul=1
helem    = 0.05             !Characteristic mesh size (maximum element size among the mesh)
Cu       = 100.0              !Algorithmic constant
ell      = 2.0              !Constante de longitud  
1/mu=λ   = 795774.7155      !Reluctivity of the medium 1/µ0=795774.71545 [T•m•A^-1]

# > > > > > > > Time Discretization
theta    = 2                !BDF1=2 ;CN=3; BDF2=4
time_ini = 0.00000          !Starting time simulation (simulation always starts at 0?)
time_fin = 150.000          !0.00001          !Time simulated (total time simulated [s])
max_time = 90              !Time steps
u0cond   = 0.0              !Value of initial condition (could be defined here or codeing at mod_timeInt.f90)

# > > > > > > > Name outPut Files
testID   = 3.10RefSol_MVAF   !data file with input parameters in each iteration Res/results
Postpro  = 3.10RefSolMVA     
Error    = erRefSolMVAF
Cordina  = coord0.2MVAF
Conecti  = conec0.2MVAF
Profile  = 0.0RefSolMVA

# > > > > > > > Physical Properties
#DIFMA_xx                   !Diffusion tensor
1.0  , 0.0  , 0.0
0.0  , 1.0 , 0.0
0.0  , 0.0  , 1.0
#DIFMA_xy
0.0  , 1.0  , 0.0
-1.0 , 0.0  , 0.0 
0.0  , 0.0  , 0.0 
#DIFMA_yx
0.0  , -1.0 , 0.0 
1.0 , 0.0  , 0.0
0.0  , 0.0  , 0.0
#DIFMA_yy
1.0 , 0.0  , 0.0
0.0  , 1.0  , 0.0
0.0  , 0.0  , 1.0
#COMAT_x                    !Convection tensor
0.0  , 0.0 , 1.0
0.0  , 0.0 , 0.0
1.0 , 0.0 , 0.0
#COMAT_y
0.0 , 0.0  , 0.0
0.0 , 0.0  , 1.0
0.0 , 1.0 , 0.0
#REAMA                      !Reaction tensor
0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0
#FORCE                      !Force tensor
0.0 , 0.0 , 0.0 

# > > > > > > > Source Configuration
#Icurr
0.0 , 1.0 , 0.0

#Nodal source Location
nodalSrc = 1                !Number of nodes will contain the source
659

#Time waveform
signal   = 1                !Signal in time: 1=step-on; 2=step-off; 3=triangular

# > > > > > > > Nodal receiver Locations
nodalRec = 1                !Number of nodes as a receiver
659                         !Directo en la fuente
