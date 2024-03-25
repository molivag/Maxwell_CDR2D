! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ∂tE + λ∇x∇xE + ∇P + β∇∇⋅E = -∂J/∂t  in  Ω 
! 2d Convection-Diffusion-Reaction simulation             !                ∇⋅E - ɣΔP  = 0       in  Ω 
! Input data file                                         !                       nxE = 0       on  ∂Ω
!                                  MAOG    Bcn, Dic. 2021 ! with:
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       λ, β, ɣ coefficients.    
 
# > > > > > > > Model Parameters
ProbType = TIME           !Problem type TIME=transient, other=static
DimPr    = 2              !Dimension del problema
ndofn    = 3              !Degrees of freedom
totGp    = 3              !1,4,9 for Q, 1,3,7 for P
exacSol  = 3              !0=None;   1=SinglrSol  ; 2=FullSpace; 3=Algebraic; 4=Double Line
srcRHS   = 0              !0=scalar; 1=SingularSol; 2=Maxwell_Polynom; 3=Lapalace_Polynom
BCsProb  = 2              !1=Ldomain; 2=Maxwell; 3=MaxwellPoly; 4=Lapalace_Poly; 5=Cavity-Driven Flow 
postpro  = 2              !Execution of post-processing routine 1=yes, 2=no 
sigma    = 1.0            !Conductivity of the medium
operator = MAXW           !PDE being solved LAPL=Laplacian; MAXW=curl-curl

#***************Geometry
#---------!File .msh that contains the mesh
meshfile = TEM_VMD_ground_air.msh  
view     = xy             !The 2D view x-y (distance) or x-z (depth)
nne      = 3              !Nodes per element Q:4-9; P:3-6
i_exp    = 0              !Exponent of characteristic mesh size 3,4,5 or 6 2^(-i)
hnatu    = 1.0            !Reference element length
refiType = NO             !NONE; PS=Powell-Sabin; CC=Criss-Cross

# > > > > > > > Stabilization
kstab    = 6              !Stabilization: 0(NONE), 1(SUPG), 2(GLS), 3/5(SGS/TG), 4(CG), 6(MVAF)
ktaum    = 1              !Tau matrix: 0, 1, 2 
patau    = 2.0            !Parameter to obtain tau
n_val    = 0.0            !n parameter in exact solution, for simul=1
helem    = 1.0            !Characteristic mesh size (maximum element size among the mesh)
Cu       = 500.0           !Algorithmic constant
ell      = 1000.0          !Constante de longitud  
1/mu=λ   = 795774.71545   !Reluctivity of the medium 1/µ0=795774.71545 [T•m•A^-1]

# > > > > > > > Fourier Transform     
TwoHalf  = N              !If it is dealing with a 2.5D modeling (Y) or not (N)
ky_min   = 1.0e-3         !Smallest wave number
ky_max   = 1.0e0          !Grater wave number
tot_ky   = 4             !Total wave numbers
splits   = N              !If the problem is splited then run as kind of parallel
y_iFT    = 0.0            !Position at y-coordinate where the IFT will be computed.
#------------!File to plot shape spectrum
specFile = xxxxxxtransformed_Voltage.dat
#------------!File name after inverse Fourier Transform 
File_iFT = Final 3D Direct Current Impresion_

# > > > > > > > Time Discretization
theta    = 2              !BDF1=2 ;CN=3; BDF2=4
time_ini = 1.42000e-7     !Starting time simulation (simulation always starts at 0?)
time_fin = 6.12800e-5     !Total time simulated in [s]  --> 1800 microseconds
t_steps  = 10             !Number of time steps
Src_ON   = 2              !Time at Source is turned ON

# > > > > > > > Name outPut Files
testID   = TEM_2point5D_test_04
postpro  = TEM_2point5D_test_04
Error    = xxxxxxxx@xxx
Cordina  = xxxxxxxxxxxx
Conecti  = xxxxxxxxxxxx
Profile  = xxxxxxxxxxxx

# > > > > > > > Physical Properties
#DIFMA_xx                  !Diffusion tensor		--     --
1.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0	| Re_ex |
0.0 , 1.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0	| Re_ey |
0.0 , 0.0  , 1.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0	| Re_ez |
0.0 , 0.0  , 0.0 , 1.0  , 0.0  , 0.0 , 0.0  , 0.0	| Re_p  |
0.0 , 0.0  , 0.0 , 0.0  , 1.0  , 0.0 , 0.0  , 0.0	| Im_ex |
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 1.0 , 0.0  , 0.0	| Im_ey |
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 1.0  , 0.0	| Im_ez |
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 1.0	| Im_p  |
#DIFMA_xz						--     --
0.0 , -1.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
1.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
#DIFMA_zx
0.0 , 1.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
-1.0 , 0.0 , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
#DIFMA_zz
1.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 1.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 1.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
0.0 , 0.0  , 0.0 , 0.0  , 0.0  , 0.0 , 0.0  , 0.0
#CONMAT_x                   !Convection tensor
0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
#CONMAT_z
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
#REAMA                     !Reaction tensor
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
#FORCE                     !Force tensor
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

# > > > > > > > Source Configuration
#Icurr
1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

#Nodal source Location
nodalSrc = 2            !Number of nodes will contain the source
19
13

#Time waveform
signal   = 1            !Signal in time: 1=step-on; 2=step-off; 3=triangular

# > > > > > > > Receiver Locations
nodalRec = 1            !Number of nodes as a receiver
0.0  , 0.0

