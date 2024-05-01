! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ∂tE + λ∇x∇xE + ∇P + β∇∇⋅E = -∂J/∂t  in  Ω 
! 2d Convection-Diffusion-Reaction simulation             !                ∇⋅E - ɣΔP  = 0       in  Ω 
! Input data file                                         !                       nxE = 0       on  ∂Ω
!                                  MAOG    Bcn, Dic. 2021 ! with:
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       λ, β, ɣ coefficients.    
 
# > > > > > > > Model Parameters
ProbType = TIME           !Problem type TIME=transient, other=static
totGp    = 3              !1,4,9 for Q, 1,3,7 for P
exacSol  = 3              !0=None;   1=SinglrSol  ; 2=FullSpace; 3=Algebraic; 4=Double Line
srcRHS   = 0              !0=scalar; 1=SingularSol; 2=Maxwell_Polynom; 3=Lapalace_Polynom
BCsProb  = 0              !1=Ldomain; 2=Maxwell; 3=MaxwellPoly; 4=Lapalace_Poly; 5=Cavity-Driven Flow 
postpro  = 2              !Execution of post-processing routine 1=yes, 2=no 
sigma    = 1.0            !Conductivity of the medium

$********************************************************************************
PHYSICAL_PROBLEM
$----------------------------------------------------------------------------------
  Maxwell_In_Non_Convex_Domain
  Cavity_Driven_Flow
  Direct_Current_Electrical_Resistivity_in_2.5-D
  Electric_Field_Excited_By_A_Double_Line_Source
 >Horizontal_Electric_Dipole_in_3-D_A_Wrong_capture_of_solution
  Transient_Electromagnetic_in_2.5-D
$----------------------------------------------------------------------------------
END_PHYSICAL_PROBLEM
$*********************************************************************************

#***************Geometry
#---------!File .msh that contains the mesh
meshfile = TEM_HED_ground_air.msh  
view     = xz             !The 2D view x-y (distance) or x-z (depth)
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
tot_ky   = 3              !Total wave numbers
splits   = N              !If the problem is splited then run as kind of parallel
y_iFT    = 0.0            !Position at y-coordinate where the IFT will be computed.
#------------!File to plot shape spectrum
specFile = TEM_3D2D_test05_transformed field.dat
#------------!File name after inverse Fourier Transform 
File_iFT = 3D TEM response

# > > > > > > > Transient Parameters
theta    = 2              !BDF1=2 ;CN=3; BDF2=4
time_ini = 1.42000e-9     !Starting time simulation (simulation always starts at 0?)
time_fin = 6.12800e-5     !Total time simulated in [s]  --> 1800 microseconds
t_steps  = 10              !Number of time steps
Src_ON   = 4              !Time at Source is turned ON KIZA ESTO NO SE USE Y BASTA CON SIGNAL Y CAMBIO POR tw
signal   = 1              !Source waveform: 1=step-on; 2=step-off; 3=triangular
initCond = 0              !Initial Condition: 0=None; 1=Double-Line; 2=VMD

# > > > > > > > Name outPut Files
testID   = No_Es_TEM_3D2D_HED_test_01
postpro  = No_Es_TEM_3D2D_HED_test_01
Error    = xxxxxxxx@xxx
Cordina  = xxxxxxxxxxxx
Conecti  = xxxxxxxxxxxx
Profile  = TEM__HED__testFile

# > > > > > > > Source Configuration
#electric Current Vectror
8.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0

#Location
nodalSrc = 1            !Number of nodes will contain the source
13

#Receivers
nodalRec = 1            !Number of nodes as a receiver
500.0 ,0.0, 0.0

