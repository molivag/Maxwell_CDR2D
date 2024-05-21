! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ∂tE + λ∇x∇xE + ∇P + β∇∇⋅E = -∂J/∂t  in  Ω 
! 2d Convection-Diffusion-Reaction simulation             !                ∇⋅E - ɣΔP  = 0       in  Ω 
! Input data file                                         !                       nxE = 0       on  ∂Ω
!                                  MAOG    Bcn, Dic. 2021 ! with:
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !       λ, β, ɣ coefficients.    

$********************************************************************************
PHYSICAL_PROBLEM
$----------------------------------------------------------------------------------
  Maxwell_In_Non_Convex_Domain
 >Cavity_Driven_Flow
  Direct_Current_Electrical_Resistivity_in_2.5-D
  Electric_Field_Excited_By_A_Double_Line_Source
  Horizontal_Electric_Dipole_in_3-D_A_Wrong_capture_of_solution
  Transient_Electromagnetic_in_2.5-D
$----------------------------------------------------------------------------------
END_PHYSICAL_PROBLEM
$*********************************************************************************

# > > > > > > > Model Parameters
ProbType = STAT           !Problem type TIME=transient, other=static
totGp    = 3              !1,4,9 for Q, 1,3,7 for P
exacSol  = 3              !0=None;   1=SinglrSol  ; 2=FullSpace; 3=Algebraic; 4=Double Line
srcRHS   = 0              !0=scalar; 1=SingularSol; 2=Maxwell_Polynom; 3=Lapalace_Polynom
BCsProb  = 2              !1=Ldomain; 2=Maxwell; 3=MaxwellPoly; 4=Lapalace_Poly; 5=Cavity-Driven Flow; 6=Resistivity; 7=Douible-Line
postpro  = 2              !Execution of post-processing routine 1=yes, 2=no 
sigma    = 1.0            !Conductivity of the medium

#***************Geometry
#---------!File .msh that contains the mesh
meshfile = 02gmsh_EM.msh 
view     = xy             !The 2D view x-y (distance) or x-z (depth)
nne      = 3              !Nodes per element Q:4-9; P:3-6
i_exp    = 0              !Exponent of characteristic mesh size 3,4,5 or 6 2^(-i)
hnatu    = 1.0            !Reference element length
refiType = NO             !NONE; PS=Powell-Sabin; CC=Criss-Cross

# > > > > > > > Stabilization
kstab    = 0              !Stabilization: 0(NONE), 1(SUPG), 2(GLS), 3/5(SGS/TG), 4(CG), 6(MVAF)
ktaum    = 1              !Tau matrix: 0, 1, 2 
patau    = 2.0            !Parameter to obtain tau
n_val    = 0.0            !n parameter in exact solution, for simul=1
helem    = 1.0            !Characteristic mesh size (maximum element size among the mesh)
Cu       = 500.0           !Algorithmic constant
ell      = 1000.0          !Constante de longitud  

# > > > > > > > Fourier Transform     
ky_min   = 1.0e-3         !Smallest wave number
ky_max   = 1.0e0          !Grater wave number
tot_ky   = 4              !Total wave numbers
splits   = N              !If the problem is splited then run as kind of parallel
y_iFT    = 0.0            !Position at y-coordinate where the IFT will be computed.
#------------!Filenames for plot the spectrums and store the field after inverse FT
s_spectr = TEM_vs_wavenumber
File_iFT = Real_and_Imaginary_TEM_in_3D 

# > > > > > > > Transient Parameters
theta    = 2              !BDF1=2 ;CN=3; BDF2=4
time_ini = 1.42000e-8     !Starting time simulation (simulation always starts at 0?)
time_fin = 6.12800e-2     !Total time simulated in [s]  --> 1800 microseconds
t_steps  = 300            !Number of time steps
twindow  = 20             !Time window in waveform signal--> ON->small, OFF->large
signal   = 0              !Source waveform: 0=None; 1=step-on; 2=step-off; 3=triangular
initCond = 1              !Initial Condition: 0=None; 1=Double-Line; 2=VMD

# > > > > > > > Name outPut Files
testID   = data_3D2D_TEM_test_05
postpro  = CavityDrivenFlow 
Error    = xxxxxxxx@xxx
Cordina  = xxxxxxxxxxxx
Conecti  = xxxxxxxxxxxx
t_profi  = xx
s_profi  = xx

# > > > > > > > Source Configuration
#electric Current Vectror
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

#Location
nodalSrc = 2            !Number of nodes will contain the source
116
119

#Receivers
nodalRec = 1            !Number of nodes as a receiver
300.0 ,-20.0, 0.0

