module param
  
  implicit none
  
  character(len=14), parameter  :: ElemType = 'Quadrilateral'
  integer, parameter :: DimPr    = 2                        !Dimension del problema 
  integer, parameter :: nelem    = 12                      !Number of lnods
  integer, parameter :: nnodes   = 20                      !Total number of nodal points
  integer, parameter :: nne      = 4                        !Number of coord in the element
  integer, parameter :: ndofn    = 1                        !Degrees of freedom
  integer, parameter :: totGp    = 4                        !1,4,9 for Q, 1,3,4,6 for P 
  integer, parameter :: nevab    = ndofn*nne                !Number of element variables
  integer, parameter :: ntotv    = ndofn*nnodes             !Total number of variables
  integer, parameter :: kstab    = 3                        !Type of stabilization 
  integer, parameter :: ktaum    = 0                        !Type of tau matrix
  real*8,  parameter :: hnatu    = 2.0                      !Reference element length    
  real*8,  parameter :: patau    = 1.0                      !Parameter to obtain tau
  integer, parameter :: maxband  = 6                       !Maximo ancho de banda    
  
  character(len=20), parameter :: File_element  = 'lnods.dat'
  character(len=20), parameter :: File_coord    = 'coord.dat'
  character(len=20), parameter :: File_tensors  = 'tensors.dat'
  !character(len=17), parameter :: File_PostMsh  = 'CDR3d.post.msh'
  !character(len=17), parameter :: File_PostRes  = 'CDR3d.post.res'
  
  double precision, allocatable, dimension(:,:)        :: ngaus, weigp !Verificar si debe ser global
  double precision, dimension(ndofn,ndofn,DimPr,DimPr) :: difma
  double precision, dimension(ndofn,ndofn,DimPr)       :: conma
  double precision, dimension(ndofn,ndofn)             :: reama !Tensor materials
  double precision, dimension(ndofn)                   :: force !Force vector 
  real, dimension(nnodes, DimPr+1)     :: coord
  integer, dimension(nelem, nne + 1)   :: lnods
  integer                              :: nBVs, nBVscol, nband  !Variables usadas en VinculBVs y ApplyBVal si guardan el valor.

end module param
