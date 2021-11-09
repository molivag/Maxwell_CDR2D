module param
  
  implicit none
  
  character(len=14), parameter  :: ElemType = 'Quadrilateral'
  integer, parameter :: DimPr    = 2    !Dimension del problema 
  integer, parameter :: nelem    = 100  !Number of lnods
  integer, parameter :: nnodes   = 121  !Total number of nodal points
  integer, parameter :: nne      = 4    !Number of coord in the element
  integer, parameter :: ndofn    = 3    !Degrees of freedom
  integer, parameter :: totGp    = 4    ! 1,4,9 for Q, 1,3,4,6 for P 


  character(len=20), parameter :: File_element  = 'lnods.dat'
  character(len=20), parameter :: File_coord    = 'coord.dat'
  character(len=29), parameter :: File_PostMsh  = 'CDR3d.post.msh'
  character(len=29), parameter :: File_PostRes  = 'CDR3d.post.res'
  character(len=20), parameter :: File_tensors  = 'tensors.dat'

  double precision, allocatable, dimension(:,:) :: ngaus, weigp !Verificar si debe ser global
  real, dimension(nnodes, DimPr+1)              :: coord
  integer, dimension(nelem, nne + 1)            :: lnods





end module param

