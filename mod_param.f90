module param
  
  implicit none
  
  character(len=14), parameter  :: ElemType = 'Quadrilateral'
  integer, parameter :: DimPr     = 2     !Dimension del problema 
  integer, parameter :: Nelem     = 100  !Number of elements
  integer, parameter :: n_nodes   = 121 !Total number of velocity nodes
  integer, parameter :: n_pnodes  = 121  !Total number of preasure nodes MAXVAL(pnodes,2)
  integer, parameter :: nUne      = 4     !Number of velocity nodes in the element
  !integer, parameter :: nPne      = 4     !Number of preasure nodes in the element
  integer, parameter :: Dof       = 3     !Degrees of fredoom: 2 of velocity + 1 of preasure
  integer, parameter :: totGp     = 4     ! 1,4,9 for Q, 1,3,4,6 for P 
  real(4), parameter :: materials = 1.0 


  character(len=20), parameter :: File_element  = 'elements.dat'
  character(len=20), parameter :: File_nodes    = 'nodes.dat'
  character(len=29), parameter :: File_PostMsh  = 'CDR3d.post.msh'
  character(len=29), parameter :: File_PostRes  = 'CDR3d.post.res'

  double precision, allocatable, dimension(:,:) :: gauss_points, gauss_weights !Verificar si debe ser global
  integer, dimension(Nelem, nUne + 1)    :: elements
  real,    dimension(n_nodes, DimPr + 1) :: nodes

end module param

