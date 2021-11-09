module library
  use param
  use biunit
  
  
  contains
   
    subroutine GeneralInfo( ) 
      print*, ' '
      print*, '- - - - 2D Cavity Driven Flow Simulation - - - - '
      print*, ' '
      print*,'!==================== GENERAL INFO ===============!'
      write(*,"(A29,8X,A13,3X,A2)") ' 1.- Element type:           ', ElemType,' |'
      write(*,"(A29,8X,I6,1X,A11)") ' 2.- Problem dimension:      ', DimPr, '  |'
      write(*,"(A29,8X,I6,1X,A11)") ' 3.- Total elements:         ', nelem,'   |'
      write(*,"(A29,8X,I6,1X,A11)") ' 4.- Total nodal points:     ', nnodes, ' |'
      write(*,"(A29,8X,I6,1X,A11)") ' 5.- DoF per element:        ', ndofn, '  |'
      write(*,"(A29,8X,I6,1X,A11)") ' 6.- Nodes per element:      ', nne, '    |'
      write(*,"(A29,8X,I6,1X,A11)") ' 7.- Total Gauss points:     ', totGp,'   |'
      write(*,*)' '
      print*,'!============== FILE READING STATUS ============!'
     
    endsubroutine GeneralInfo
    
    subroutine ReadRealFile(UnitNum, FileName, NumRows, NumCols, Real_Array)
      implicit none
      
      integer :: i, j, status, UnitNum, NumRows, NumCols
      character (len=*), intent (in) :: FileName
      character(len=:), allocatable :: fileplace
      real, dimension (1:NumRows, 1:NumCols), intent (out) :: Real_Array
     
      fileplace = "~/Dropbox/1.Doctorado/1.Research/Computing/Fortran/ConDifRea/Geo/"
      
      open (unit = UnitNum, file =fileplace//FileName, status='old', action='read' , iostat = status)
      
      read(UnitNum,*) ((Real_Array(i,j), j=1,NumCols), i=1,NumRows)
      if (status.ne.0) then
        print *, "Status_Real_File ", status
      else
        continue
      end if
      
      close (UnitNum)
      
    end subroutine
    
    subroutine ReadIntegerFile(UnitNum, FileName, NumRows, NumCols, IntegerArray)
     
      integer :: i, j, status
      integer, intent(in)            :: UnitNum, NumRows, NumCols
      character(len=*), parameter    :: fileplace = "~/Dropbox/1.Doctorado/1.Research/Computing/Fortran/ConDifRea/Geo/"
      character (len=*), intent (in) :: FileName
      integer, dimension (1:NumRows, 1:NumCols), intent (out) :: IntegerArray
      
      
      open (unit = UnitNum, file =fileplace//FileName, status='old', action='read' , iostat = status)
      
      read(UnitNum,*) ((IntegerArray(i,j), j=1,NumCols), i=1,NumRows)
      if (status.ne.0) then
        print *, "Status_Int_File  ", status
      else
        continue
      end if
      close (UnitNum)
      
    end subroutine ReadIntegerFile
   ! 
   ! subroutine ReadTensors(UnitNum, FileName, value)
   !   
   !   character(len=*), parameter    :: fileplace = "~/Dropbox/1.Doctorado/1.Research/Computing/Fortran/ConDifRea/Geo/"
   !   character (len=*), intent (in) :: FileName
   !   double precision               :: difma(3,3,2,2), conma(3,3,2), reama(3,3), force(3) !tensor materials
   !   integer :: status, UnitNum, nevab, nr
   !   
   !   
   !   open (unit = nr, file =fileplace//FileName, status='old', action='read' , iostat = status)
   !   
   !   if (status.ne.0) then
   !     print *, "Status_Single_Val", status
   !   else
   !     continue
   !   end if
   !   
   !  
   !   nr = 8
   !   difma = 0.0
   !   conma = 0.0
   !   reama = 0.0
   !   force = 0.0
   !   
   !   !read(nr,1) npoin,nelem,nnode,ngaut,ndofn
   !   
   !   if(ndofn.eq.2) then
   !     read(nr,2) difma(1,1,1,1),difma(1,2,1,1), difma(2,1,1,1),difma(2,2,1,1)
   !     read(nr,2) difma(1,1,1,2),difma(1,2,1,2), difma(2,1,1,2),difma(2,2,1,2)
   !     read(nr,2) difma(1,1,2,2),difma(1,2,2,2), difma(2,1,2,2),difma(2,2,2,2)
   !     
   !     read(nr,2) conma(1,1,1), conma(1,2,1), conma(2,1,1), conma(2,2,1)
   !     read(nr,2) conma(1,1,2), conma(1,2,2), conma(2,1,2), conma(2,2,2)
   !     
   !     read(nr,2) reama(1,1), reama(1,2), reama(2,1), reama(2,2)
   !     
   !     read(nr,3) force(1)      ,force(2)
   !     
   !   else if(ndofn.eq.3) then                              
   !     read(nr,5)
   !     difma(1,1,1,1),difma(1,2,1,1),difma(1,3,1,1),&
   !     difma(2,1,1,1),difma(2,2,1,1),difma(2,3,1,1),&
   !     difma(3,1,1,1),difma(3,2,1,1),difma(3,3,1,1)
   !     
   !     read(nr,5)
   !     difma(1,1,1,2),difma(1,2,1,2),difma(1,3,1,2),&
   !     difma(2,1,1,2),difma(2,2,1,2),difma(2,3,1,2),& 
   !     difma(3,1,1,2),difma(3,2,1,2),difma(3,3,1,2)
   !     
   !     read(nr,5)
   !     difma(1,1,2,2),difma(1,2,2,2),difma(1,3,2,2),&
   !     difma(2,1,2,2),difma(2,2,2,2),difma(2,3,2,2),&
   !     difma(3,1,2,2),difma(3,2,2,2),difma(3,3,2,2)
   !     
   !     read(nr,5)
   !     conma(1,1,1), conma(1,2,1), conma(1,3,1),& 
   !     conma(2,1,1), conma(2,2,1), conma(2,3,1),&
   !     conma(3,1,1), conma(3,2,1), conma(3,3,1)
   !     read(nr,5)
   !     conma(1,1,2)  ,conma(1,2,2)  ,conma(1,3,2),&
   !     conma(2,1,2)  ,conma(2,2,2)  ,conma(2,3,2),&
   !     conma(3,1,2)  ,conma(3,2,2)  ,conma(3,3,2)
   !     read(nr,5)
   !     reama(1,1), reama(1,2), reama(1,3),&
   !     reama(2,1), reama(2,2), reama(2,3),&
   !     reama(3,1), reama(3,2), reama(3,3)
   !     read(nr,3)
   !     force(1), force(2), force(3)
   !     
   !   end if
   !   
   !  ! plate = 0
   !  ! if(ndofn.eq.-3) then
   !  !   read(nr,3) young,poiss,thick,force(3)
   !  !   ndofn=3
   !  !   plate=1
   !  ! end if                                                
   !   
   !  ! read(nr,4) ksoty,kprec,hnatu,kstab,ktaum,patau,iout
   !  ! call geodat(coord,ifpre,lnods,posgx,posgy,weigp,unkno)
   !   
   !   close (UnitNum)
   !   
   !  ! if(plate.eq.1)
   !  !   call plamat(young,poiss,thick,difma,conma,reama)
   !  ! endif
   !   
   !   nevab = ndofn*nne
   !   
   !   !The slash / descriptor begins a new line (record) on output and skips to the next line on input, ignoring any unread information on the current record format(6/) o 6/
   !   1 format(6/),5(39x,i10,/))
   !   2 format(39x,2(e15.5),/,39x,2(e15.5))
   !   3 format(39x,4(e15.5))
   !   4 format(/,2(39x,i10,/),(39x,e15.5,/),2(39x,i10,/),(39x,e15.5,/),(39x,i10,/)/)
   !   5 format(39x,3(e15.5),/,39x,3(e15.5),/,39x,3(e15.5))
   !   
   !   difma(1,1,2,1)=difma(1,1,1,2)
   !   difma(1,2,2,1)=difma(2,1,1,2)
   !   difma(2,1,2,1)=difma(1,2,1,2)
   !   difma(2,2,2,1)=difma(2,2,1,2)
   !   
   !   if(ndofn.eq.3) then
   !     difma(1,3,2,1)=difma(3,1,1,2)
   !     difma(2,3,2,1)=difma(3,2,1,2)
   !     difma(3,1,2,1)=difma(1,3,1,2)
   !     difma(3,2,2,1)=difma(2,3,1,2)
   !     difma(3,3,2,1)=difma(3,3,1,2)
   !   end if
   !   
   !   ntotv=ndofn*npoin
   !   
   ! end subroutine ReadTensors
    
    subroutine ReadMixFile(UnitNum, FileName, NumRows, NumCols, Real_Array)
      implicit none
      
      ! - - - - - - - - - - * * * * * * * * * * - - - - - - - - - -
      ! Rutina que lee un conjunto de datos en el formato indicado
      ! en la etiqueta 22
      !- - - - - - - - - - * * * * * * * * * * - - - - - - - - - -
      
      integer :: i, j, status, UnitNum, NumRows, NumCols
      character(len=*), parameter    :: fileplace = "~/Dropbox/1.Doctorado/1.Research/Computing/Fortran/ConDifRea/Geo/"
      character (len=*), intent (in) :: FileName
      real, dimension (1:NumRows, 1:NumCols), intent (out) :: Real_Array
      
      open (unit = UnitNum, file =fileplace//FileName, status='old', action='read' , iostat = status)
      
      ! read in values
      read(UnitNum,22) ((Real_Array(i,j), j=1,NumCols), i=1,NumRows)
      if (status.ne.0) then
        print *, "Status_Mix_File  ", status
      else
        continue
      end if
      
      22 format(3F13.10)
      close (UnitNum)
      
    end subroutine
    
    subroutine SetElementNodes(elm_num, element_nodes, node_id_map)
      
      implicit none
      
      integer,intent(in)                       :: elm_num ! number of element for each elemental integral in do of K global
      real, dimension(nne,DimPr), intent(out) :: element_nodes
      integer, dimension(nne,1), intent(out)  :: node_id_map
      integer                                  :: i,j, global_node_id
      
      
      element_nodes = 0.0
      node_id_map = 0.0
      
      do i = 1, nne
        global_node_id = lnods(elm_num,i+1)
        do j=1 ,DimPr
          element_nodes(i,j) = coord(global_node_id,j+1)
        end do
        node_id_map(i,1) = global_node_id
      end do
      
    end subroutine SetElementNodes
    
    !subroutine PreassureElemNods(elm_num, pelement_nodes, pnode_id_map)
    !  
    !  implicit none
    !  
    !  integer,intent(in)                       :: elm_num ! number of element for each elemental integral in do of K global
    !  real, dimension(nPne,DimPr), intent(out) :: pelement_nodes
    !  integer, dimension(nPne,1), intent(out)  :: pnode_id_map
    !  integer                                  :: i,j, global_node_id
    !  
    !  
    !  pelement_nodes = 0.0
    !  pnode_id_map = 0.0
    !  
    !  do i = 1, nPne
    !    global_node_id = plnods(elm_num,i+1)
    !    do j=1 ,DimPr
    !      pelement_nodes(i,j) = nodes(global_node_id,j+1)
    !    end do
    !    pnode_id_map(i,1) = global_node_id
    !  end do
    !  
    !end subroutine PreassureElemNods
    
    function J2D( element_nodes, dN_dxi, dN_deta, Gp)
      implicit none
      
      integer, intent(in)                      :: Gp !esta variable se usara en el lazo principal con el numero de punto de gauss para evaluar las integrales elementales
      real, dimension(nne,DimPr), intent(in)  :: element_nodes
      double precision, dimension(nne,totGp), intent(in) :: dN_dxi, dN_deta
      double precision, dimension(DimPr,nne)  :: Basis2D
      double precision, dimension(1,nne)      :: Nxi, Neta
      double precision, dimension(DimPr,DimPr) :: J2D
      
      !con estas instrucciones extraigo la columna de Nx como renglon y lo guardo en Nxi, Gp se
      !ira moviendo conforme la funcion J2D sea llamada en el lazo principal para cada elemento lo mismo para Neta con dN_deta
      Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)
      Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)
      
      !Las siguientes tres lineas realizan de forma implicita el calculo de las derivadas
      !espaciales es decir dN/dx and dN/dy (eq. 5.114 - 5.117). Las derivadas espaciales 
      !no se calcula explicitamente, en su lugar se usa:
        
      !            d/dy = (d/deta)J^-1      (ver eq. 5.77 y 5.109)
      
      Basis2D(1,:) = Nxi(1,:)
      Basis2D(2,:) = Neta(1,:)
      J2D = matmul(Basis2D,element_nodes)
      
      return
    end function J2D

    !function JP2D( element_nodes, dN_dxi, dN_deta, Gp)
    !  implicit none
    !  
    !  integer, intent(in)                      :: Gp !esta variable se usara en el lazo principal con el numero de Gauss points
    !  real, dimension(nPne,DimPr), intent(in)  :: element_nodes
    !  double precision, dimension(nPne,totGp), intent(in) :: dN_dxi, dN_deta
    !  double precision, dimension(DimPr,nPne)  :: Basis2D
    !  double precision, dimension(1,nPne)      :: Nxi, Neta
    !  double precision, dimension(DimPr,DimPr) :: JP2D
    !  
    !  Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)
    !  Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)
    !  
    !  ! Esta forma de 
    !  Basis2D(1,:) = Nxi(1,:)
    !  Basis2D(2,:) = Neta(1,:)
    !  JP2D = matmul(Basis2D,element_nodes) !Aqui se usan directamente las derivadas (eqs 5.114-5.117) de las coordenadas f
    !  
    !  
    !  return
    !end function JP2D
    
    
    function inv2x2(A)
      
      implicit none
      
      double precision, dimension(DimPr, DimPr), intent(in) :: A
      double precision, dimension(DimPr, DimPr)             :: inv2x2
      double precision, dimension(DimPr,DimPr)              :: cofactor
      double precision, parameter :: EPS = 1.0E-10
      double precision            :: det
      
      
      det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      
      if (abs(det) .le. EPS) then
        inv2x2 = 0.0D0
        return
      end if
      
      cofactor(1,1) = +A(2,2)
      cofactor(1,2) = -A(2,1)
      cofactor(2,1) = -A(1,2)
      cofactor(2,2) = +A(1,1)
      
      inv2x2 = transpose(cofactor) / det
      
      return
      
    end function inv2x2
    
    function buildJb(A)
      !Funcion que construye una matriz de 4 x 4 en bloques de 2 para el caso 2D
      
      implicit none
      
      double precision, dimension(DimPr,DimPr), intent (in)   :: A
      double precision, dimension(2*DimPr, 2*DimPr)           :: buildJb
      
      buildJb(1:2,1:2) = A
      buildJb(3:4,3:4) = A
      
    end function
    
    function m22det(A)
      
      implicit none
      double precision :: m22det
      double precision, dimension(2,2), intent(in)  :: A
      
      
      
      m22det =   A(1,1)*A(2,2) - A(1,2)*A(2,1)
      
      return
      
    end function m22det
    
    function CompH()
      implicit None
      
      ! integer :: CompH
      ! integer, dimension(3,4) :: H
      integer, dimension(ndofn ,2*DimPr) :: CompH
      CompH = 0
      CompH(1,1)=1;
      CompH(2,4)=1;
      CompH(3,2)=1;
      CompH(3,3)=1;
      
      ! CompH = H
      ! - - - * * * D U D A * * *
        !no puedo colocar direwctamente el nombre d ela funcion (la funcion misma) como variable global y debo pasarselo a otra variable y esa si ponerla como
        !vbariable global por eso hago el cambio de CompH = H y H esta como variable global. Es Asi?
      ! - - - * * * D U D A * * *
      
      return
      
    end function CompH

    function compBmat(dN_dxi, dN_deta, Gp)

      implicit none

      double precision, dimension(nne,totGp), intent(in) :: dN_dxi, dN_deta
      integer, intent(in) :: Gp

      double precision, dimension(2*DimPr, DimPr*nne) :: compBmat
      double precision, dimension(1, nne)             :: Nxi, Neta
      integer ::  i

      compBmat = 0.0
      Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)
      Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)


      do i=1, nne
        compBmat(1,2*i-1)= Nxi(1,i)
        compBmat(3,2*i)  = Nxi(1,i)
        compBmat(2,2*i-1)= Neta(1,i)
        compBmat(4,2*i)  = Neta(1,i)
      end do

      ! compBmat = B
      return
      ! - - - * * * D U D A * * * - - -
        !En matlab basta con  Nxi(i) aqui quneuq es posible indicar un vector solo con una dimension, no sirve para multiplicarlo.
        !Siempre se debe indicar matriz como un vector fila o vector columna?
      ! - - - * * * D U D A * * * - - -

    end function compBmat

    !function BPmat(dNp_dxi, dNp_deta, Gp)                                               
    !  !Computation of the pressure Strain-Displacement Matrix
    !
    !  implicit none                                                                                      
    !
    !  integer, intent(in) :: Gp                                              
    !  double precision, dimension(nPne,totGp), intent(in) :: dNp_dxi, dNp_deta
    !  double precision, dimension(1*DimPr, 1*nPne) :: Bpmat ! 1 grado de libertad por nodo para los elementos de presion
    !  double precision, dimension(1, nPne)         :: Npxi, Npeta
    !  integer :: i                                                                                                   
    !
    !  Bpmat = 0.0                                                                                                                        
    !  Npxi  = spread(dNp_dxi(:,Gp),dim = 1, ncopies= 1)        
    !  Npeta = spread(dNp_deta(:,Gp),dim = 1, ncopies= 1)     
    !
    !  !Npxi  = Nx(:,Gp)                                                                         
    !  !Npeta = Ny(:,Gp)                                                              
    !
    !  do i=1, nPne                                                                                      
    !    Bpmat(1,i) = Npxi(1,i)                                  
    !    Bpmat(2,i) = Npeta(1,i)                                 
    !  end do                                                    
    !
    !  return                                                    
    !
    !end function BPmat  

    subroutine AssembleK(K, ke, node_id_map, ndDOF)

      implicit none
      real(8), dimension(2*nnodes, 2*nnodes),intent(in out)  :: K !Global Stiffnes matrix debe 
      !                                                                               llevar inout por que entra como variable (IN) 
      !                                                                                pero en esta funcion se modifica (out)
      real(8), dimension(2*nne, 2*nne), intent(in)   :: ke
      integer, dimension(nne,1), intent(in)           :: node_id_map
      integer, intent(in)                              :: ndDOF 
      integer :: i, j, row_node, row, col_node, col !nodal Degrees of Freedom
      
      do i = 1, nne
        row_node = node_id_map(i,1)
        row = ndDOF*row_node - (ndDOF-1)
        
        do j = 1, nne
          col_node = node_id_map(j,1)
          col = ndDOF*col_node - (ndDOF-1)
          K(row:row+ndDOF-1, col:col+ndDOF-1) =  K(row:row+ndDOF-1, col:col+ndDOF-1) + &
          ke((i-1)*ndDOF+1:i*ndDOF,(j-1)*ndDOF+1:j*ndDOF)
        enddo
        
      enddo
      
      return
      
    end subroutine AssembleK
    
    
    !subroutine AssemblyStab(ke, node_id_map, K)
    !  
    !  implicit none
    !  double precision, dimension(,), intent(in out)  :: K !Global Stiffnes matrix debe 
    !  !                                                           llevar inout por que entra como variable (IN) 
    !  !                                                            pero en esta funcion se modifica (out)
    !  double precision, dimension(nPne, nPne), intent(in) :: ke
    !  integer, dimension(nPne,1), intent(in)              :: node_id_map
    !  integer :: i, j, row_node, row, col_node, col, pnode_id !nodal Degrees of Freedom
    !  
    !  !K 
    !  
    !  do i = 1, nPne
    !    row_node = node_id_map(i,1)
    !    pnode_id = pnodes(row_node,2)
    !    row = pnode_id !ndDOF*col_node - (ndDOF-1)
    !    
    !    do j = 1, nPne
    !      col_node = node_id_map(j,1)
    !      pnode_id =  pnodes(col_node,2)
    !      col = pnode_id !ndDOF*col_node - (ndDOF-1)
    !      
    !      K(row,col) =  K(row , col) + ke(i,j)
    !    enddo
    !    
    !  enddo
    !  
    !  
    !  return
    !  
    !end subroutine AssemblyStab
    
    
    subroutine GlobalK( A_K, dN_dxi, dN_deta) !Al tener un solo parametro de salida puedo declararla como funcion
      
      implicit none
      
      double precision, dimension(2*nnodes, 2*nnodes), intent(out) :: A_K  !Global Stiffnes matrix
      double precision, dimension(nne,TotGp), intent(in) :: dN_dxi, dN_deta
      double precision, dimension(2*nne, 2*nne)       :: ke
      double precision, dimension(DimPr, DimPr)       :: Jaco, Jinv!, JinvP, JacoP
      double precision                                :: detJ!, detJP
      double precision, dimension(2*DimPr, 2*DimPr)   :: Jb ! aqui tmb es ndofn no DimPr pero 2 para vel y dos para P
      double precision, dimension(2*DimPr, DimPr*nne) :: B  !no es DimPr es ndofn del elemento en cuestion
      double precision, dimension(ndofn,DimPr*DimPr)  :: HJ
      double precision, dimension(ndofn,2*nne)        :: HJB
      double precision, dimension(2*nne,ndofn)        :: HJB_T !Todos estos dos, hablan de los DoF de la velocidad 
      double precision, dimension(2*nne,ndofn)        :: part1 !Todos estos dos, hablan de los DoF de la velocidad
      double precision, dimension(2*nne,2*nne)        :: part2 !Todos estos dos, hablan de los DoF de la velocidad
      real, dimension(ndofn,2*DimPr)                  :: H
      real, dimension(ndofn,ndofn)                    :: cc, C !Derived from elasticity formulation as Matertial matrix of Hook's Law
      real, dimension(nne,DimPr)   :: element_nodes
      integer, dimension(nne,1)    :: node_id_map
      double precision             :: dvolu  
      integer                      :: gp, e
      
      
      
      A_K  = 0.0
      cc = reshape([2, 0, 0, 0, 2, 0, 0, 0, 1],[ndofn,ndofn])
      C  = 1.0 * cc
      H  = CompH()
      
      !Setup for K11 block or Kuu
      do e = 1, nelem    !lnods loop for K11 block Global K
        ke = 0.0
        Jb = 0.0
        call SetElementNodes(e, element_nodes, node_id_map)
        !do-loop: compute element (velocity-velocity) stiffness matrix ke
        do gp = 1, TotGp
          Jaco = J2D(element_nodes, dN_dxi, dN_deta, gp)
          detJ = m22det(Jaco)
          Jinv = inv2x2(Jaco)
          Jb   = buildJb (Jinv)
          B    = compBmat( dN_dxi, dN_deta, gp)
          HJ   = matmul(H,Jb)
          HJB  = matmul(HJ,B)
         HJB_T = transpose(HJB)
         part1 = matmul(HJB_T,C)
         part2 = matmul(part1,HJB)
         dvolu = detJ *  weigp(gp,1) 
         ke    = ke + part2 * dvolu !
        end do
        
        call AssembleK(A_K, ke, node_id_map, 2) ! assemble global K
        
      end do
      
      !Setup for K12 block or KuP
      !allocate (K12(nnodes*2,n_pnodes),K12_T(n_pnodes,nnodes*2))
      !allocate (K22(,))
      !K12 = 0.0
      !K22 = 0.0
      
      !Tau = (0.99**2 / 4.0 * 1.0)
      !print"(A10,f10.5)",'𝜏 =  ', Tau
      !print*, ' '
      !call ShapeFunctions(gauss_points, nPne, Np, dNp_dxi, dNp_deta)
      !!for-loop: compute K12 block of K
      !do e = 1, nelem
      !  kep = 0.0
      !  Stab = 0.0
      !  call SetElementNodes(e, element_nodes,  node_id_map)
      !  call PreassureElemNods(e, pelement_nodes, pnode_id_map) !-Arreglar esto para q'sea con p en todos argumen
      !  ! for-loop: compute element stiffness matrix kup_e
      !  do gp = 1, TotGp
      !    Jaco  = J2D(element_nodes, dN_dxi, dN_deta, gp)
      !    JacoP = JP2D(pelement_nodes, dNp_dxi, dNp_deta, gp)
      !    detJ  = m22det(Jaco)
      !    detJP = m22det(JacoP)
      !    Jinv  = inv2x2(Jaco)
      !    JinvP = inv2x2(JacoP)
      !    nabP  = Bpmat(dNp_dxi, dNp_deta, gp)  !∇P 
      !    JnabP = matmul(JinvP,nabP) !J^-1 * ∇P 
      !    JP_T  = transpose(JnabP)   !(J^-1 * ∇P)^T
      !    dn    = 0.0
      !    do j = 1, nne
      !      part4(j,:) = [ dN_dxi(j,gp), dN_deta(j,gp) ]  
      !      part5 = reshape([part4(j,:)],[2,1]) 
      !      A =  matmul(Jinv,part5)           
      !      dN(2*j-1:2*j ,1)= A(:,1)         
      !    end do
      !    part6(:,1) = Np(:,gp)
      !    part7 = transpose(part6)
      !    part8 = matmul(dn,part7)
      !    nabTPnabP = matmul(JP_T,JnabP) !∇'δP · ∇P 
      !    kep  = kep + part8 * (detJ*weigp(gp,1)) 
      !    Stab = Stab + nabTPnabP * detJP * weigp(gp,1) ! ∫ (∇'δP : ∇P) dΩ  
      !  end do  
      !    Stab = Stab * Tau 
      !   
      !  ! for-loop: assemble ke into global KP (it mean K12)
      !  do i = 1, nne
      !    row_node = node_id_map(i,1)
      !    row = 2*row_node - 1
      !    do j = 1, nPne
      !      col_node = pnode_id_map(j,1)
      !      pnode_id = pnodes(col_node,2)
      !      col = pnode_id !Aqui puedo sustituir directamente y evito esta misma linea
      !      K12(row:row+1, col) = K12(row:row+1, col) + kep(2*i-1:i*2, j)
      !    end do
      !  end do 
      !  
      !  call AssemblyStab(Stab, pnode_id_map, K22) ! assemble global K
      !  
      !end do
      
      !dimAK = size(A_K,1)
      !symmetric = dimAK - 
      !do i = 1, 2*nnodes
      !  do j = 2*nnodes+1, (2*nnodes)
      !    A_K(i, j) = -K12(i,j-symmetric)
      !  end do
      !end do
      !!========== Lower
      !K12_T = transpose(K12)
      !do i = 2*nnodes+1, (2*nnodes)
      !  do j = 1, 2*nnodes 
      !    A_K(i, j) = -K12_T(i-symmetric,j)
      !  end do
      !end do
      !!========== Filling the stabilization global matrix into A_K ==========
      !do i = 1,  
      !  do j = 1, 
      !    k = 2*nnodes + i
      !    l = 2*nnodes + j
      !    A_K(k, l) = -K22(i,j)
      !  end do
      !end do
      
      !DEALLOCATE(K12)
      !DEALLOCATE(K12_T)
      !DEALLOCATE(Np)
      !DEALLOCATE(dNp_dxi)
      !DEALLOCATE(dNp_deta)
      !DEALLOCATE(K22)
    end subroutine GlobalK
    
    
    
    subroutine SetBounCond( nBVs, nBVscol )
      !========================================================================
      !Esta subroutina revisa todos los nodos de la malla y define el tipo de
      !nodo en la frontera. Abre un archivo en donde comenzara a escribir, 
      ! en la primer columna: el numero de nodo. 
      ! La segunda columna tendra el tipo de nodo
      ! 1 = ux (componente x de la velocidad) 
      ! 2 = uy (componente y de la velocidad) 
      ! 3 = para la presion 
      !La tercera columna asigna el valor correspondiente de la condicion de forntera
      !=========================================================================
      implicit none
                                                     !"/home/maoliva/Codes/ConDifRea_Aca/Geo/"
      character(len=*), parameter :: fileplace ="~/Dropbox/1.Doctorado/1.Research/Computing/Fortran/ConDifRea/Geo/"
      integer, intent(out) :: nBVs, nBVscol
      integer :: ierror, a ,c, i !,b
      real    :: x, y, xmin, xmax, ymin, ymax, xhalf
      
      ! call ReadRealFile(10,"nodes.dat", 341,3, nodes) inicializamos los contadores. Los contadores son para que cada vez
      ! que un if se cumpla, se sume el numero equivalente a los renglones escritos en archivo de texto que se esta creando
      ! y asi se tenga el numero total de nodos en las fronterasi
      
      open(unit=100, file=fileplace//'BVs.dat',Status= 'replace', action= 'write',iostat=ierror)
      
      a = 0
      !b = 0
      c = 0 
      
      xmin = minval(coord(:,2)) !the smallest number in y column
      xmax = maxval(coord(:,2)) !the smallest number in y column
      ymin = minval(coord(:,3)) !the smallest number in y column
      ymax = maxval(coord(:,3)) !the smallest number in y column
      xhalf = xmax/2.0
      
      
      !print*, ' '
      !print*, 'xmin= ', xmin
      !print*, 'xmax= ', xmax
      !print*, 'ymin= ', ymin
      !print*, 'ymax= ', ymax
      !print*, 'xhalf= ', xhalf
      !print*, ' '
      
      
      nBVscol = size(coord,2)     
     
      do i =1, nnodes
        x=coord(i,2)
        y=coord(i,3)
        if(y.eq.ymax) then !top edge: velocity boundary condition
          write(100,50) i, 1, real(1)
          write(100,50) i, 2, real(0)
          a=a+2
          !if(x.eq.xhalf)then !center zero pressure
          !  write(100,50) i,3, real(0)
          !  b=b+1
          !end if
        else if (x.eq.xmin .or. y.eq.ymin .or. x.eq.xmax)then !The other 3 edges
          write(100,50) i, 1, real(0) !x-velocity
          write(100,50) i, 2, real(0) !y-velocity
          c=c+2
        end if
        nBVs = a+c!+b+c
      end do
      
      close(100)
      
      50 format(2I6,f10.3)
      
      
    end subroutine SetBounCond  
    
    
    subroutine ApplyBoundCond( nBVs, BVs, A_K, rhsgl )
      ! - - - - - - - - - - * * * * * * * * * * - - - - - - - 
      ! Set velocity (u) and pressure (p) boundary condition by penalty method
      ! - - - - - - - - - - * * * * * * * * * * - - - - - - - - - -
      implicit none
                          !ndofn
      integer , dimension(nBVs,3), intent(in) :: BVs
      double precision, dimension(2*nnodes, 2*nnodes),intent(in out) :: A_K  !Global Stiffnes matrix
      double precision, dimension(2*nnodes, 1), intent(in out) :: rhsgl
      double precision :: param, coeff
      integer          :: nBVs, i, component, node_id !, pressure_row
      
      !Esencialmente la siguiente instruccion hace: A_K(1*2-1,:) = A_K(1,:) Es decir, obtene el valor maximo de
      !la primera fila de la matriz global K (A_K). No le veo el caso pero lo dejamos asi.
      param = maxval(A_K(int(BVs(1,1))*2-1,:))
      coeff = abs(param) * 1.0E7
      
      print*, 'param', param
      print*, 'coeff', coeff
      
      !pressure_row = 2*nnodes
      
      do i =1, nBVs
        node_id   = BVs(i,1) !se pone este int() pq la 1a y 2a col de BVs esta leida como integer pero 
        component = BVs(i,2)!la matriz completa esta declarada como real en esta subroutina y en el main.
        if ( component .le. 2 ) then
          A_K(2*node_id-2+component, 2*node_id-2 +component) = coeff
          rhsgl( 2*node_id-2+component, 1) = BVs(i,3)*coeff 
        !else                                                     
          !pnode_id = pnodes(node_id,2)
          !A_K(pressure_row+pnode_id, pressure_row + pnode_id) = coeff
          !rhsgl(pressure_row+pnode_id,1) = BVs(i,3)*coeff !3 por que la columna 3 esta el valor de la condicon de forntera
        end if
      end do
      
    end subroutine ApplyBoundCond
    
    subroutine MKLfactoResult( value )
      implicit none
      
      integer :: value, val
      character(len=34) :: text
      character(len=48) :: text2
      
      text  = '   *FACTORIZATION DONE WITH STATUS'
      text2 = '   *THE FACTORIZATION HAS BEEN COMPLETED, BUT U('
      if ( value .eq. 0 ) then
        write(*, 101) text, value, ', THE EXECUTION IS SUCCESSFUL.'
      elseif(value .lt. 0 )then
        val = abs(value)
        write(*, 102) '    THE',val,'-TH PARAMETER HAD AN ILLEGAL VALUE.'
      elseif(value .gt. 0 )then
        write(*, 103) text2, value,',',value,') IS EXACTLY SINGULAR.'
        print*,'   DIVISION BY 0 WILL OCCUR IF YOU USE THE FACTOR U FOR SOLVING A SYSTEM'
        print*,'   OF LINEAR EQUATIONS.'
      endif
      print*, ' '
      
      101 format (A, 1x, I1, A)
      102 format (A, I4, A)
      103 format (A, I3, A, I3, A)
      
    end subroutine MKLfactoResult
    
    subroutine MKLsolverResult( value )
      implicit none
      
      integer :: value, val
      character(len=30) :: text
      character(len=35) :: text2
      text =  '   *SYSTEM SOLVED WITH STATUS'
      text2 = '-TH PARAMETER HAD AN ILLEGAL VALUE.'
      
      if ( value .eq. 0 ) then
        write(*,101) text, value, ', THE EXECUTION IS SUCCESSFUL.'
      elseif(value .lt. 0 )then
        val = abs(value)
        write(*,102) '    THE',val, text2
      endif
      print*,' '
      
      101 format (A, 1x, I1, A)
      102 format (A, I3, A)
      
    end subroutine MKLsolverResult
    
    subroutine writeMatrix(Matrix, unit1, name1, Vector, unit2, name2)
      implicit none
      
      character(len=*), parameter    :: fileplace = "~/Dropbox/1.Doctorado/1.Research/Computing/Fortran/ConDifRea/Res/"
      character(*) :: name1, name2
      integer :: i, j, mrow, ncol, unit1, unit2
      double precision, dimension(2*nnodes ,2*nnodes ), intent(in) :: Matrix
      double precision, dimension(2*nnodes ,1), intent(in) :: Vector
      
      100 format (900E20.12)
      
      mrow = 2*nnodes 
      ncol = 2*nnodes
      open(unit=unit1, file= fileplace//name1, ACTION="write", STATUS="replace")
      
      do i=1,2*nnodes 
        write(unit1, 100)( Matrix(i,j) ,j=1,2*nnodes)
      end do
      close(unit1)
      
      open(unit=unit2, file= fileplace//name2, ACTION="write", STATUS="replace")
      do i=1,2*nnodes 
        write(unit2, 100) Vector(i,1)
      end do
      close(unit2)
      
    end subroutine writeMatrix
    
    subroutine PosProcess(solution, nameFile1, activity)
      
      implicit none
      
      character(len=*), parameter    :: fileplace = "~/Dropbox/1.Doctorado/1.Research/Computing/Fortran/ConDifRea/Pos/"
      real*8, dimension(2*nnodes, 1), intent(in) :: solution
      character(*), intent(in)                             :: nameFile1, activity
      double precision, dimension(1, 2*nnodes)   :: solution_T
      double precision, dimension(1,nnodes)               :: xcor, ycor
      integer      :: ipoin   
      
      solution_T = transpose(solution)
      xcor  = spread(coord(:,2),dim = 1, ncopies= 1)
      ycor  = spread(coord(:,3),dim = 1, ncopies= 1)
      
      
      open(unit=555, file= fileplace//nameFile1, ACTION="write", STATUS="replace")
      
      if(activity == "msh")then !quitar este if y acomodar el numero de unidad
        
        write(555,902) 'MESH', '"Cavity"', 'dimension', DimPr, 'ElemType', ElemType, 'Nnode', nne
        write(555,"(A)") '#2D Cavity Driven Flow Results' 
        write(555,900) '#Element tipe: ', ElemType,'/',ElemType 
        write(555,"(A)")'Coordinates'
        write(555,"(A)") '#   No        X           Y'
        do ipoin = 1, nnodes
          write(555,906) ipoin, xcor(1,ipoin), ycor(1,ipoin)
        end do
        write(555,"(A)") 'End Coordinates'
        write(555,"(A)") 'Elements'
        do ipoin = 1, nelem
          write(555,908) lnods(ipoin,:) 
        end do
        write(555,"(A)") 'End Elements'
        close(555)
        
      elseif(activity == "res")then
        write(555,"(A)") 'GiD Post Results File 1.0'
        write(555,"(A)") '#2D Cavity Driven Flow Results' 
        write(555,900) '#Element tipe: ', ElemType,'/',ElemType 
        write(555,"(A)") 'Result "Velocity Components" "Velocity" 0 Vector OnNodes'
        write(555,"(A)") 'ComponentNames "Ux" "Uy" "Uz" "" '
        write(555,"(A)") 'Values'
        ! se escribe el res de las componentes de la velocidad
        write(555,910) 
        do ipoin = 1, nnodes
          write(555,912) ipoin, solution_T(1, 2*ipoin-1), solution_T(1,2*ipoin)
        end do
        write(555,"(A)") 'End Values'
       ! write(555,"(A)") 'Result "Pressure" "Pressure" 0 Scalar OnNodes'
       ! write(555,"(A)") 'ComponentNames "" '
       ! write(555,"(A)") 'Values'
       ! ! se escribe el res de la presion 
       ! write(555,914)
       ! do ipoin = 1, nnodes
       !   pnode_id = pnodes(ipoin,2)
       !   write(555,916) ipoin, solution(prow+pnode_id, 1)  
       ! end do
       ! write(555,"(A)") 'End Values'
        close(555)
      else
        write(*,"(A)") ' "Activity" must be "msh" or "res" '
        close(555)
        stop
      end if
      
      
      900 format(A15, A13, A1, A13)
      902 format(A4,1x,A8,1X,A9,1X,I1,1X,A8,1X,A13,A6,1X,I1)
      906 format(I7,2(3x,f9.4)) !format for msh           
      908 format(9(2x,I7) )
      910 format('#',3x,'No    ' 3x, ' Ux ', 8x, ' Uy')
      912 format(I7,2x,2f12.5) !format for res velocity
      914 format('#',3x,'No'     9x, 'P')
      916 format(I7,2x,f12.5)  !format for res pressure
      
    end subroutine PosProcess
    
   
    
  !Fin de contains


end module library
