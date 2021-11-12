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
      write(*,"(A29,8X,I6,1X,A11)") ' 3.- Elements:               ', nelem,'   |'
      write(*,"(A29,8X,I6,1X,A11)") ' 4.- Nodal points:           ', nnodes, ' |'
      write(*,"(A29,8X,I6,1X,A11)") ' 5.- DoF per element:        ', ndofn, '  |'
      write(*,"(A29,8X,I6,1X,A11)") ' 6.- Nodes per element:      ', nne, '    |'
      write(*,"(A29,8X,I6,1X,A11)") ' 7.- Total Gauss points:     ', totGp,'   |'
      write(*,"(A29,8X,I6,1X,A11)") ' 8.- Element variabless:     ', nevab,'   |'
      write(*,"(A29,8X,I6,1X,A11)") ' 9.- Total unknowns:         ', ntotv,'   |'
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
    
    subroutine ReadTensors(nr, FileName, difma, conma, reama, force)
      
      character(len=*), parameter    :: fileplace = "~/Dropbox/1.Doctorado/1.Research/Computing/Fortran/ConDifRea/Geo/"
      character (len=*), intent (in) :: FileName
      integer :: status, nr
      double precision, intent(out)  :: difma(3,3,2,2), conma(3,3,2), reama(3,3), force(3) !tensor materials
      
      open (unit = nr, file = fileplace//FileName, status='old', iostat = status)
      
      difma = 0.0
      conma = 0.0
      reama = 0.0
      force = 0.0
      
      !read(nr,1) npoin,nelem,nnode,ngaut,ndofn
      
      if(ndofn.eq.2) then
        read(nr,2) difma(1,1,1,1),difma(1,2,1,1), difma(2,1,1,1),difma(2,2,1,1)
        read(nr,2) difma(1,1,1,2),difma(1,2,1,2), difma(2,1,1,2),difma(2,2,1,2)
        read(nr,2) difma(1,1,2,2),difma(1,2,2,2), difma(2,1,2,2),difma(2,2,2,2)
        
        read(nr,2) conma(1,1,1), conma(1,2,1), conma(2,1,1), conma(2,2,1)
        read(nr,2) conma(1,1,2), conma(1,2,2), conma(2,1,2), conma(2,2,2)
        
        read(nr,2) reama(1,1), reama(1,2), reama(2,1), reama(2,2)
        
        read(nr,3) force(1), force(2)
        
        print*, force 
        
      else if(ndofn.eq.3) then                              
        print*, 'test if'
        read(nr,5) difma(1,1,1,1),difma(1,2,1,1),difma(1,3,1,1)
        read(nr,5) difma(2,1,1,1),difma(2,2,1,1),difma(2,3,1,1)
        read(nr,5) difma(3,1,1,1),difma(3,2,1,1),difma(3,3,1,1)
        
        read(nr,5) difma(1,1,1,2),difma(1,2,1,2),difma(1,3,1,2)
        read(nr,5) difma(2,1,1,2),difma(2,2,1,2),difma(2,3,1,2) 
        read(nr,5) difma(3,1,1,2),difma(3,2,1,2),difma(3,3,1,2)
        
        read(nr,5) difma(1,1,2,2),difma(1,2,2,2),difma(1,3,2,2)
        read(nr,5) difma(2,1,2,2),difma(2,2,2,2),difma(2,3,2,2)
        read(nr,5) difma(3,1,2,2),difma(3,2,2,2),difma(3,3,2,2)
        
        read(nr,5) conma(1,1,1), conma(1,2,1), conma(1,3,1)
        read(nr,5) conma(2,1,1), conma(2,2,1), conma(2,3,1)
        read(nr,5) conma(3,1,1), conma(3,2,1), conma(3,3,1)
        
        read(nr,5) conma(1,1,2), conma(1,2,2), conma(1,3,2)
        read(nr,5) conma(2,1,2), conma(2,2,2), conma(2,3,2)
        read(nr,5) conma(3,1,2), conma(3,2,2), conma(3,3,2)
       
        read(nr,5) reama(1,1), reama(1,2), reama(1,3)
        read(nr,5) reama(2,1), reama(2,2), reama(2,3)
        read(nr,5) reama(3,1), reama(3,2), reama(3,3)
        
        read(nr,3) force(1), force(2), force(3)
        
        print*, force 
        
      end if
      
     ! plate = 0
     ! if(ndofn.eq.-3) then
     !   read(nr,3) young,poiss,thick,force(3)
     !   ndofn=3
     !   plate=1
     ! end if                                                
      
     ! read(nr,4) ksoty,kprec,hnatu,kstab,ktaum,patau,iout
     ! call geodat(coord,ifpre,lnods,posgx,posgy,weigp,unkno)
      
      close (nr)
      
     ! if(plate.eq.1)
     !   call plamat(young,poiss,thick,difma,conma,reama)
     ! endif
      
      !The slash / descriptor begins a new line (record) on output and skips to the next line on input, ignoring any unread information on the current record format(6/) o 6/
      1 format((6/),5(39x,i10,/))
      2 format(39x,2(e15.5),/,39x,2(e15.5))
      3 format(39x,3(e15.5))
      !3 format((12/),39x,3(e15.5))
      4 format(/,2(39x,i10,/),(39x,e15.5,/),2(39x,i10,/),(39x,e15.5,/),(39x,i10,/)/)
      5 format(39x,3(E15.5),/,39x,3(E15.5),/,39x,3(E15.5))
      
      difma(1,1,2,1) = difma(1,1,1,2)
      difma(1,2,2,1) = difma(2,1,1,2)
      difma(2,1,2,1) = difma(1,2,1,2)
      difma(2,2,2,1) = difma(2,2,1,2)
      
      if(ndofn.eq.3) then
        difma(1,3,2,1)=difma(3,1,1,2)
        difma(2,3,2,1)=difma(3,2,1,2)
        difma(3,1,2,1)=difma(1,3,1,2)
        difma(3,2,2,1)=difma(2,3,1,2)
        difma(3,3,2,1)=difma(3,3,1,2)
      end if
      
      
    end subroutine ReadTensors
    
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
    
    function J2D( element_nodes, dN_dxi, dN_deta, Gp)
      implicit none
      
      integer, intent(in)                      :: Gp !esta variable se usara en el lazo principal con el punto de Gauss
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
    
    subroutine DerivativesXY(Gp, Jaco, InvJaco, dN_dxi, dN_deta, Hesxieta, dN_dxy, HesXY)
      
      implicit none
     
      double precision, dimension(DimPr,DimPr),intent(in):: Jaco, InvJaco     
      double precision, dimension(nne,totGp), intent(in) :: dN_dxi, dN_deta
      double precision, dimension(3,nne), intent(in)     :: Hesxieta
      integer, intent(in)                                :: Gp !esta variable se usara en el lazo principal
      double precision, dimension(1,nne) :: Nxi, Neta
      double precision, dimension(2,nne) :: derst
      integer                            :: idime, inode, jdime
      double precision, dimension(2,nne), intent(out)    :: dN_dxy
      double precision, dimension(3,nne), intent(out)    :: HesXY
      
      
      Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)     
      Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)   
      
      
      do idime=1,2
        do inode=1,nne
          dN_dxy = 0.0
          derst(1,inode) = Nxi(1,inode)
          derst(2,inode) = Neta(1,inode)
          do jdime=1,2
            dN_dxy(idime,inode)  = dN_dxy(idime,inode) + Jaco(idime,jdime) * derst(jdime,inode)
          end do
        end do
      end do
      
      !The Hessian matrix
      HesXY = 0.0 
      do inode=1,nnode                                      
        HesXY(1,inode) = InvJaco(1,1)*InvJaco(1,1)*Hesxieta(1,inode)+&
          2.0*InvJaco(1,1)*InvJaco(2,1)*Hesxieta(2,inode)+InvJaco(2,1)*InvJaco(2,1)*Hesxieta(3,inode)
        
        HesXY(2,inode) = InvJaco(1,1)*InvJaco(1,2)*Hesxieta(1,inode) + (InvJaco(1,1)*InvJaco(2,2)+&
          InvJaco(2,1)*InvJaco(1,2)) * Hesxieta(2,inode) + InvJaco(2,1)*InvJaco(2,2)*Hesxieta(3,inode)
        
        HesXY(3,inode) = InvJaco(2,2)*InvJaco(2,2)*Hesxieta(3,inode)+& 
          2.0*InvJaco(2,2)*InvJaco(1,2)*Hesxieta(2,inode)+InvJaco(1,2)*InvJaco(1,2)*Hesxieta(1,inode)
      end do
      
      
    end subroutine DerivativesXY 
    
    
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
    
    function elemSize(Jacobian)
      implicit none
      
      double precision, dimension(DimPr,DimPr), intent(in) :: Jacobian
      double precision :: hx, hy, elemSize
      
      
      hx    = sqrt(Jacobian(1,1)**2+Jacobian(2,1)**2)
      hy    = sqrt(Jacobian(1,2)**2+Jacobian(2,2)**2)
      
      elemSize = hnatu/(min(hx,hy))
     
      return
      
    end function elemsize
    
    
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

    subroutine  GalCon(nne, dvol, basis, dNdxy, ndofn, nevab, amate, rhslo)
      
      implicit none

      double precision, intent(in) :: basis(nnode), derxy(2,nnode)
      double precision, intent(in) :: dvol
      integer, intent(in) :: nne
      integer :: inode, idofn, jevab, jnode, jdofn, i, j
      double precision ::  prod1, prod2, prod3
      double precision, intent(out) :: amate(nevab,nevab), rhslo(nevab)
      ievab=0
      do inode=1,nnode
        do idofn=1,ndofn
          ievab=ievab+1
          jevab=0
          do jnode=1,nnode
            do jdofn=1,ndofn
              jevab=jevab+1
              prod1=0.0
              do i=1,2
                do j=1,2
                  prod1=prod1+ derxy(i,inode) * difma(idofn,jdofn,i,j)* derxy(j,jnode)
                end do
              end do
              prod2=0.0
              do i=1,2
                prod2=prod2+basis(inode)*conma(idofn,jdofn,i)*derxy(i,jnode)
              end do
              prod3=basis(inode)*reama(idofn,jdofn)*basis(jnode)
              amate(ievab,jevab) = amate(ievab,jevab) + (prod1 + prod2 + prod3) * dvol
            end do
          end do
          rhslo(ievab) = rhslo(ievab) + basis(inode) * force(idofn) * dvol
        end do
      end do
      
    end subroutine GalCon 
    
    subroutine AssembleK(K, ke, node_id_map, ndDOF)
      
      implicit none
      real(8), dimension(2*nnodes, 2*nnodes),intent(in out)  :: K !Global Stiffnes matrix debe 
      !                                                         llevar inout por que entra como variable (IN) 
      !                                                         pero en esta funcion se modifica (out)
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
    
    
    subroutine GlobalK( A_K, N, dN_dxi, dN_deta) !Al tener un solo parametro de salida puedo declararla como funcion
      
      implicit none
      
      double precision, dimension(nne,TotGp), intent(in) :: N, dN_dxi, dN_deta
      declarar aqui el Hesxieta
      double precision, dimension(nne)                :: basis
      !double precision, dimension(2*nne, 2*nne)       :: ke
      double precision, dimension(nevab, nevab)       :: Ke
      double precision, dimension(nevab)              :: rhslo
      double precision, dimension(DimPr, dimPr)       :: Jaco, Jinv!, JinvP, JacoP
      double precision                                :: detJ!, detJP
      double precision, dimension(2*DimPr, 2*DimPr)   :: Jb ! aqui tmb es ndofn no DimPr pero 2 para vel y dos para P
      double precision, dimension(2*DimPr, DimPr*nne) :: B  !no es DimPr es ndofn del elemento en cuestion
      double precision, dimension(ndofn,DimPr*DimPr)  :: HJ
      double precision, dimension(ndofn,2*nne)        :: HJB
      double precision, dimension(2*nne,ndofn)        :: HJB_T !Todos estos dos, hablan de los DoF de la velocidad 
      double precision, dimension(2*nne,ndofn)        :: part1 !Todos estos dos, hablan de los DoF de la velocidad
      double precision, dimension(2*nne,2*nne)        :: part2 !Todos estos dos, hablan de los DoF de la velocidad
      double precision, dimension(nevab,nebav)        :: amate
      double precision, dimension(3,nne), intent(in)  :: Hesxieta
      real, dimension(nne,DimPr)   :: element_nodes
      integer, dimension(nne,1)    :: node_id_map
      double precision             :: dvol, hmaxi
      integer                      :: igaus, ielem
      
      
      
      A_K  = 0.0
      cc = reshape([2, 0, 0, 0, 2, 0, 0, 0, 1],[ndofn,ndofn])
      C  = 1.0 * cc
      H  = CompH()
      
      !Setup for K11 block or Kuu
      do ielem = 1, nelem    !lnods loop for K11 block Global K
        !gather
        ke = 0.0
        Jb = 0.0
        amate = 0.0             !amate(nevab*nevab)
        !rhslo = 0.0             !rhslo(nevab)
        call SetElementNodes(ielem, element_nodes, node_id_map)
        !do-loop: compute element (velocity-velocity) stiffness matrix ke
        do igaus = 1, TotGp
          Jaco = J2D(element_nodes, dN_dxi, dN_deta, igaus)
          detJ = m22det(Jaco)
          Jinv = inv2x2(Jaco)
          dvol = detJ *  weigp(igaus,1) 
          call DerivativesXY(igaus, Jaco, Jinv, dN_dxi, dN_deta, Hesxieta, dN_dxy, HesXY)
          hmaxi = elemSize(Jaco) 
          basis = spread(N(:,igaus),dim = 1, ncopies= 1)     
          call galcon(nne, ndofn, nevab, dvol, basis, dN_dxy, Ke, rhslo)
          call taumat(tauma,hmaxi,ndofn) 


          
         part1 = matmul(HJB_T,C)
         part2 = matmul(part1,HJB)
         ke    = ke + part2 * dvolu !
        end do
        
        call AssembleK(A_K, ke, node_id_map, 2) ! assemble global K
        
      end do
     
      
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
