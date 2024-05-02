module geometry
use param
  
  implicit none

  integer :: mxnod, mxnow, mxelw, mxpow
  
  contains
    
    subroutine readMesh(file_mesh)
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      !                                                                                   !
      ! subrutina que lee todos los parametros de entrada para la simulacion,             !
      ! la geometria, lista de nodos, coordenadas y parametros de estabilizacion          !
      !                                                                                   !
      ! The MSH file format version 1 is Gmsh’s original native mesh file format,         !
      ! now superseded by the format described in MSH file format.                        !
      ! It is defined as follows:                                                         !
      !                                                                                   !
      !                                                                                   !
      ! $ENDNOD                                                                           !
      !    $ELM                                                                           !
      !    number-of-elements                                                             !
      !    elm-number elm-type reg-phys reg-elem number-of-nodes node-number-list         !
      !    …                                                                              !
      ! $ENDELM                                                                           !
      !                                                                                   !
      ! -- reg-phys                                                                       !
      !         is the tag of the physical entity to which the element belongs;           !
      !         reg-phys must be a positive integer, or zero. With this will be           !
      !         easily identified what element belong to what latyer or surface           !
      !                                                                                   !
      !   --reg-elem                                                                      !
      !         is the tag of the elementary entity to which the element belongs;         !
      !         reg-elem must be a positive (non-zero) integer.                           !
      !         2 mean: element belong to surface (surface identify as entity 2)          !
      !                                                                                   !
      !                                                                                   !
      !                                                                                   !
      !                                                                                   !
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      
      implicit none
      
      character(len=*), parameter                   :: fileplace = "Msh/"
      character(len=:), allocatable, intent(in)     :: file_mesh
      character(len=180)                            :: msg
      double precision, allocatable, dimension(:,:) :: cord3D,coorw
      integer,          allocatable, dimension(:,:) :: lnodw, list_new_nodes
      integer,          allocatable, dimension(:)   :: mesh_conduc, chk_physical_region
      integer                                       :: npoiw,nelew,nnodw, npoif,ielem, jpoin, idime
      integer                                       :: initOrderElem, reg_phys, i,j, stat, dmy, gmsh_nne
      
      
      open(5, file=fileplace//file_mesh, status='old', action='read',IOSTAT=stat, IOMSG=msg)
      
      !do i=1,skipline
      !read(5,*) !se salta todas las lineas del input file hasta donde comienza la malla
      !end do
      read(5,*,IOSTAT=stat, IOMSG=msg) !se salta la lineas del archivo .msh
      call checkStatus(4,stat,msg)
      
      read(5,*,IOSTAT=stat, IOMSG=msg) initNodes
      call checkStatus(5,stat,msg)
      
      nnodes = initNodes
      allocate(coord(Dimpr,initNodes))
      coord = 0.0 
      if(view == 'xz')then
        print*, 'plano xz'
        allocate(cord3D(3,initNodes))
        coord = 0.0
        do i=1,nnodes !number of total nodes
          read(5,*,iostat=stat,iomsg=msg) jpoin,(cord3D(idime,jpoin), idime =1,3 )
        end do
        do j =1 , initNodes
          coord(1,j) = cord3D(1,j) !coordenada x
          coord(2,j) = cord3D(3,j) !coordenada z
        enddo
        deallocate(cord3D)
      else
        print*, 'plano xy'
        read_coordinates:do i=1,nnodes !number of total nodes
          read(5,*,iostat=stat,iomsg=msg) jpoin,(coord(idime,jpoin), idime =1,DimPr )
        end do read_coordinates
      endif
      do i=1,2
        read(5,*) !se salta todas las lineas entre nodes y elements y comienza a leer los elementos 
      end do
      
      !do jpoin = 1, nnodes
      !  print'(i5, 2x,2(F16.9))', jpoin, (coord(idime,jpoin),idime=1,DimPr)
      !end do
      
      read(5,*) initElem 
      call checkStatus(8,stat,msg)
      nelem = initElem
      print*, 'Se asigna mesh_conductivity'
      allocate( lnods(initElem,initnne), mesh_conductivity(nelem),  chk_physical_region(nelem) )
      lnods = 0.0
      read_conectivity:do i=1,nelem
        read(5,*,iostat=stat,iomsg=msg) ielem, initOrderElem, reg_phys, dmy, gmsh_nne, (lnods(ielem,j), j =1,initnne)
        if(gmsh_nne.ne.initnne)then
          print*," error in Number of Nodes in the element, while reading lnodes in geometry module "
          stop
        endif
        call checkStatus(9,stat,msg)
        chk_physical_region(i) = reg_phys
      end do read_conectivity
      read(5,*) !se salta todas las lineas entre nodes y elements y comienza a leer los elementos 
      !Aqui fallo al leer triangulos cuando debio decirme que los nne no coincidian o que los GP no coincidian
      !Por lo tanto, agregar un algo que en lugar de que falle aqui, verifique que son triangulos o cuadrados
      !lo que se esta leyendo
      
      !do i = 1, nelem
      !  print'(i5, 2x,5(I6))', i, (lnods(i,j),j=1,initnne)
      !end do
      
      close(5)
      
      select case(initOrderElem)
      case(2)
        initElemType = 'TRIA'
        gmsh_nne = 3
        if(initnne.eq.gmsh_nne)then
          continue
        else
          write(*,'(A)') 'Error in nne, program stop inside Geometry module'
          write(*,'(A40, I6)') 'Number of nodes in the element must be: ', gmsh_nne
          stop
        endif
      case(9)
        initElemType = 'TRIA' !2nd Order (6node) Triangle 
        gmsh_nne = 6
        if(initnne.eq.gmsh_nne)then
          continue
        else
          write(*,'(A)') 'Error in nne, program stop inside Geometry module'
          write(*,'(A40, I6)') 'Number of nodes in the element must be: ', gmsh_nne
          stop
        endif
      case(3)
        initElemType = 'QUAD'
        gmsh_nne = 4
        if(initnne.eq.gmsh_nne)then
          continue
        else
          write(*,'(A)') 'Error in nne, program stop inside Geometry module'
          write(*,'(A40, I6)') 'Number of nodes in the element must be: ', gmsh_nne
          stop
        endif
      case(10)
        initElemType = 'QUAD' !2nd Order (9node) Quadrilateral 
        gmsh_nne = 9
        if(initnne.eq.gmsh_nne)then
          continue
        else
          write(*,'(A)') 'Error in nne, program stop inside Geometry module'
          write(*,'(A40, I6)') 'Number of nodes in the element must be: ', gmsh_nne
          stop
        endif
      case default
        write(*,'(A)') 'error >>>> The initElemType must be 2, 9, 3 or 10'
        stop
      end select
      
      initnevab = ndofn*initnne
      initntotv = ndofn*initNodes
      
      if(refiType.eq.'NO')then
        ElemType = initElemType
        nevab    = initnevab
        ntotv    = initntotv
        nne      = initnne
        
        !Sino hay refinado entonces el reg_phys es el proveniente de gmsh directamente
        mesh_conductivity = chk_physical_region
        goto 101
        
      elseif(refiType.eq.'PS'.or.refitype.eq.'CC')then
        if(refiType.eq.'PS')then
          mxnow = 5 
          mxelw = nelem*6
          mxpow = initNodes+nelem*4
        elseif(refiType.eq.'CC')then
          mxnow = 5 
          mxelw = nelem*4
          mxpow = initNodes+nelem+1 ! el +1 es si acaso
        else
          write(*,*) 'No Refinement type defined'
          write(*,*) '>>>>>STOP on L191 GEOMETRY'
          stop
        endif
        deallocate(mesh_conductivity)
        
        !***  Undertakes the mesh change
        call AddNodes(refiType,npoiw,nnodw,coorw,lnodw)
        
        !***  Checks if there are repeated nodes and output of results
        call SplitElem(refiType,chk_physical_region, nelew,lnodw, mesh_conduc, list_new_nodes) 
        
        !***  Checks if there are repeated nodes and reallocate coord and lnods
        call checkMesh(coorw,list_new_nodes,nnodw,nelew,npoiw,mesh_conduc, npoif)
        deallocate( chk_physical_region, mesh_conduc, coorw,list_new_nodes)
        
        !* Recounting of nodes and elements after the refination
        nnodes = npoif
        nelem  = nelew 
        nne    = nnodw
        nevab  = ndofn*nne
        ntotv  = ndofn*nnodes
        
      end if
      
      101 continue
      
      ! open(135, file='./test_asign_mediumPS.dat', status='old', action='write',IOSTAT=stat, IOMSG=msg)
      ! do i=1,nelem
      !   if(i.lt.10)then
      !     write(135,"(I2,2x,4(I5,1x))") i, (lnods(i,j), j=1,nne)
      !   else
      !     write(135,"(I0,2x,4(I5,1x))") i, (lnods(i,j), j=1,nne)
      !   endif
      ! end do
      ! close(135)
      
      ! open(115, file='./media_meshPS.dat', status='old', action='write',IOSTAT=stat, IOMSG=msg)
      ! do i=1,nelem
      !   if(i.lt.10)then
      !     write(115,"(I2,2x,I5)") i, mesh_conductivity(i)
      !   else
      !     write(115,"(I0,2x,I5)") i, mesh_conductivity(i)
      !   endif
      ! end do
      ! close(135)
      
      
    end subroutine readMesh
    !
    !*************************************************************************     
    !
    subroutine AddNodes(refiType,npoiw,nnodw,coorw,lnodw)
      !******************************************************************************
      !
      !**** This program changes the element type of a finite element mesh.
      !**** The possibilities are the following:      
      !
      !**** For NDIME = 2
      !      
      !****     NNODE = 3     ---->    NNODE = 7    (more nodes are created)
      !****     NNODE = 9     ---->    NNODE = 3    (same number of nodes, more elements)
      !
      !****     NNODE = 4     ---->    NNODE = 5    (more nodes are created)
      !****     NNODE = 5     ---->    NNODE = 3    (triangular elements of 3nodes)
      !
      !******************************************************************************
      !       
      character(len=2),intent(in)                                :: refitype
      integer                                                    :: idime, inode, ielem, ipoin 
      ! double precision, intent(out), dimension(DimPr,mxpow)    :: coorw
      ! integer         , intent(out), dimension(mxelw,mxnow)    :: lnodw
      double precision, allocatable, dimension(:,:), intent(out) :: coorw
      integer         , allocatable, dimension(:,:), intent(out) :: lnodw
      integer         , intent(out)                              :: npoiw, nnodw 
      
      allocate( coorw(DimPr,mxpow),lnodw(mxelw,mxpow) )
      
      !***  initializations
      do idime=1,DimPr
        do ipoin=1,nnodes
          coorw(idime,ipoin)=coord(idime,ipoin)
        end do
      end do
      npoiw=nnodes
      nnodw=3 !is a fix value because both PS & CB, has as a final element a 3-node ∆ element
      
      do ielem=1,nelem
        !*** 2D: NNODE = 3 --> NNODE = 4 6 & 7
        if(refiType.eq.'PS')then
        !Si quiero que el refinado PS parta tambien de cuadrados entonces primero debo poner las lineas
        !de dividir un cuadrado en 4 triangulos (es decir abajo la que crea el criss-cross pero en lugar
        !de dividir en 4, divido en 2 y con ello ya tendre triangulos que luego seguira en el ciclo !here
        !que agregara 3 nodos mas y luego uno al centro es decir todo queda igual solo agregar aqui previ
        !al primer ciclo do, el ciclo de dividir un cuadrado en dos elementos y ya
          if(ielem.eq.1)write(*,'(a)') 'Powell-Sabin refinement type selectd use triangle as a initial element'
          if(initnne.ne.3)then 
            write(*,'(a)') 'PS refinment not compatible with quadrilateral element'
            write(*,'(a)') '>>> Verify PS must be nne=3'
            stop
          else
            continue
            !if (nnodw.ge.6) then 
              !Para construir el elemento de 7 nodos. Primero hay que construir el de 6
              !por eso la codicion dice mayor o igual, si pongo 7 entonces 
              !entrara al de 6 y luego al de 7
              do idime=1,2 !here
                coorw(idime,npoiw+1)= 0.5*(coord(idime,lnods(ielem,1))+coord(idime,lnods(ielem,2)))
                coorw(idime,npoiw+2)= 0.5*(coord(idime,lnods(ielem,2))+coord(idime,lnods(ielem,3)))
                coorw(idime,npoiw+3)= 0.5*(coord(idime,lnods(ielem,3))+coord(idime,lnods(ielem,1)))
              end do
              do inode=1,3
                lnodw(ielem,  inode)=lnods(ielem,inode)
                lnodw(ielem,3+inode)=npoiw+inode
              end do
              npoiw=npoiw+3
            !if(nnodw.eq.7) then
              do idime=1,2
                coorw(idime,npoiw+1)= (1.0/3.0)*(coord(idime,lnods(ielem,1)) + &
                &                     coord(idime,lnods(ielem,2))+coord(idime,lnods(ielem,3)))
              end do
              lnodw(ielem,7)=npoiw+1
              npoiw=npoiw+1
            !endif
          endif
          
        elseif(refiType.eq.'CC')then
          if(initnne.ne.4)then
            
            write(*,'(a)') 'CC refinment not compatible with triangular element'
            write(*,'(a)') '>>> Verify CC must be nne=4'
            stop
          else
            continue
            !  !Creation of one more node at the center of 4 node quadrilateral
            do idime=1,2
              coorw(idime,npoiw+1)= 0.25*(coord(idime,lnods(ielem,1))+coord(idime,lnods(ielem,2))&
              &                         +coord(idime,lnods(ielem,3))+coord(idime,lnods(ielem,4)))
            end do
            
            lnodw(ielem  ,1)=lnods(ielem,1)
            lnodw(ielem  ,2)=lnods(ielem,2)
            lnodw(ielem  ,3)=lnods(ielem,3)
            lnodw(ielem  ,4)=lnods(ielem,4)
            lnodw(ielem  ,5)=npoiw+1
            npoiw=npoiw+1
            
          endif
        else
          write(*,'(A)')' --No refinment selected-- '
          goto 101
        end if
        
      end do
      
      if(npoiw.gt.mxpow) then
        write(*,'(a,i5)') 'Increase MXPOW to ' , npoiw
        stop
      end if
      
      101 continue
      
    end subroutine AddNodes 
    !
    !*************************************************************************     
    !
    subroutine SplitElem(refiType,chk_physical_region, nelew,lnodw, mesh_conduc, list_new_nodes)

      implicit none 
      
      character(len=2), intent(in)                          :: refitype
      integer, allocatable, dimension(:)  , intent(in)      :: chk_physical_region
      integer, allocatable, dimension(:,:), intent(in)      :: lnodw! -------> ya esta allocate en AddNodes
      integer                                               :: ielem, reg_phys
      integer, allocatable, dimension(:,:)                  :: lnod_add! --------> Se asigna aqui 
      ! integer   , dimension(mxelw,mxnow) , intent(out)    :: lnodw
      ! integer   , dimension(mxelw)       , intent(out)    :: mesh_conduc
      integer, allocatable, dimension(:,:) , intent(out)    :: list_new_nodes! -------> se asigna aqui
      integer, allocatable, dimension(:)   , intent(out)    :: mesh_conduc
      integer                              , intent(out)    :: nelew 
      
      
      allocate( mesh_conduc(mxelw) )
      allocate( lnod_add(mxelw,mxpow), list_new_nodes(mxelw,mxpow) )
     
      lnod_add = lnodw 
      list_new_nodes = 0
      
      nelew=nelem
      do ielem=1,nelem
        !*** 2D: NNODE = 3 --> NNODE = 4 6 & 7
        if(refiType.eq.'PS')then
          if(initnne.ne.3)then 
            write(*,'(a)') 'PS refinment not compatible with quadrilateral element'
            write(*,'(a)') '>>> Verify PS must be nne=3'
            stop
          else
            !identifico que reg_phys es respecto del elemento en cuestion
            reg_phys = chk_physical_region(ielem) 
            
            continue
            !Split of a 7-node ∆ element into 6 ∆ elements of 3-nodes. 
            !  ^
            !  |        3
            !  |        o
            !  |       / \
            !  |      /   \
            !  Y    6o  7  o5
            !  |    /   o   \
            !  |   /         \
            !  |  o-----o-----o
            !  |  1     4     2
            !  |
            !  +--------X-------->
            list_new_nodes(ielem  ,1)=lnod_add(ielem,1)
            list_new_nodes(ielem  ,2)=lnod_add(ielem,4)
            list_new_nodes(ielem  ,3)=lnod_add(ielem,7)!1er elemento generado 
            list_new_nodes(nelew+1,1)=lnod_add(ielem,1)
            list_new_nodes(nelew+1,2)=lnod_add(ielem,7)
            list_new_nodes(nelew+1,3)=lnod_add(ielem,6)!2º elemento generado
            list_new_nodes(nelew+2,1)=lnod_add(ielem,7)
            list_new_nodes(nelew+2,2)=lnod_add(ielem,3)
            list_new_nodes(nelew+2,3)=lnod_add(ielem,6)!3er elemento generado 
            list_new_nodes(nelew+3,1)=lnod_add(ielem,4)
            list_new_nodes(nelew+3,2)=lnod_add(ielem,2)
            list_new_nodes(nelew+3,3)=lnod_add(ielem,7)!4º elemento generado 
            list_new_nodes(nelew+4,1)=lnod_add(ielem,2)
            list_new_nodes(nelew+4,2)=lnod_add(ielem,5)
            list_new_nodes(nelew+4,3)=lnod_add(ielem,7)!5º elemento generado 
            list_new_nodes(nelew+5,1)=lnod_add(ielem,5)
            list_new_nodes(nelew+5,2)=lnod_add(ielem,3)
            list_new_nodes(nelew+5,3)=lnod_add(ielem,7)!6º elemento generado
            
            mesh_conduc(ielem  ) = reg_phys
            mesh_conduc(ielem  ) = reg_phys
            mesh_conduc(ielem  ) = reg_phys
            mesh_conduc(nelew+1) = reg_phys
            mesh_conduc(nelew+1) = reg_phys
            mesh_conduc(nelew+1) = reg_phys
            mesh_conduc(nelew+2) = reg_phys
            mesh_conduc(nelew+2) = reg_phys
            mesh_conduc(nelew+2) = reg_phys
            mesh_conduc(nelew+3) = reg_phys
            mesh_conduc(nelew+3) = reg_phys
            mesh_conduc(nelew+3) = reg_phys
            mesh_conduc(nelew+4) = reg_phys
            mesh_conduc(nelew+4) = reg_phys
            mesh_conduc(nelew+4) = reg_phys
            mesh_conduc(nelew+5) = reg_phys
            mesh_conduc(nelew+5) = reg_phys
            mesh_conduc(nelew+5) = reg_phys
            nelew=nelew+5
            
            ElemType = initElemType
           
          endif
          
        elseif(refiType.eq.'CC')then
          if(initnne.ne.4)then
            
            write(*,'(a)') 'CC refinment not compatible with triangular element'
            write(*,'(a)') '>>> Verify CC must be nne=4'
            stop
          else
            continue
            !  !Creation of one more node at the center of 4 node quadrilateral
            !  
            !  !Split of a 5-node element into four 3 nodes ∆ elements. 
            !  !  ^
            !  !  | 4       3
            !  !  | o-------o
            !  !  | | \   / |
            !  !  | |  \5/  |
            !  !  Y |   o   |        Is a square compund by 4 three-nodes triangle      
            !  !  | |  / \  |
            !  !  | | /   \ |
            !  !  | o-------o
            !  !  | 1       2
            !  !  +------X------->
            !  
            list_new_nodes(ielem  ,1)=lnod_add(ielem,1)
            list_new_nodes(ielem  ,2)=lnod_add(ielem,2)
            list_new_nodes(ielem  ,3)=lnod_add(ielem,5)
            list_new_nodes(nelew+1,1)=lnod_add(ielem,1)
            list_new_nodes(nelew+1,2)=lnod_add(ielem,5)
            list_new_nodes(nelew+1,3)=lnod_add(ielem,4)
            list_new_nodes(nelew+2,1)=lnod_add(ielem,5)
            list_new_nodes(nelew+2,2)=lnod_add(ielem,3)
            list_new_nodes(nelew+2,3)=lnod_add(ielem,4)
            list_new_nodes(nelew+3,1)=lnod_add(ielem,2)
            list_new_nodes(nelew+3,2)=lnod_add(ielem,3)
            list_new_nodes(nelew+3,3)=lnod_add(ielem,5)
            nelew=nelew+3
            
            ElemType = 'TRIA'
            
          endif
        else
          write(*,'(A)')' --No refinment selected-- '
          goto 101
        end if                                            ! nne
        
        
      end do                                              ! ielem=1,nelem
      deallocate(lnod_add)
      
      if(nelew.gt.mxelw) then
        write(*,'(a,i5)') 'Increase MXELW to ' , nelew
        stop
      end if
      
      101 continue
      
    end subroutine SplitElem
    !
    !*************************************************************************     
    !
    subroutine checkMesh(coorw,list_new_nodes,nnodw,nelew,npoiw, mesh_conduc, npoif)
      
      !implicit     real*8 (a-h,o-z)
      
      implicit none !integer :: (i,j)
      
      ! integer         , intent(in) :: list_new_nodes(mxelw,mxnow), mesh_conduc(mxelw)
      ! double precision, intent(in)  :: coorw(DimPr,mxpow)
      
      integer         , allocatable, dimension(:,:),intent(in) :: list_new_nodes
      integer         , allocatable, dimension(:)  ,intent(in) :: mesh_conduc
      double precision, allocatable, dimension(:,:),intent(in) :: coorw
      integer                                      ,intent(in) :: nnodw, nelew, npoiw
      integer                                                  :: ipoiw, ielew, inodw, idime, kpoiw, jpoiw 
      integer         , allocatable, dimension(:)              :: lpoiw
      double precision                                         :: x0, x1, y0, y1, z1, z0, dist
      double precision, allocatable, dimension(:,:)            :: tempo!Total EleMent POint
      integer         , allocatable, dimension(:,:)            :: tlnod!Total List of NODes
      integer                                      ,intent(out):: npoif
      
      integer :: ielem, ipoin, inode
      allocate( tempo(DimPr,npoiw), tlnod(mxelw,mxnow), lpoiw(npoiw))
      
      
      do ipoiw=1,npoiw
        lpoiw(ipoiw)=ipoiw
      end do
     
      z1=0.0
      z0=0.0
      do ipoiw=nnodes+2,npoiw
        x1=coorw(1,ipoiw)
        y1=coorw(2,ipoiw)
        if(DimPr.eq.3) z1=coorw(3,ipoiw)
        do jpoiw=nnodes+1,ipoiw-1
          x0=coorw(1,jpoiw)
          y0=coorw(2,jpoiw)
          if(DimPr.eq.3) z0=coorw(3,jpoiw)
          dist=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
          if(dist.lt.1.0e-7) then
            lpoiw(ipoiw)=-abs(lpoiw(jpoiw))
            do kpoiw=ipoiw+1,npoiw
              lpoiw(kpoiw)=lpoiw(kpoiw)-1
            end do
          end if
        end do
      end do
      
      npoif=0
      do ipoiw=1,npoiw
        if(lpoiw(ipoiw).gt.0) then
          npoif=npoif+1
          do idime=1,DimPr
            tempo(idime,npoif)=coorw(idime,ipoiw)
          end do
        end if
      end do
      
      do ielew=1,nelew
        do inodw=1,nnodw
          tlnod(ielew,inodw)=abs(lpoiw(list_new_nodes(ielew,inodw)))
        end do
      end do
      
      !* Setting the refinment mesh in previous coord and lnods (global)) variables 
      deallocate(coord, lnods)
      allocate(coord(DimPr,npoif), lnods(nelew,nnodw), mesh_conductivity(nelew))
     
      do ielem=1,nelew
        do inode=1,nnodw
          lnods(ielem,inode)=tlnod(ielem,inode)
        end do
      end do
      
      do ielem=1,nelew
        mesh_conductivity(ielem)=mesh_conduc(ielem)
      end do
      
      
      do ipoin=1,npoif
          do idime=1,DimPr
            coord(idime,ipoin)=tempo(idime,ipoin)
          end do
      end do
      
      !write(*,'(a)') 'ELEMENTS'
      !do ielem=1,nelew
      !  write(*,10) ielem,(lnods(ielem,inode),inode=1,nnodw)
      !end do
      !write(*,'(a)') 'END_ELEMENTS'
     
      !write(*,'(a)') 'COORDINATES'
      !do ipoin=1,npoif
      !  write(*,20) ipoin,(coord(idime,ipoin),idime=1,DimPr)
      !end do
      !write(*,'(a)') 'END_COORDINATES'
      ! 
      !10 format(1x,i5,10(1x,i5))
      !20 format(5x,i6,3(2x,f16.9))
      deallocate(tempo, tlnod, lpoiw )
      
    end subroutine checkMesh
    
  !end contains 
    
end module geometry
