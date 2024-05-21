module param
  
  implicit none

  integer, parameter :: max_filename_length = 256
  double precision   :: viscocity = 0.5
  double precision   :: lambda = 1.0/1.256637061D-6 
  integer, parameter :: DimPr = 2
  ! character(len=20) :: shape_spec_file
  ! character(len=14) :: testID
  ! character(len=12) :: File_Nodal_Vals, error_name, coord_name, conec_name, time_profile_name, mesh_file
  ! character(len=12) :: File_3DNodal_Vals
  
  character(len=:), allocatable  :: shape_spec_file, testID, File_Nodal_Vals, error_name, coord_name
  character(len=:), allocatable  :: conec_name, time_profile_name, mesh_file, spatial_profile_name,spectrum_FileNames
  character(len=:), allocatable  :: File_3DNodal_Vals, File_Nodal_Vals_ky, PHYSICAL_PROBLEM
  
  
  character(len=4)  :: ProbType, ElemType, initElemType, ky_id, oper
  character(len=2)  :: refiType, splits, TwoHalf, view
  integer           :: nBVs, nBVscol, nband, t_steps, exacSol, BCsProb
  integer           :: upban, lowban, totban, ldAKban  !variables defined in GlobalSystem
  integer           :: initnne, nne, ndofn, totGp, kstab, ktaum, maxband, theta, twindow
  integer           :: nelem, nnodes, nevab, ntotv, initnevab, initntotv, initNodes, initElem, tot_ky, idk_y
  integer           :: i_exp, nodalSrc, nodalRec, postpro, signal, srcType, srcRHS, initCond
  real              :: hnatu, patau
  double precision  :: Cu, ell, helem, n_val, time_ini, time_fin, delta_t
  double precision  :: ky_min, ky_max, y_iFT,  k_y, sigma
  double precision  , allocatable, dimension(:,:)     :: ngaus, weigp

  double precision  , allocatable, dimension(:,:)     :: coord 
  integer           , allocatable, dimension(:,:)     :: lnods
  double precision  , allocatable, dimension(:)       :: Icurr, WaveNumbers !Force and Current vector (respectively)
  character(len=8)  , allocatable, dimension(:)       :: files_ky        
  character(len=4)  , allocatable, dimension(:)       :: nodal_ky 
  character(len=99) , allocatable, dimension(:)       :: PHYSICAL_PROBLEM_OPTIONS 
  integer           , allocatable, dimension(:)       :: receivers, srcLoc, mesh_conductivity

  contains
    
    subroutine inputData(name_inputFile, mesh_file)
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      !                                                                                   !
      ! subrutina que lee todos los parametros de entrada para la simulacion,             !
      ! la geometria, lista de nodos, coordenadas y parametros de estabilizacion          !
      !                                                                                   !
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      implicit none
      
      character(len=*), parameter                   :: fileplace = "./"
      character(len=19)                             :: name_inputFile
      ! character(len=5)                             :: obj
      character(len=180)                            :: msg
      double precision, allocatable, dimension(:)   :: node_found_it
      integer         , allocatable, dimension(:,:) :: recLoc
      integer                                       :: stat, ii, idime, filename_length
      character(len=:),allocatable, intent(out)     :: mesh_file
      
      !Asignar memoria para el nombre del archivo
      allocate(character(max_filename_length) :: shape_spec_file)
      allocate(character(max_filename_length) :: spectrum_FileNames)
      allocate(character(max_filename_length) :: testID)
      allocate(character(max_filename_length) :: File_Nodal_Vals)
      allocate(character(max_filename_length) :: error_name)
      allocate(character(max_filename_length) :: conec_name)
      allocate(character(max_filename_length) :: coord_name)
      allocate(character(max_filename_length) :: time_profile_name)
      allocate(character(max_filename_length) :: spatial_profile_name)
      allocate(character(max_filename_length) :: mesh_file)
      allocate(character(max_filename_length) :: File_3DNodal_Vals)
      allocate(character(max_filename_length) :: PHYSICAL_PROBLEM)
      
      allocate( PHYSICAL_PROBLEM_OPTIONS(6) )
      
      ! DimPr = 2
      open(5, file=fileplace//name_inputFile, status='old', action='read',IOSTAT=stat, IOMSG=msg)
      
      
      !read(5, 100,iostat=stat,iomsg=msg) &
      !ProbType,DimPr,ndofn,totGp,exacSol, srcRHS, BCsProb, postpro, sigma, oper,PHYSICAL_PROBLEM_OPTIONS(1),&
      !PHYSICAL_PROBLEM_OPTIONS(2), PHYSICAL_PROBLEM_OPTIONS(3), PHYSICAL_PROBLEM_OPTIONS(4),&
      !PHYSICAL_PROBLEM_OPTIONS(5), PHYSICAL_PROBLEM_OPTIONS(6),
      !
      !
      do ii =1, 9
        read(5, *) !Salto de lineas iniciales
      end do
      do ii =1, 6
        read(5, *,iostat=stat,iomsg=msg) PHYSICAL_PROBLEM_OPTIONS(ii)
      end do
      do ii =1, 6
        print'(A)', PHYSICAL_PROBLEM_OPTIONS(ii)
      end do
      read(5, 100,iostat=stat,iomsg=msg) &
      ProbType,totGp,exacSol, srcRHS, BCsProb, postpro, sigma
      read(5, 98,iostat=stat,iomsg=msg) &
      mesh_file, view ,initnne, i_exp, hnatu, refiType,&
      kstab, ktaum, patau, n_val, helem, Cu, ell, &
      ky_min, ky_max, tot_ky, splits, y_iFT, shape_spec_file ,File_3DNodal_Vals, &
      ! TwoHalf, ky_min, ky_max, tot_ky, splits, y_iFT, shape_spec_file ,File_3DNodal_Vals, &
      theta, time_ini, time_fin, t_steps, twindow, signal, initCond, & 
      testID, File_Nodal_Vals, error_name, coord_name, conec_name,&
      time_profile_name, spatial_profile_name
      call checkStatus(0,stat,msg)

      ! Iterar sobre las cadenas de caracteres y procesarlas
      do ii = 1, 6
          ! Buscar el carácter '>'
          if (index(PHYSICAL_PROBLEM_OPTIONS(ii), '>') /= 0) then
              ! Extraer la cadena de caracteres después del '>'
              PHYSICAL_PROBLEM = trim(adjustl(PHYSICAL_PROBLEM_OPTIONS(ii)(index(PHYSICAL_PROBLEM_OPTIONS(ii), '>')+1:)))
          else
              ! Código para el caso en que no se encuentre '>'
          end if
      end do
      
      select case (PHYSICAL_PROBLEM)
        case ('Maxwell_In_Non_Convex_Domain')
          ! Código para el primer problema
          oper  = 'MAXW'
          TwoHalf = 'N'
          ndofn = 3
          kstab = 6
          
        case ('Cavity_Driven_Flow')
          ! Código para el segundo problema
          mesh_file = 'Stokes_Flow.msh'
          oper      = 'LAPL'
          TwoHalf   = 'N'
          BCsProb   = 5
          ndofn     = 3
          kstab     = 3     !0(NONE), 1(SUPG), 2(GLS), 3/5(SGS/TG), 4(CG), 6(MVAF)
          
          
        case ('Direct_Current_Electrical_Resistivity_in_2.5-D')
          ! Código para el tercer problema
          oper  = 'LAPL'
          mesh_file = 'resistivity_half_space.msh'
          TwoHalf = 'Y'
          ndofn = 1
          BCsProb = 6
          kstab = 0
          ProbType = 'STAT'
          spectrum_FileNames = 'DC_Spectrum'
          
          
        case ('Electric_Field_Excited_By_A_Double_Line_Source')
          ! Código para el tercer problema
          oper  = 'LAPL'
          mesh_file = 'Double_Line_2_EM.msh'
          ! nodalSrc = 2
          ! srcLoc = (/116, 119 /)
          TwoHalf = 'N'
          ndofn = 1
          BCsProb = 7
          kstab = 0
          ProbType = 'TIME'
          
        case ('Horizontal_Electric_Dipole_in_3-D_A_Wrong_capture_of_solution')
          ! Código para el quinto problema
          oper  = 'MAXW'
          TwoHalf = 'N'
          BCsProb = 2! Revisar pues este caso creo qu eno esta definido depende de si son 8 DoF 
          ndofn = 3
          srcRHS = 0
          kstab = 6
          
        case ('Transient_Electromagnetic_in_2.5-D')
          ! Código para el sexto problema
          ProbType = 'TIME'
          TwoHalf = 'Y'
          BCsProb = 2! Revisar pues este caso creo qu eno esta definido depende de si son 8 DoF 
          oper    = 'MAXW'
          ndofn   = 8
          kstab   = 6
          spectrum_FileNames = 'TEM_Spectrum'
          
          
      end select
      
      
      allocate( Icurr(ndofn) )
      
      if(ndofn.eq.1)then
        read(5,105,iostat=stat,iomsg=msg) Icurr(1)
        if(stat.ne.0)then
          print'(A14,I3)','iostat_Icurr= ',stat
          print'(A8,1x,A180)','iomsg = ',msg
        end if
      elseif(ndofn.eq.3)then
        read(5,107,iostat=stat,iomsg=msg) Icurr(1), Icurr(2), Icurr(3)
        if(stat.ne.0)then
          print'(A14,I3)','iostat_Icurr= ',stat
          print'(A8,1x,A180)','iomsg = ',msg
        end if
      elseif(ndofn.eq.8)then
        read(5,112,iostat=stat,iomsg=msg)&
        Icurr(1), Icurr(2), Icurr(3), Icurr(4), Icurr(5), Icurr(6), Icurr(7), Icurr(8)
        if(stat.ne.0)then
          print'(A14,I3)','iostat_Icurr= ',stat
          print'(A8,1x,A180)','iomsg = ',msg
          print*,'Se detiene en la lectura de Icurr' 
          stop
        endif
      else
        write(*,*) 'Source must be 1,3 or 8 DoF'
      endif
      
      read(5,104,iostat=stat,iomsg=msg) nodalSrc
      allocate( srcLoc(nodalSrc) )
      do ii =1, nodalSrc
        read(5,*,iostat=stat,iomsg=msg) srcLoc(ii)
      end do
      
      ! read(5,109,iostat=stat,iomsg=msg) signal
      read(5,109,iostat=stat,iomsg=msg) nodalRec
      allocate( recLoc(DimPr,nodalRec) )
      !Aqui deberia poner un if para verificar que los nodos tanto de la fuente como del receptor
      !No sobrepasan los nodos maximos disponibles en la malla
      
      !Reading receivers location x,z from input file
      do ii=1,nodalRec !number of total nodes
        read(5,*,iostat=stat,iomsg=msg) (recLoc(idime,ii), idime =1,DimPr )
        call checkStatus(15,stat,msg)
      end do
      close(5)
      
      File_3DNodal_Vals    = trim(File_3DNodal_Vals )
      spectrum_FileNames   = trim(spectrum_FileNames)
      shape_spec_file      = trim(shape_spec_file   )
      File_Nodal_Vals      = trim(File_Nodal_Vals   )
      time_profile_name    = trim(time_profile_name )
      spatial_profile_name = trim(spatial_profile_name )
      error_name           = trim(error_name        )
      conec_name           = trim(conec_name        )
      coord_name           = trim(coord_name        )
      mesh_file            = trim(mesh_file         )
      testID               = trim(testID            )
      
      allocate(receivers(nodalRec), node_found_it(nodalRec))
      !Searching the nearest node to the coordinate pair read
      call SearchingNodes(mesh_file, recLoc, node_found_it)
      !Store the nodes at receivers variable
      receivers = node_found_it
      
      if(TwoHalf == 'Y')then 
        ! Asignar memoria para el vector
        allocate(WaveNumbers(tot_ky))
        filename_length = len(File_Nodal_Vals)
        allocate(character(filename_length) :: File_Nodal_Vals_ky)
        
        ! Generar el vector logarítmico
        do ii = 1, tot_ky 
          WaveNumbers(ii) = 10.0**(log10(ky_min)+(log10(ky_max/ky_min)/real(tot_ky-1))*real(ii-1) )
        end do
        
        ! WaveNumbers=(/1.0D-6, 5.0D-6, 1.0D-5, 5.5D-5, 5.0D-4, 2.5D-3, 5.0D-3,&
        ! &             1.0D-2, 2.5D-2, 5.0D-2, 8.0D-2, 1.0D-1, 2.5D-1, 5.0D-1/)
        ! ky_min = WaveNumbers(1) 
        ! ky_max = WaveNumbers(tot_ky)
        
        nodal_ky=(/'_k01','_k02','_k03','_k04','_k05','_k06','_k07',&
        &          '_k08','_k09','_k10','_k11','_k12','_k13','_k14'/)
        
        files_ky = (/'ky01.dat','ky02.dat','ky03.dat','ky04.dat','ky05.dat',&
        &            'ky06.dat','ky07.dat','ky08.dat','ky09.dat','ky10.dat',&
        &            'ky11.dat','ky12.dat','ky13.dat','ky14.dat'/)

        ! File_3DNodal_Vals = "Final_3D_DC_"
                           ! 3D_Potential
        ! shape_spec_file = "transformed_Voltage.dat"
        ! File_3DNodal_Vals
        File_Nodal_Vals_ky = File_Nodal_Vals
      else
        tot_ky = 1
        k_y    = 1.0
      endif
      != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

      !delta_t  = time_ini
      !time_fin = 20*delta_t
      
      delta_t = (time_fin/ t_steps)
      print*,signal
      print*,initCond
      
      
      
      
      ! delta_t = (time_fin - time_ini)/t_steps
      
      !print*, 'delta_t', delta_t
      !t_steps = floor(tsteps) !redondeo al numero inmediato superior 
      !print*, 'time stpes ', t_steps
      
      !if(kstab.eq.6)then
      !  !helem = 2.0**(-(i_exp)) 
      !  print*,'Cu    : ', Cu
      !  print*,'lambda: ', lambda
      !  print*,'helem : ', helem
      !  print*,'ℓ     : ', ell
      !  !Parameter of stabilization in augmented formulation
      !  print*,' '  
      !  Su = Cu*lambda*(helem**2/ell**2)
      !  Sp = ell**2 / lambda
      !  
      !  print*,'Su : ', Su
      !  print*,'Sp : ', SP
      !  
      !  difma(1,1,1,1) = difma(1,1,1,1)*Su
      !  difma(1,2,1,2) = difma(1,2,1,2)*Su
      !  difma(2,1,2,1) = difma(2,1,2,1)*Su
      !  difma(2,2,2,2) = difma(2,2,2,2)*Su
      !  
      !  difma(3,3,1,1) = difma(3,3,1,1)*Sp
      !  difma(3,3,2,2) = difma(3,3,2,2)*Sp
      !  
      !  difma(2,2,1,1) = difma(2,2,1,1)*lambda
      !  difma(2,1,1,2) = difma(2,1,1,2)*lambda
      !  difma(1,2,2,1) = difma(1,2,2,1)*lambda
      !  difma(1,1,2,2) = difma(1,1,2,2)*lambda
      !end if
      
      !Initial elemental and global variables, it will changes if refination is selected.
      
      !100 format(7/ ,11x, A4,/, 7(11x,I5,/), 11x,F15.7,/, 11x,A4,/,                    3/,&  !model parameters
      !&          11x,A,/, 11x,A2,/, 2(11x,I7,/), 11x,F7.2,/, 11x,A2,/,                 2/,&  !geometry
      !&          2(11x,I5,/), 3(11x,F10.5,/), 3(11x,F15.5,/),                          2/,&  !stabi
      !&          11x,A,/, 2(11x,e12.5,/), 11x,I3,/,11x,A,/, 11x,F10.3,/, 2(/,11x,A,/), 2/,&  !Fourier Transform
      !&          11x,I1,/, 2(11x,e14.7,/), 4(11x,I5,/),                                2/,&  !time
      !&
      !
      !
      100 format(5/ ,11x, A4,/, 5(11x,I5,/), 11x,F15.7,/,                           2/)   !model parameters
      98 format(11x,A,/, 11x,A2,/, 2(11x,I7,/), 11x,F7.2,/, 11x,A2,/,            2/,&  !geometry
      &          2(11x,I5,/), 3(11x,F10.5,/), 2(11x,F15.5,/),                       2/,&  !stabi
      &          2(11x,e12.5,/), 11x,I3,/,11x,A,/, 11x,F10.3,2/, 2(11x,A,/),        2/,&  !Fourier Transform
      &          11x,I1,/, 2(11x,e14.7,/), 4(11x,I5,/),                             2/,&  !time
      &          11x,A,/, 6(11x,A,/),                                               1/ )  !output files

      !**Format for 1Dof tensors
      101 format(1/,F12.5,7/)
      !**Format for 2Dof tensors
      102 format(1/,e15.5, e15.5,/, e15.5,e15.5,/)
      122 format(2(e15.5),5/)
      !**Format for 3Dof tensors
      103 format(1/,3(e15.5))
      133 format(3(e15.5),5/)
      !**Format for 8Dof tensors
      111 format(1/,8(e15.5)) 
      112 format(1/,8(f15.5),2/) 
      
      104 format(1(11x,I4))
      105 format(1/,F12.5,2/)
      106 format(1/,e15.5,e15.5,2/) 
      107 format(1/,3(f15.5),2/) 
      108 format(1/,F10.5,F10.5,F10.5,2/) 
      109 format(2/,11x,I2)
      110 format(2/,11x,I2,/)
      !108 format(11x,I3,I3,/) 
      
      
      ! 103 format(1/,3(e15.5))
      ! 107 format(1/,3(f15.5),2/) 
      
      
    end subroutine inputData
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    ! 
    subroutine SearchingNodes(file_mesh, recLoc, node_found_it)
      ! subroutine SearchingNodes(recLoc, node_found_it)
      implicit none

      character(len=*), parameter :: fileplace = "Msh/"
      integer         , parameter :: max_nodos = 10000
      character(len=:), allocatable, intent(in) :: file_mesh
      integer, dimension(DimPr,nodalRec), intent(in) :: recLoc
      character(len=180)              :: msg
      double precision                :: coord_x(max_nodos), coord_y(max_nodos)
      double precision                :: x, dummy, y
      double precision                :: distancia_minima, distancia_actual
      integer                         :: nodo(max_nodos)
      integer                         :: nodo_mas_cercano
      integer                         :: i, ireceiver, num_nodos, stat
      double precision, dimension(nodalRec), intent(out) :: node_found_it

      ! Leer el archivo de coordenadas (coor.dat)
      open(1, file=fileplace//file_mesh, status='old', action='read',IOSTAT=stat, IOMSG=msg)
      ! open(unit=1, file='coor.dat', status='old', action='read')
      num_nodos = 0
      read(1,*)
      read(1,*) nnodes
      if(view.eq.'xy')then
        do i = 1, nnodes 
          read(1, *, iostat=stat, iomsg=msg) nodo(i), coord_x(i), coord_y(i), dummy  
          num_nodos = i
          call checkStatus(6,stat,msg)
        end do
      elseif(view.eq.'xz')then
        do i = 1, nnodes 
          read(1, *, iostat=stat, iomsg=msg) nodo(i), coord_x(i), dummy , coord_y(i)
          num_nodos = i
          call checkStatus(6,stat,msg)
        end do


      endif
      close(1)

      loop_receiver: do ireceiver =1, nodalRec
        ! Leer las coordenadas a buscar
        ! write(*,*) "Ingrese la coordenada x:"
        ! read(*,*) x = recLoc(1, ireceiver)
        x = recLoc(1, ireceiver)
        ! write(*,*) "Ingrese la coordenada y:"
        ! read(*,*) y = recLoc(2, ireceiver)
        y = recLoc(2, ireceiver)

        ! Inicializar la distancia mínima
        distancia_minima = sqrt((x - coord_x(1))**2 + (y - coord_y(1))**2)
        nodo_mas_cercano = nodo(1)
        ! Buscar la coordenada más cercana
        do i = 2, num_nodos
           distancia_actual = sqrt((x - coord_x(i))**2 + (y - coord_y(i))**2)
           if (distancia_actual < distancia_minima) then
              distancia_minima = distancia_actual
              nodo_mas_cercano = nodo(i)
           end if
        end do
        node_found_it(ireceiver) = nodo_mas_cercano
        ! write(*,125) "El nodo: ", nodo_mas_cercano, " es el más cercano a (", x, " , ", y, ")"
      end do loop_receiver
      
      ! 125 format(A26,f10.1,A3,f3.1,A14,I0)
      125 format(A9,I5,A24,f6.1,A3,f3.1,A1)
      ! Mostrar el resultado
      
    end subroutine SearchingNodes
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    ! 
    !
    subroutine checkStatus(id_reading, stat, msg)
      implicit none
      character(len=150) , intent(in)  :: msg
      integer, intent(in)              :: stat, id_reading
      ! character(len=5)                 :: object
      
      select case(id_reading)
        case(0)
        if (stat.ne.0)then
          print*, ' '
          print*,'!============ STATUS READING INPUT FILE. ============!'
          print "(A38,I0)", "- Status while reading input file is: ", stat
          ! print*, ' '
          print'(A9,A)','- iomsg: ',msg
          print*, ' >>>>> if status = 64 = Preconnected file comprises unformatted records'
          print*, ' >>>>> mean there is a <tab> instead of <BS> in input file'
          print*, ' >>>>> program stopped <<<<'
          stop
        else
          continue
        end if
        case(1)
          if(stat.ne.0)then
            print'(A13,1x,A180)','In input file error while reading tensor properties= ',msg
            print*, ' >>>>> program stopped <<<<'
          end if
        case(2)
          if(stat.ne.0)then
            print'(A13,1x,A)','In input file error while reading tensor properties= ',msg
            print*, ' >>>>> program stopped <<<<'
          end if
        case(3)
          if(stat.ne.0)then
            print'(A13,1x,A)','In input file error while reading tensor properties= ',msg
            print*, ' >>>>> program stopped <<<<'
            stop
          end if
        case(4)
          if ( stat /= 0 )then
            print*, ' ' 
            print*, 'error in read mesh file first line, module geometry' 
            print'(A9,I3)','iostat= ',stat
            print'(A8,1x,A)','iomsg= ',msg
            print'(A55,A)', 'error >>>> Something wrong during reading of mesh file ', mesh_file
            print*, ' ' 
            stop
          end if
          
        case(5)
          if( stat /= 0 )then
            print*, ' ' 
            print*, 'error while reading init nodes in mesh file at, module geometry' 
            print'(A9,I3)','iostat= ',stat
            print'(A8,1x,A180)','iomsg= ',msg
            print'(A55,A)', 'error >>>> Something wrong during reading of mesh file ',mesh_file
            print*, ' ' 
            stop
          end if
        case(6)
          if(stat.ne.0)then
            print*, ' '
            print "(A46,I4)", " -Status while reading receiver locations is: ", stat
            print'(A74,1x,A)',' >>>>>>>Error in reading receivers coordinate location in SearchingNodes = ',msg
            print*, ' '
            stop
          end if
        case(7)
          IF( stat /= 0 )then
            print*,'error while reading mesh file in geometry module cord3D, iostat= ',stat
            print*, msg
          end if
        case(8)
          IF( stat /= 0 )then
            print*,'error while reading initELem in geometry module iostat= ',stat
            print*, msg
          END IF
        case(9)
          IF( stat /= 0 )then
            print*,'error while reading lnodes in geometry module iostat= ',stat
            print*, msg
          END IF
        case(10)
        if (stat.ne.0)then
          print*, ' '
          print "(A9,I2,A16,I1,A)", "- Status ", stat, " while reading ",ndofn," DoF in tensor properties"
          ! print*, ' '
          print'(A9,A)','- iomsg: ',msg
          ! print*, ' >>>>> if status = 64 = Preconnected file comprises unformatted records'
          ! print*, ' >>>>> mean there is a <tab> instead of <BS> in input file'
          ! print*, ' >>>>> program stopped <<<<'
          ! stop
        else
          continue
        end if
        case(15)
          IF( stat /= 0 )then
            print*,'iostat= ',stat
            print'(A8,1x,A180)','iomsg = ',msg
            print*, 'In recLoc there is an error in reading input file'
            stop
          end if
      end select
      
    end subroutine checkStatus
    
    
  !end contains
    
end module param
