module param
  
  implicit none

  character(len=12) :: File_Nodal_Vals, error_name, coord_name, conec_name, profile_name, mesh_file
  character(len=14) :: testID
  character(len=4)  :: ProbType, ElemType, initElemType
  character(len=2)  :: refiType
  integer           :: nBVs, nBVscol, nband, t_steps, exacSol, BCsProb
  integer           :: upban, lowban, totban, ldAKban !variables defined in GlobalSystem
  integer           :: DimPr, initnne, nne, ndofn, totGp, kstab, ktaum, maxband, theta, Src_ON
  integer           :: i_exp, nodalSrc, nodalRec, postpro, signal, srcType, srcRHS!, srcLoc
  integer           :: nelem, nnodes, nevab, ntotv, initnevab, initntotv,initNodes, initElem
  real              :: hnatu, patau
  double precision  :: Cu,lambda, ell, helem, n_val, time_ini, time_fin, delta_t
  double precision, allocatable, dimension(:,:)     :: ngaus, weigp

  double precision, allocatable, dimension(:,:)     :: coord !, coordRef
  integer,          allocatable, dimension(:,:)     :: lnods !, lnodsRef
  double precision, allocatable, dimension(:,:,:,:) :: difma
  double precision, allocatable, dimension(:,:,:)   :: conma
  double precision, allocatable, dimension(:,:)     :: reama !Tensor materials
  double precision, allocatable, dimension(:)       :: force, Icurr !Force and Current vector (respectively)
  integer         , allocatable, dimension(:)       :: srcLoc, recLoc

  contains
    
    subroutine inputData(name_inputFile, mesh_file)
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      !                                                                                   !
      ! subrutina que lee todos los parametros de entrada para la simulacion,             !
      ! la geometria, lista de nodos, coordenadas y parametros de estabilizacion          !
      !                                                                                   !
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      implicit none
      
      character(len=19)               :: name_inputFile
      double precision                :: tsteps
      character(len=180)              :: msg
      character(len=*), parameter     :: fileplace = "./"
      integer                         :: stat, ii
      character(len=13), intent(out)  :: mesh_file
      
      
      open(5, file=fileplace//name_inputFile, status='old', action='read',IOSTAT=stat, IOMSG=msg)
      !open(5, file=fileplace//'inputCDR.dsc',status='old', action='read',IOSTAT=stat, IOMSG=msg)
      
      read(5, 100,iostat=stat,iomsg=msg) &
      ProbType,DimPr,ndofn,totGp,exacSol, srcRHS, BCsProb, postpro,&
      mesh_file, initnne, i_exp, hnatu, refiType,&
      kstab, ktaum, patau, n_val, helem, Cu, ell, lambda,&
      theta, time_ini, time_fin, t_steps, Src_ON,&
      testID, File_Nodal_Vals, error_name, coord_name, conec_name, profile_name
      if (stat.ne.0)then
        print*, ' '
        print*,'!============ STATUS READING INPUT FILE. ============!'
        print "(A38,I4)", "- Status while reading input file is: ", stat
        print*, ' '
        print*, ' '
        print'(A8,1x,A180)','iomsg = ',msg
        print*, ' '
        print*, ' '
      else
        continue
      end if
     
      allocate( difma(ndofn,ndofn,DimPr,DimPr) )
      allocate( conma(ndofn,ndofn,DimPr) )
      allocate( reama(ndofn,ndofn) )
      allocate( force(ndofn) )
      allocate( Icurr(ndofn) )
      allocate( srcLoc(nodalSrc) )
      allocate( recLoc(nodalRec) )
      
      difma = 0.0
      conma = 0.0
      reama = 0.0
      force = 0.0
      
      if(ndofn.eq.1)then
        read(5,101) difma(1,1,1,1), difma(1,1,1,2)
        read(5,101) difma(1,1,2,1), difma(1,1,2,2)
       
        read(5,101) conma(1,1,1)
        read(5,101) conma(1,1,2)
        
        read(5,101) reama(1,1)
        
        read(5,105) force(1)
        
      elseif(ndofn.eq.2) then
        read(5,102) difma(1,1,1,1), difma(1,2,1,1), difma(2,1,1,1),difma(2,2,1,1)
        read(5,102) difma(1,1,1,2), difma(1,2,1,2), difma(2,1,1,2),difma(2,2,1,2)
        read(5,102) difma(1,1,2,1), difma(1,2,2,1), difma(2,1,2,1),difma(2,2,2,1)
        read(5,102) difma(1,1,2,2), difma(1,2,2,2), difma(2,1,2,2),difma(2,2,2,2)
        
        read(5,102) conma(1,1,1), conma(1,2,1), conma(2,1,1), conma(2,2,1)
        read(5,102) conma(1,1,2), conma(1,2,2), conma(2,1,2), conma(2,2,2)
        
        read(5,102) reama(1,1), reama(1,2), reama(2,1), reama(2,2)
        
        read(5,106) force(1), force(2)
        
      elseif(ndofn.eq.3)then
        
        read(5,103,iostat=stat,iomsg=msg) &
        difma(1,1,1,1), difma(1,2,1,1), difma(1,3,1,1), &
        difma(2,1,1,1), difma(2,2,1,1), difma(2,3,1,1), &
        difma(3,1,1,1), difma(3,2,1,1), difma(3,3,1,1)
        if(stat.ne.0)then
          print'(A17,I4)','iostat_DIFMA_xx= ',stat
          print'(A13,1x,A180)','iomsg_DIFMA= ',msg
        end if

        read(5,103,iostat=stat,iomsg=msg) &
        difma(1,1,1,2), difma(1,2,1,2), difma(1,3,1,2), &
        difma(2,1,1,2), difma(2,2,1,2), difma(2,3,1,2), &
        difma(3,1,1,2), difma(3,2,1,2), difma(3,3,1,2)
        if(stat.ne.0)then
          print'(A17,I4)','iostat_DIFMA_xy= ',stat
          print'(A13,1x,A180)','iomsg_DIFMA= ',msg
        end if
        
        read(5,103,iostat=stat,iomsg=msg) &
        difma(1,1,2,1), difma(1,2,2,1), difma(1,3,2,1), &
        difma(2,1,2,1), difma(2,2,2,1), difma(2,3,2,1), &
        difma(3,1,2,1), difma(3,2,2,1), difma(3,3,2,1)
        if(stat.ne.0)then
          print'(A17,I4)','iostat_DIFMA_yx= ',stat
          print'(A13,1x,A180)','iomsg_DIFMA= ',msg
        end if
        
        read(5,103,iostat=stat,iomsg=msg) &
        difma(1,1,2,2), difma(1,2,2,2), difma(1,3,2,2), &
        difma(2,1,2,2), difma(2,2,2,2), difma(2,3,2,2), &
        difma(3,1,2,2), difma(3,2,2,2), difma(3,3,2,2)
        if(stat.ne.0)then
          print'(A17,I4)','iostat_DIFMA_yy= ',stat
          print'(A13,1x,A180)','iomsg_DIFMA= ',msg
        end if
        
        read(5,103,iostat=stat,iomsg=msg) &
        conma(1,1,1), conma(1,2,1), conma(1,3,1), &
        conma(2,1,1), conma(2,2,1), conma(2,3,1), &
        conma(3,1,1), conma(3,2,1), conma(3,3,1)
        if(stat.ne.0)then
          print'(A9,I3)','iostat= ',stat
          print'(A8,1x,A180)','iomsg_CONMA= ',msg
        end if
        
        read(5,103,iostat=stat,iomsg=msg) &
        conma(1,1,2), conma(1,2,2), conma(1,3,2), &
        conma(2,1,2), conma(2,2,2), conma(2,3,2), &
        conma(3,1,2), conma(3,2,2), conma(3,3,2)
        if (stat.ne.0)print'(A8,1x,A180)','iomsg_CONMA = ',msg
        
        read(5,103,iostat=stat,iomsg=msg) &
        reama(1,1), reama(1,2), reama(1,3), &
        reama(2,1), reama(2,2), reama(2,3), &
        reama(3,1), reama(3,2), reama(3,3)
        if (stat.ne.0)print'(A8,1x,A180)','iomsg = ',msg
        
       read(5,107,iostat=stat,iomsg=msg) force(1), force(2), force(3)
        if(stat.ne.0)then
          print'(A9,I3)','iostat= ',stat
          print'(A8,1x,A180)','iomsg = ',msg
        end if
        
      end if
      
      if(ndofn.eq.1)then
        read(5,105,iostat=stat,iomsg=msg) Icurr(1)
        if(stat.ne.0)then
          print'(A9,I3)','iostat= ',stat
          print'(A8,1x,A180)','iomsg = ',msg
        end if
      elseif(ndofn.eq.3)then
        read(5,107,iostat=stat,iomsg=msg) Icurr(1), Icurr(2), Icurr(3)
        if(stat.ne.0)then
          print'(A9,I3)','iostat= ',stat
          print'(A8,1x,A180)','iomsg = ',msg
        end if
      else
        write(*,*) 'Source must be 1 or 3 DoF'
      endif
      
      read(5,104,iostat=stat,iomsg=msg) nodalSrc
      do ii =1, nodalSrc
        read(5,*,iostat=stat,iomsg=msg) srcLoc(ii)
      end do
      
      read(5,109,iostat=stat,iomsg=msg) signal
      read(5,109,iostat=stat,iomsg=msg) nodalRec
      do ii =1, nodalRec
        read(5,*,iostat=stat,iomsg=msg) recLoc(ii)
      end do
      
      close(5)
      
      !delta_t  = time_ini
      !time_fin = 20*delta_t
      
      delta_t = (time_fin/ t_steps)
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
      
      100 format(7/ ,11x, A4,/, 7(11x,I5,/),                    2/,&  !model parameters
      &          11x,A13,/, 2(11x,I7,/), 11x,F7.2,/, 11x,A2,/,  2/,&  !geometry
      &          2(11x,I5,/), 3(11x,F10.5,/), 3(11x,F15.5,/),   2/,&  !stabi
      &          11x,I1,/, 2(11x,e15.7,/),2(11x,I5,/),          2/,&  !time
      &          11x,A14,/, 5(11x,A12,/),1/ )              !output files
     
      101 format(1/,F12.5,2/)
      102 format(1/,e15.5, e15.5,/, e15.5,e15.5,/)
      103 format(1/,3(e15.5))
      104 format(1(11x,I2))
      105 format(1/,F12.5,2/)
      106 format(1/,e15.5,e15.5,2/) 
      107 format(1/,e15.5,e15.5,e15.5,2/) 
      108 format(1/,F10.5,F10.5,F10.5,2/) 
     
      109 format(2/,11x,I2)
      110 format(2/,11x,I2,/)
      !108 format(11x,I3,I3,/) 
      
    end subroutine inputData
    
  !end contains
    
end module param
