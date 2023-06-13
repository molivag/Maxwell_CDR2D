module param
  
  implicit none

  character(len=12) :: File_Nodal_Vals, error_name, coord_name, conec_name, profile_name
  character(len=14) :: testID
  character(len=4)  :: InitElemType, ProbType
  character(len=2)  :: refiType
  integer           :: initnevab, initntotv, nBVs, nBVscol, nband, max_time, simul
  integer           :: upban, lowban, totban, ldAKban !variables defined in GlobalSystem
  integer           :: initElem, initNodes, DimPr, nne, ndofn, totGp, kstab, ktaum, maxband, theta
  integer           :: i_exp, nodalSrc, nodalRec, skipline, postpro, signal, srcType!, srcLoc
  real              :: hnatu, patau
  double precision  :: Cu,lambda, ell, helem, n_val, time_ini, time_fin, u0cond
  double precision, allocatable, dimension(:,:)     :: ngaus, weigp

  double precision, allocatable, dimension(:,:,:,:) :: difma
  double precision, allocatable, dimension(:,:,:)   :: conma
  double precision, allocatable, dimension(:,:)     :: reama !Tensor materials
  double precision, allocatable, dimension(:)       :: force, Icurr !Force and Current vector (respectively)
  integer         , allocatable, dimension(:)       :: srcLoc, recLoc

  contains
    
    subroutine inputData(name_inputFile)
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      !                                                                                   !
      ! subrutina que lee todos los parametros de entrada para la simulacion,             !
      ! la geometria, lista de nodos, coordenadas y parametros de estabilizacion          !
      !                                                                                   !
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      implicit none
      
      character(len=19)             :: name_inputFile
      double precision  :: param_stab1, param_stab2 
      character(len=180)            :: msg
      character(len=*), parameter   :: fileplace = "./"
      integer                       ::  stat, ii
      
      
      open(5, file=fileplace//name_inputFile, status='old', action='read',IOSTAT=stat, IOMSG=msg)
      !open(5, file=fileplace//'inputCDR.dsc',status='old', action='read',IOSTAT=stat, IOMSG=msg)
      
      read(5, 100,iostat=stat,iomsg=msg) &
      InitElemType,ProbType,DimPr,ndofn,totGp,simul,postpro,&
      initElem, initNodes, nne, i_exp, hnatu, refiType,&
      kstab, ktaum, patau, n_val, Cu, ell, lambda,&
      theta, time_ini, time_fin, max_time, u0cond,&
      testID, File_Nodal_Vals, error_name, coord_name, conec_name, profile_name
      if (stat.ne.0) then
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
      helem = 0d0
      param_stab1 = 0.0
      param_stab2 = 0.0

      
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
        if (stat.ne.0)print'(A8,1x,A180)','iomsg = ',msg
        
        read(5,103,iostat=stat,iomsg=msg) &
        difma(1,1,1,2), difma(1,2,1,2), difma(1,3,1,2), &
        difma(2,1,1,2), difma(2,2,1,2), difma(2,3,1,2), &
        difma(3,1,1,2), difma(3,2,1,2), difma(3,3,1,2)
        if (stat.ne.0)print'(A8,1x,A180)','iomsg = ',msg
        
        read(5,103,iostat=stat,iomsg=msg) &
        difma(1,1,2,1), difma(1,2,2,1), difma(1,3,2,1), &
        difma(2,1,2,1), difma(2,2,2,1), difma(2,3,2,1), &
        difma(3,1,2,1), difma(3,2,2,1), difma(3,3,2,1)
        if (stat.ne.0)print'(A8,1x,A180)','iomsg = ',msg
        
        read(5,103,iostat=stat,iomsg=msg) &
        difma(1,1,2,2), difma(1,2,2,2), difma(1,3,2,2), &
        difma(2,1,2,2), difma(2,2,2,2), difma(2,3,2,2), &
        difma(3,1,2,2), difma(3,2,2,2), difma(3,3,2,2)
        if (stat.ne.0)print'(A8,1x,A180)','iomsg = ',msg
        
        read(5,103,iostat=stat,iomsg=msg) &
        conma(1,1,1), conma(1,2,1), conma(1,3,1), &
        conma(2,1,1), conma(2,2,1), conma(2,3,1), &
        conma(3,1,1), conma(3,2,1), conma(3,3,1)
        if (stat.ne.0)print'(A8,1x,A180)','iomsg = ',msg
        
        read(5,103,iostat=stat,iomsg=msg) &
        conma(1,1,2), conma(1,2,2), conma(1,3,2), &
        conma(2,1,2), conma(2,2,2), conma(2,3,2), &
        conma(3,1,2), conma(3,2,2), conma(3,3,2)
        if (stat.ne.0)print'(A8,1x,A180)','iomsg = ',msg
        
        read(5,103,iostat=stat,iomsg=msg) &
        reama(1,1), reama(1,2), reama(1,3), &
        reama(2,1), reama(2,2), reama(2,3), &
        reama(3,1), reama(3,2), reama(3,3)
        if (stat.ne.0)print'(A8,1x,A180)','iomsg = ',msg
        
       read(5,107,iostat=stat,iomsg=msg) force(1), force(2), force(3)
        if (stat.ne.0)print'(A8,1x,A180)','force iomsg = ',msg
        
      end if
      
      if(ndofn.eq.1)then
        read(5,105,iostat=stat,iomsg=msg) Icurr(1)
        if (stat.ne.0)print'(A8,1x,A180)','force iomsg = ',msg
      elseif(ndofn.eq.3)then
        read(5,107,iostat=stat,iomsg=msg) Icurr(1), Icurr(2), Icurr(3)
        if (stat.ne.0)print'(A8,1x,A180)','force iomsg = ',msg
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
      
      !print*, 'Icurr', Icurr(1), Icurr(2), Icurr(3)
      !print*, 'nodalSrc', nodalSrc 
      !print*, 'SrcLoc', srcLoc(1)
      !print*, 'signal', signal 
      !print*, 'nodalRec', nodalRec
      !print*, 'RecLoc', recLoc(1)
      
      close(5)
      !stop 
      
      !print*, srcLoc(1), srcLoc(2), srcType, signal 
      
      
      if(kstab.eq.6)then
        !helem = 2.0**(-(i_exp)) 
        helem = 0.2
        print*,'helem: ', helem
        !Parameter of stabilization in augmented formulation
        param_stab1 = Cu*lambda*(helem**2/ell**2)
        param_stab2 = ell**2 / lambda
        
        difma(1,1,1,1) = difma(1,1,1,1)*param_stab1
        difma(2,2,1,1) = difma(2,2,1,1)*lambda
        difma(3,3,1,1) = difma(3,3,1,1)*param_stab2
        
        difma(1,2,1,2) = difma(1,2,1,2)*param_stab1
        difma(2,1,1,2) = difma(2,1,1,2)*lambda
        
        difma(1,2,2,1) = difma(1,2,2,1)*lambda
        difma(2,1,2,1) = difma(2,1,2,1)*param_stab1
        
        difma(1,1,2,2) = difma(1,1,2,2)*lambda
        difma(2,2,2,2) = difma(2,2,2,2)*param_stab1
        difma(3,3,2,2) = difma(3,3,2,2)*param_stab2
      end if
      
      !Initial elemental and global variables, it will changes if refination is selected.
      initnevab = ndofn*nne
      initntotv = ndofn*initNodes
      
      100 format(7/ 11x, A4,/ ,11x, A4,/, 5(11x,I5,/),         2/,&  !model parameters
      &          4(11x,I7,/), 11x,F7.2,/, 11x,A2,/,            2/,&  !geometry
      &          2(11x,I5,/), 5(11x,F15.5,/),                   2/,&  !stabi
      &          11x,I1,/, 2(11x,f15.5,/), 11x,I7,/,11x,f7.2,/, 2/,&  !time
      &          11x,A14,/, 5(11x,A12,/), / )              !output files
     
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
