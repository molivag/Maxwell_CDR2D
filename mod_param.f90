module param
  
  implicit none

  character(len=10) :: File_PostProcess, error_name, coord_name, conec_name
  character(len=4)  :: InitElemType, ProbType
  character(len=2)  :: refiType, testNo
  integer           :: initnevab, initntotv, nBVs, nBVscol, nband, max_time, simul
  integer           :: upban, lowban, totban, ldAKban !variables defined in GlobalSystem
  integer           :: initElem, initNodes, DimPr, nne, ndofn, totGp, kstab, ktaum, maxband, theta
  integer           :: elemSour, skipline
  real              :: hnatu, patau, time_ini, time_fin, u0cond
  real              :: Cu,mu, ell, helem, i_exp, n_val
  double precision, allocatable, dimension(:,:)     :: ngaus, weigp

  double precision, allocatable, dimension(:,:,:,:) :: difma
  double precision, allocatable, dimension(:,:,:)   :: conma
  double precision, allocatable, dimension(:,:)     :: reama !Tensor materials
  double precision, allocatable, dimension(:)       :: force !Force vector 
  
  integer         , allocatable, dimension(:)     ::  sourceLocation 
  contains
    
    subroutine inputData
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      !                                                                                   !
      ! subrutina que lee todos los parametros de entrada para la simulacion,             !
      ! la geometria, lista de nodos, coordenadas y parametros de estabilizacion          !
      !                                                                                   !
      ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      implicit none
      
      character(len=180)                   :: msg
      character(len=*), parameter          :: fileplace = "./"
      double precision                     :: param_stab1, param_stab2 
      integer                              :: stat, ii, jj, kk, isourelem
      
      open(5, file=fileplace//'inputCDR.dsc',status='old', action='read',IOSTAT=stat, IOMSG=msg)
      
      read(5, 100,iostat=stat,iomsg=msg) InitElemType,ProbType,DimPr,ndofn,totGp,&
      simul,elemSour,skipline,&
      initElem, initNodes, nne, refiType, theta, time_ini, time_fin,max_time, u0cond,&
      kstab, ktaum, patau, hnatu, Cu, mu, ell,i_exp, n_val, testNo, File_PostProcess,&
      error_name, coord_name, conec_name
      
      if (stat.ne.0) then
        print*, ' '
        print*,'!============ STATUS READING INPUT FILE. ============!'
        print "(A38,I4)", "- Status while reading input file is: ", stat
        print*, ' '
        print'(A8,1x,A180)','iomsg = ',msg
      else
        continue
      end if
     
     
      allocate( difma(ndofn,ndofn,DimPr,DimPr) )
      allocate( conma(ndofn,ndofn,DimPr) )
      allocate( reama(ndofn,ndofn) )
      allocate( force(ndofn) )
      allocate( sourceLocation(elemSour) )
      
      
      difma = 0.0
      conma = 0.0
      reama = 0.0
      force = 0.0
      
      helem = 2**(-i_exp)
     
      !Parameter of stabilization in augmented formulation
      param_stab1 = Cu*mu*(helem**2/ell**2)
      param_stab2 = ell**2 / mu
      
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
        !if (stat.ne.0)print'(A8,1x,A180)','force iomsg = ',msg
       
        difma(1,1,1,1) = difma(1,1,1,1)*param_stab1
        difma(2,2,1,1) = difma(2,2,1,1)*mu
        difma(3,3,1,1) = difma(3,3,1,1)*param_stab2
        
        difma(1,2,1,2) = difma(1,2,1,2)*param_stab1
        difma(1,2,1,2) = difma(1,2,1,2)*mu
        
        difma(2,1,2,1) = difma(2,1,2,1)*mu
        difma(2,1,2,1) = difma(2,1,2,1)*param_stab1
        
        difma(1,1,2,2) = difma(1,1,2,2)*mu
        difma(2,2,2,2) = difma(2,2,2,2)*param_stab1
        difma(3,3,2,2) = difma(3,3,2,2)*param_stab2
        
      end if
      
      do ii =1, elemSour
        read(5,*,iostat=stat,iomsg=msg) sourceLocation(ii)
      end do
      
      close(5)
      
      !Initial elemental and global variables, it will changes if refination is selected.
      initnevab = ndofn*nne
      initntotv = ndofn*initNodes
      
      100 format(7/ 11x, A4,/ ,11x, A4,/, 6(11x,I5,/), 2/,&                      !model parameters
      &          3(11x,I7,/), 11x,A2,/      ,2/,&                                !geometry
      &          11x,I1,/, 2(11x,f7.2,/),11x,I7,/,11x,f7.2,/, 2/,&               !time
      &          2(11x,I5,/), 3(11x,F7.2,/), 1(11x,e15.5,/), 3(11x,F7.2,/), 2/,& !stabi
      &          11x,A2,/, 4(11x,A10,/), / )                     !source location
     
      101 format(1/,F12.5,2/)
      102 format(1/,e15.5, e15.5,/, e15.5,e15.5,/)
      103 format(1/,3(e15.5))
      105 format(1/,F12.5,2/)
      106 format(1/,e15.5,e15.5,2/) 
      107 format(1/,e15.5,e15.5,e15.5,2/) 
     
    end subroutine inputData
    
  !end contains
    
end module param
