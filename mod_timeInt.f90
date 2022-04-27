module timeInt
  use param
  use library!, only: ApplyBVs, GlobalSystem_Time, file_name_inc, GID_PostProcess, MKLsolverResult, MKLfactoResult

  contains

    ! subroutine initialCondition(A_K, A_F, time_ini, time_fin, max_time, u0_cond, U_init) 
      
    !   !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    !   ! Routine to proyect the initial condition into the FE space:
    !   ! 
    !   !           u(0) = u^0         ; continuous form of the initial condition
    !   !      (uh(0),v) = (u^0,v)     ; variational form
    !   !  di(0)*(Ni,Nj) = (u^0,Nj)   
    !   !         C*d(0) = U^0         ; matrix form
    !   ! where:
    !   !   C = (c_ij) ; c_ij = ∫(Ni*Nj)dΩ ; capacity matrix
    !   ! U^0 = (Ui^0) ; Ui^0 = ∫(u0*Nj)dΩ ; initial condition vector
    !   ! d(0)= di(0)                      ; dof at time = 0
    !   !     
    !   !     from: Johnson, C. (2009). Numerical Sol of PDE by the FEM. Dover, 149-150
    !   !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      
    !   implicit none
      
    !   double precision, allocatable, dimension(:,:) :: A_C0, A_U0, U_init
    !   integer, intent(in)                           :: time_ini, time_fin
    !   integer                                       :: ielem, igaus, ibase, Euler_method, time, max_time
    !   double precision                              :: delta_t
    !   real                                          :: u0_cond, time_ini, time_fin
    !   !=============== S O L V E R ===============              
    !   external                             :: dgetrf, dgetrs
    !   integer, allocatable, dimension(:)   :: S_ipiv
    !   character(len=1)                     :: S_trans
    !   integer                              :: S_m, S_n, S_lda, S_info, S_nrhs
    !   !SE DEBERIA CAMBIAR LAS VARIABLES DE TIEMPO A VARIABLES  LOCALES PONIENDOLAS COMO SALIDA DUMMY VARIABLE
    !   ! EN EL CALL DE INPUT DATA PUES OSLO SE USAN EN DOS SUBRUTINAS Y NO TIENE CASO QUE ESTEN EN TODO EL CODIGO
      
    !   Euler_method = 2
    !   select case(Euler_method)
    !     case(1)
    !       !Delta t for implicit time
    !       delta_t  = 0.001 
         
    !     case(2)
    !       !Delta t for implicit time
    !       delta_t  = ( time_fin - time_ini ) / (max_time + 1.0)   !Step size
    !     case default
    !       print*, 'Not Euler Method defined, must be 1 or 2. Current value ', Euler_type
    !   end select
      
    !   write(*,*) ' '
    !   print*,'!============ TIME DISCRETIZATION ============!'
    !   write(*,"(A19,4X,F10.3,1X,A10)") ' - Initial time:          ', time_ini,' '
    !   write(*,"(A19,4X,F10.3,1X,A10)") ' - Final time:            ', time_fin,' '
    !   write(*,"(A19,8X,I3,1X,A10)")    ' - Number of steps:       ', max_time,' '
    !   write(*,"(A21,7X,F10.6,1X,A10)") ' - Step size(∆t):         ', delta_t,' '
    !   write(*,*) ' '
      
      
    !   call GID_PostProcess(Uprev, File_PostMsh, 'msh', time, 0.0, time_fin)
      
    !   allocate( U_init(ntotv), A_U0(ntotv), A_C0(ntotv,ntotv) )
      
    !   do ielem = 1, nelem 
    !     !gather
    !     u_init = 0.0    !Fe(nevab)
    !     call SetElementNodes(ielem, element_nodes, nodeIDmap)
    !     !do-loop: compute element capacity and stiffness matrix Ke Ce and element vector Fe
    !     do igaus = 1, TotGp
    !       !Jaco = J2D(element_nodes, dN_dxi, dN_deta, igaus)
    !       !detJ = m22det(Jaco)
    !       !Jinv = inv2x2(Jaco)
    !       !dvol = detJ *  weigp(igaus,1)
    !       !hmaxi = elemSize(Jinv)
    !       do ibase = 1, nne
    !         basis(ibase) = N(ibase,igaus)
    !       end do
    !       call Galerkin_Init_Cond(dvol, basis, u0_cond, C0e, u0e) 
    !       !call TauMat(hmaxi,tauma)
    !       !call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, Ke, Fe)
    !     end do
        
    !     call Assemb_Glob_Mat(nodeIDmap, C0e, A_C0)      !Assemble Global Capacity Matrix C 
    !     call Assemb_Glob_Mat(nodeIDmap, u0e, A_U0)      !Assemble Global Capacity Matrix C 
        
    !   end do
      
    !   print*,'!=============== SOLVER (LAPACK) ===============!'
    !   S_m   = size(A_C0,1)
    !   S_n   = size(A_C0,2)
    !   S_lda = max(1,size(A_C0,1)) ! lda ≥ max(1, n).
    !   S_trans = 'N'
    !   S_nrhs  = 1
    !   S_ldb   = max(1,size(A_U0,1))
      
    !   allocate( S_ipiv(max(1,min(S_m, S_n)) ) )
      
    !   print*,'  •INITIALIZING LU FACTORIZATION A = P*L*U.....'
    !   call dgetrf( S_m, S_n, A_C0, S_lda, S_ipiv, S_infoLU )
    !   call MKLfactoResult( S_infoLU )
      
    !   U_ini = A_U0
    !   print*, ' '  
    !   print*,'  •SOLVING SYSTEM OF EQUATIONS..... '
    !   call dgetrs( S_trans, S_n, S_nrhs, A_C0, S_lda, S_ipiv, U_ini, S_ldSol, S_infoSOL )
    !   call MKLsolverResult( S_infoSOL )
      
      
      
    ! end subroutine initialCondition
    
    
    !subroutine ForwardEuler(N, dN_dxi, dN_deta, Hesxieta, time_ini, time_fin, max_time, u0_cond,&
    !  &                      nofix, ifpre, presc, S_m, S_n, S_trans, S_nrhs, S_ipiv, S_ldSol, workdim)
    !  
    !  implicit none
    !  
    !  external :: dgbtrf, dgbtrs, dgbrfs
    !  
    !  double precision, dimension(nne,TotGp), intent(in)     :: N, dN_dxi, dN_deta
    !  double precision, dimension(3,nne), intent(in)         :: Hesxieta
    !  double precision, dimension(ndofn,nBVs), intent(in)    :: presc
    !  integer, dimension(ndofn,nBVs), intent(in)             :: ifpre
    !  integer, dimension(nBVs), intent(in)                   :: nofix
    !  character(len=1), intent(in)                           :: S_trans
    !  integer,intent (in)                                    :: S_m, S_n, S_nrhs, S_ldSol, max_time
    !  real, intent(in)                                       :: time_ini, time_fin, u0_cond
    !  ! - - Local Variables - -!
    !  double precision, allocatable, dimension(:,:) :: A_K, A_C, A_F
    !  double precision, allocatable, dimension(:,:) :: AK_time, AK_LU, Uprev, Unext, rhs_time, u_init
    !  double precision, allocatable, dimension(:)   :: S_ferr, S_berr, S_work
    !  integer, allocatable, dimension(:)            :: S_ipiv, S_iwork
    !  !character(len=20)                             :: solution_file_name
    !  double precision :: delta_t
    !  integer          :: time, info, workdim, ans, i, j
    !  real             :: nt 
    !  
    !  
    !  allocate( AK_time(ldAKban,ntotv), AK_LU(ldAKban,ntotv))
    !  allocate( rhs_time(ntotv,1), u_init(ntotv,1) )
    !  allocate( Uprev(S_ldSol, 1))
    !  allocate( Unext(S_ldSol, 1))
    !  allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))
    !  
    !  allocate( S_work(workdim), S_iwork(S_ldSol), S_ferr(S_nrhs), S_berr(S_nrhs) )
    !  
    !  
    !  Uprev    = 0.0
    !  Unext    = 0.0
    !  time     = 0                                            !initializing the time
    !  delta_t  =  
    !  write(*,"(A21,7X,F10.6,1X,A10)") ' - Step size(∆t):         ', delta_t,' '
    !  write(*,*) ' '
    !  
    !  call GlobalSystem(N, dN_dxi, dN_deta, Hesxieta, A_K, A_F)
    !  call ApplyBVs(nofix,ifpre,presc,A_K,A_F)
    !  
    !  DO J=1,ntotv
    !    DO I=MAX(J-upban,1),MIN(J+lowban, ntotv)
    !      iband = I-J+totban
    !      A_KB(iband,J)=A_K(I,J)
    !    end do
    !  end do
    !  
    !  call initial_condition(A_K, A_F, S_ldSol, u0_cond, u_init) 
    !  Uprev  = u_init                                   !u in present time 
    !  
    !  write(*,*) ' '
    !  print'(A11,I3,A3,F8.3,A)',' time step:',time,'  = ',time_ini,' is the value of u by the initial condiction'
    !  call GID_PostProcess(Uprev, File_PostRes, 'res', time, 0.0, time_fin)
    !  print*, 'Starting time integration. . . . .'
    !  write(*,*) ' '
    !  
    !  nt = 0.0 
    !  do time = 1, max_time +1
    !    nt = nt + delta_t
    !   
    !    !-------- Explicit Scheme --------!
    !                       
    !    1.0/delta_t * A_C*Unext = A_F + ( 1.0/delta_t*A_C - A_K)* Uprev
    !    
    !    AK_time  = (1.0/delta_t)*A_C + A_K !A_C + delta_t*A_K
    !    rhs_time = A_F 
    !    
    !    Uprev = Unext
    !    
    !    !---------- Printing and writing results -----------!
    !    !call file_name_inc(solution_file_name)
    !    print'(A11,I3,A3,F8.3,A5,F8.3,A5)',' time step:',time,' =',nt,'   of',time_fin,' seg'
    !    call GID_PostProcess(Uprev, File_PostRes, 'res', time, nt, time_fin)
    !    
    !    
    !  end do
    !  
    !  print*, ' '
    !  print*, 'Shape of Global K: ',shape(AK_time)
    !  print*, 'Shape of Global F: ',shape(rhs_time)
    !  print*, 'Shape of Solution: ',shape(Uprev)
    !  write(*,*)
    !  DEALLOCATE( AK_time, rhs_time, Uprev, Unext)
    !  
    !end subroutine ForwardEuler
    
    
    
    subroutine BackwardEuler(N, dN_dxi, dN_deta, Hesxieta, time_ini, time_fin, max_time, u0_cond,&
      &                      nofix, ifpre, presc, S_m, S_n, S_trans, S_nrhs, S_ipiv, S_ldSol, workdim)
      
      implicit none
      
      external :: dgbtrf, dgbtrs, dgbrfs
      
      double precision, dimension(nne,TotGp), intent(in)     :: N, dN_dxi, dN_deta
      double precision, dimension(3,nne), intent(in)         :: Hesxieta
      double precision, dimension(ndofn,nBVs), intent(in)    :: presc
      integer, dimension(ndofn,nBVs), intent(in)             :: ifpre
      integer, dimension(nBVs), intent(in)                   :: nofix
      character(len=1), intent(in)                           :: S_trans
      integer,intent (in)                                    :: S_m, S_n, S_nrhs, S_ldSol, max_time
      real, intent(in)                                       :: time_ini, time_fin, u0_cond
      ! - - Local Variables - -!
      double precision, allocatable, dimension(:,:) :: A_K, A_C, A_F, AK_2ord, rhs_2ord
      double precision, allocatable, dimension(:,:) :: AK_time, Uprev, Ucurr, Unext, rhs_time, Uinit, F_plus_MU
      double precision, allocatable, dimension(:)   :: S_ferr, S_berr, S_work
      integer, allocatable, dimension(:)            :: S_ipiv, S_iwork
      integer                                       :: time, info, workdim
      double precision                              :: delta_t
      real                                          :: nt 
      
      
      allocate( AK_time(ldAKban,ntotv), AK_2ord(ldAKban,ntotv))
      allocate( rhs_time(ntotv,1), rhs_2ord(ntotv,1) )
      allocate( Uinit(ntotv,1), Uprev(S_ldSol, 1), Ucurr(S_ldSol, 1), Unext(S_ldSol, 1))
      allocate( F_plus_MU(S_ldSol, 1))
      allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))
      allocate( S_work(workdim), S_iwork(S_ldSol), S_ferr(S_nrhs), S_berr(S_nrhs) )
      
      !prev   ;       curr    ;       futu
      Ucurr    = 0.0
      Unext    = 0.0
      rhs_time = 0.0
      time     = 0                                            !initializing the time
      call GID_PostProcess(Ucurr, File_PostMsh, 'msh', time, 0.0, time_fin)
      
      !call initial_condition(0_cond, u_init) 
      Uinit = u0_cond                                  !Para prueba lo dejo sin subroutina
      Uprev = Uinit                                   !u in present time 
      
      write(*,*) ' '
      print'(A11,I3,A3,F8.3,A)',' time step:',time,'  = ',time_ini,' is the value of u by the initial condiction'
      call GID_PostProcess(Uprev, File_PostRes, 'res', time, 0.0, time_fin)
      print*, 'Starting time integration. . . . .'
      write(*,*) ' '
      
      nt = 0.0 
      do time = 1, max_time +1
        nt = nt + delta_t!,time_fin,delta_t
       
        call GlobalSystem_Time(N, dN_dxi, dN_deta, Hesxieta, delta_t, Uprev, A_K, A_C, A_F)
        !call GlobalSystem_Time(N, dN_dxi, dN_deta, Hesxieta, delta_t, Ucurr, glbCondMa, glbCapaMa, glbSrcVec)
        
        !-------- Implicit Scheme --------!
        
        !---1st-Order Backward Euler   
        AK_time  = (1.0/delta_t)*A_C + A_K !A_C + delta_t*A_K
        rhs_time = A_F 
       
        call ApplyBVs(nofix,ifpre,presc,AK_time,rhs_time)
        
        !------------- Solver -------------!
        Ucurr = rhs_time   !here mkl will rewrite Unext by the solution vector  
        
        !--- Factorizing Matrix
        call dgbtrf(S_m, S_n, lowban, upban, AK_time, ldAKban, S_ipiv, info)
        if(info.ne.0)then
          call MKLfactoResult('dgbtrf',info) 
          print'(A32,I3)', '<<<Error in factorization at time: ', time
        endif
        !--- Solving System of Eqns
        call dgbtrs(S_trans, S_n, lowban, upban, S_nrhs, AK_time, ldAKban, S_ipiv, Ucurr, S_ldSol, info )
        
        if(info.ne.0)then
          call MKLsolverResult('dgbtrs',info) 
          print'(A32,I3)', '<<<Error in solving system of equation at time: ', time
        endif
        
        !glbCapaMa - - - M
        !glbCondMa - - - K
        !glbSrcVec - - - F
        
        ! - - - 2nd-Order Backward Euler
        AK_2ord  = (3/2*delta_t)*A_C + A_K 
        ! Ucurr(Un) --> computed throug 1st order Euler
        ! Uprev(Un-1) --> computed throug the initial cond for time 1   
        !M*(Un - Un-1)
        call TimeLevels(N, dN_dxi, dN_deta, Hesxieta, delta_t, Ucurr, Uprev, F_plus_MU)
        !By the last subroutine we compute the product of right-hand side of 2nd order backward Euler 
        !done in elemental form to match the dimension btwn band matrix and global vecotr
        rhs_2ord = F_plus_MU 
        
        call ApplyBVs(nofix,ifpre,presc, AK_2ord, rhs_2ord)
        !------------- Solver for future time -------------!
        Unext = rhs_2ord   !here mkl will rewrite Unext by the solution vector  
        
        !--- Factorizing Matrix
        call dgbtrf(S_m, S_n, lowban, upban, AK_2ord, ldAKban, S_ipiv, info)
        if(info.ne.0)then
          call MKLfactoResult('dgbtrf',info) 
          print'(A32,I3)', '<<<Error in factorization at time: ', time
        endif
        !--- Solving System of Eqns
        call dgbtrs(S_trans, S_n, lowban, upban, S_nrhs, AK_2ord, ldAKban, S_ipiv, Unext, S_ldSol, info )
        if(info.ne.0)then
          call MKLsolverResult('dgbtrs',info) 
          print'(A32,I3)', '<<<Error in solving system of equation at time: ', time
        endif          
        !---------- Printing and writing results -----------!
        print'(A11,I3,A3,F8.3,A5,F8.3,A5)',' time step:',time,' =',nt,'   of',time_fin,' seg'
        call GID_PostProcess(Unext, File_PostRes, 'res', time, nt, time_fin)

        Uprev = Unext
       
      end do
      
      print*, ' '
      print*, 'Shape of Global K: ',shape(AK_time)
      print*, 'Shape of Global F: ',shape(rhs_time)
      print*, 'Shape of Solution: ',shape(Unext)
      write(*,*)
      DEALLOCATE( AK_time, rhs_time, Uprev, Unext)
      
    end subroutine BackwardEuler
    
    
    
  !end contains

end module timeInt

!Estas lineas son para revisar la solucion e implementar el refinamiento de la solucion,
!Deben agregarse despues de Uprev = u_fut
        ! print*, ' '
        ! print*,"Solution vector"
        ! do i =1, ntotv
        !   print '(F13.5)', u_fut(i,1)
        ! end do
        ! !----- Refining Solution -----!
        ! call dgbrfs(S_trans, S_n, lowban, upban, S_nrhs, A_K, ldAKban, AK_LU, ldAKban, S_ipiv,&
        !&            rhs_time, S_ldSol, u_fut, S_ldSol, S_ferr, S_berr, S_work, S_iwork, info )
        ! if(info.ne.0)then
        !   call MKLsolverResult('dgbrfs',info) 
        !   print'(A32,I3)', '<<<The parameter ', info, ' had an illegal value in time: ', time
        ! endif
        ! print*, ' '
        ! print*, 'Refined Solution'
        ! do i =1, ntotv
        !   print '(9F13.5)', u_fut(i,1)
        ! end do
        ! 
        ! print*, ' '
        ! print*, 'Avanzar? 1=yes ; other stop'
        ! read(*,*) ans
        ! if(ans.eq.1)then
        !   continue
        ! else
        !   stop
        ! endif
        
