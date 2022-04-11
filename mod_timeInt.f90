module timeInt
  use param
  use library!, only: ApplyBVs, GlobalSystem_Time, file_name_inc, GID_PostProcess, MKLsolverResult, MKLfactoResult

  contains

    subroutine initialCondition( )
      
      !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      ! Routine to proyect the initial condition into the FE space:
      ! 
      !           u(0) = u^0         ; continuous form of the initial condition
      !      (uh(0),v) = (u^0,v)     ; variational form
      !  di(0)*(Ni,Nj) = (u^0,Nj)   
      !         C*d(0) = U^0         ; matrix form
      ! where:
      !   C = (c_ij) ; c_ij = ∫(Ni*Nj)dΩ ; capacity matrix
      ! U^0 = (Ui^0) ; Ui^0 = ∫(u0*Nj)dΩ ; initial condition vector
      ! d(0)= di(0)                      ; dof at time = 0
      !     
      !     from: Johnson, C. (2009). Numerical Sol of PDE by the FEM. Dover, 149-150
      !* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      
      implicit none
      
      
      
      
      
      
      
      
      
    end subroutine initialCondition
    
    
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
      double precision, allocatable, dimension(:,:) :: A_K, A_C, A_F
      double precision, allocatable, dimension(:,:) :: AK_time, AK_LU, u_pre, u_fut, rhs_time, u_init
      double precision, allocatable, dimension(:)   :: S_ferr, S_berr, S_work
      integer, allocatable, dimension(:)            :: S_ipiv, S_iwork
      !character(len=20)                             :: solution_file_name
      double precision :: delta_t
      integer          :: time, info, workdim, ans, i, j
      real             :: nt 
      
      
      allocate( AK_time(ldAKban,ntotv), AK_LU(ldAKban,ntotv))
      allocate( rhs_time(ntotv,1), u_init(ntotv,1) )
      allocate( u_pre(S_ldSol, 1))
      allocate( u_fut(S_ldSol, 1))
      allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))
      
      allocate( S_work(workdim), S_iwork(S_ldSol), S_ferr(S_nrhs), S_berr(S_nrhs) )
      
      
      !S_ipiv   = 0
      u_pre    = 0.0
      u_fut    = 0.0
      rhs_time = 0.0
      AK_time  = 0.0
      time     = 0                                            !initializing the time
      delta_t  = ( time_fin - time_ini ) / (max_time + 1.0)   !Step size
      
      !solution_file_name = 'CDR_u0000.post.res'
      write(*,*) ' '
      print*,'!============ TIME DISCRETIZATION ============!'
      write(*,"(A19,4X,F10.3,1X,A10)") ' - Initial time:          ', time_ini,' '
      write(*,"(A19,4X,F10.3,1X,A10)") ' - Final time:            ', time_fin,' '
      write(*,"(A19,8X,I3,1X,A10)")    ' - Number of steps:       ', max_time,' '
      write(*,"(A21,7X,F10.6,1X,A10)") ' - Step size(∆t):         ', delta_t,' '
      write(*,*) ' '
      call GID_PostProcess(u_pre, File_PostMsh, 'msh', time, 0.0, time_fin)
      
      !call initial_condition(node_num, node_xy, u0_cond, u_init) 
      u_init = u0_cond                                  !Para prueba lo dejo sin subroutina
      u_pre  = u_init                                   !u in present time 
      
      write(*,*) ' '
      print'(A11,I3,A3,F8.3,A)',' time step:',time,'  = ',time_ini,' is the value of u by the initial condiction'
      call GID_PostProcess(u_pre, File_PostRes, 'res', time, 0.0, time_fin)
      print*, 'Starting time integration. . . . .'
      write(*,*) ' '
      nt = 0.0 
      do time = 1, max_time +1
        nt = nt + delta_t!,time_fin,delta_t
        
        call GlobalSystem_Time(N, dN_dxi, dN_deta, Hesxieta, S_ldsol, delta_t, u_pre, A_K, A_C, A_F)
        !-------- Implicit Scheme --------!
        AK_time  = (1.0/delta_t)*A_C + A_K !A_C + delta_t*A_K
        rhs_time = A_F 
       ! print*, 'AK_time Matrix before BV'
       ! write(*,"(I2,1x, 16F8.5)") (i, (AK_time(i,j), j=1,ntotv), i=1,ldAKban)
       ! print*," "
       ! print*, 'rhs_time before BV'
       ! do i =1, ntotv
       !   print '(I2, 1x, E13.5)', i, rhs_time(i,1)
       ! end do
        
        call ApplyBVs(nofix,ifpre,presc,AK_time,rhs_time)
        
        !------------- Solver -------------!
        u_fut = rhs_time   !here mkl will rewrite u_fut by the solution vector
        AK_LU = AK_time  
       ! print*, 'AK_LU Matrix after BV'
       ! write(*,"(I2,1x, 16F8.5)") (i, (AK_time(i,j), j=1,ntotv), i=1,ldAKban)
       ! print*," "
       ! print*, 'rhs_time after BV'
       ! do i =1, ntotv
       !   print '(I2, 1x, E13.5)', i, rhs_time(i,1)
       ! end do
        
        !--- Factorizing Matrix
        call dgbtrf(S_m, S_n, lowban, upban, AK_LU, ldAKban, S_ipiv, info)
        if(info.ne.0)then
          call MKLfactoResult('dgbtrf',info) 
          print'(A32,I3)', '<<<Error in factorization at time: ', time
        endif
        !--- Solving System of Eqns
        call dgbtrs(S_trans, S_n, lowban, upban, S_nrhs, AK_LU, ldAKban, S_ipiv, u_fut, S_ldSol, info )
        if(info.ne.0)then
          call MKLsolverResult('dgbtrs',info) 
          print'(A32,I3)', '<<<Error in solving system of equation at time: ', time
        endif
        u_pre = u_fut
        !print*, ' '
        !print*,"u_fut"
        !do i =1, ntotv
        !  print '(I2, 1x, E13.5)', i, u_fut(i,1)
        !end do
        !print*, ' '
        !print*, 'Avanzar? 1=yes ; other stop'
        !read(*,*) ans
        !if(ans.eq.1)then
        !  continue
        !else
        !  stop
        !endif
        
        !---------- Printing and writing results -----------!
        !call file_name_inc(solution_file_name)
        print'(A11,I3,A3,F8.3,A5,F8.3,A5)',' time step:',time,' =',nt,'   of',time_fin,' seg'
        call GID_PostProcess(u_pre, File_PostRes, 'res', time, nt, time_fin)
        
        
      end do
      
      print*, ' '
      print*, 'Shape of Global K: ',shape(AK_time)
      print*, 'Shape of Global F: ',shape(rhs_time)
      print*, 'Shape of Solution: ',shape(u_pre)
      write(*,*)
      DEALLOCATE( AK_time, rhs_time, u_pre, u_fut)
      
    end subroutine BackwardEuler
    
    
    
  !end contains

end module timeInt

!Estas lineas son para revisar la solucion e implementar el refinamiento de la solucion,
!Deben agregarse despues de u_pre = u_fut
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
        
