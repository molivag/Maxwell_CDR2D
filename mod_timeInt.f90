module timeInt
  use param
  use library, only: ApplyBVs, GlobalSystem_Time, file_name_inc, GID_PostProcess, MKLsolverResult, MKLfactoResult

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
      &                      nofix, ifpre, presc, S_m, S_n, S_trans, S_nrhs, S_ipiv, S_ldSol )
      
      implicit none
      
      external :: dgbtrf, dgbtrs
      
      double precision, dimension(nne,TotGp), intent(in)     :: N, dN_dxi, dN_deta
      double precision, dimension(3,nne), intent(in)         :: Hesxieta
      double precision, dimension(ndofn,nBVs), intent(in)    :: presc
      integer, dimension(ndofn,nBVs), intent(in)             :: ifpre
      integer, dimension(nBVs), intent(in)                   :: nofix
      character(len=1), intent(in)                           :: S_trans
      integer,intent (in)                                    :: S_m, S_n, S_nrhs, S_ldSol
      real, intent(in)                                       :: time_ini, time_fin, max_time, u0_cond
      ! - - Local Variables - -!
      double precision, allocatable, dimension(:,:) :: A_K, A_C, A_F
      double precision, allocatable, dimension(:,:) :: AK_time, u_pre, u_fut, rhs_time, u_init
      integer, allocatable, dimension(:)            :: S_ipiv
      character(len=20)                             :: solution_file_name
      double precision :: delta_t
      integer          :: time, info
      real             :: nt 
      
      
      allocate( AK_time(ldAKban,ntotv))
      allocate( rhs_time(ntotv,1), u_init(ntotv,1) )
      allocate(u_pre(S_ldSol, 1))
      allocate( u_fut(S_ldSol, 1))
      allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))
      
      
      delta_t  = ( time_fin - time_ini ) / max_time     !Step size
      u_init = 0.0                                      !Para prueba lo dejo sin subroutina
      time   = 0                                        !initializing the time
      solution_file_name = 'CDR_u0000.post.res'
      write(*,*) ' '
      print*,'!============ TIME DISCRETIZATION ============!'
      write(*,"(A19,4X,F10.3,1X,A10)") ' - Initial time:          ', time_ini,' '
      write(*,"(A19,4X,F10.3,1X,A10)") ' - Final time:            ', time_fin,' '
      write(*,"(A21,4X,F10.3,1X,A10)") ' - Step size(∆t):         ', delta_t,' '
      write(*,"(A19,4X,F10.3,1X,A10)") ' - Number of steps:       ', max_time,' '
      write(*,*) ' '
      call GID_PostProcess(u_pre, File_PostMsh, 'msh', time )
      
      !call initial_condition ( node_num, node_xy, u0_cond, u_init ) !could be 0
      u_pre  = u_init                                   !u in present time 
      
      !time_unit = 101
      !open ( unit = time_unit, file = time_file_name, status = 'replace' )
      !write ( time_unit, '(g14.6)' ) time
      
      
      write(*,*) ' '
      print*, 'Starting time integration. . . . .'
      write(*,*) ' '
      print'(A11,I4,A2,F10.3,A)',' time step:',time,' = ',time_ini,' is the value of the initial condiction'
      call GID_PostProcess(u_pre, Solution_File_Name, 'res', time)
      do nt = time_ini+delta_t,time_fin,delta_t
        
        time = time + 1
        
        call GlobalSystem_Time(N, dN_dxi, dN_deta, Hesxieta, S_ldsol, delta_t, u_pre, A_K, A_C, A_F)
        
        !-------- Implicit Scheme --------!
        AK_time  = 1/delta_t*A_C + A_K
        rhs_time = A_F 
        ! rhs_time = A_F +  (A_C * u_pre * 1/elta_t)
        call ApplyBVs(nofix,ifpre,presc,AK_time,rhs_time)
        u_fut  = A_F                                      !here mkl will rewrite u_fut by the solution vector
        
        !------- Factorizing Matrix -------!
        call dgbtrf(S_m, S_n, lowban, upban, AK_time, ldAKban, S_ipiv, info)
        if(info.ne.0)then
          call MKLfactoResult('dgbtrf',info)   !Aqui agregar el tiempo para en cada tiempo indicar el info de ejecucion
          print'(A32,I3)', '<<<Factorization error in time: ', time
        endif
        !----- Solving System of Eqns -----!
        call dgbtrs( S_trans, S_n, lowban, upban, S_nrhs, AK_time, ldAKban, S_ipiv, u_fut, S_ldSol, info )
        if(info.ne.0)then
          call MKLsolverResult('dgbtrs',info)  !Aqui agregar el tiempo para en cada tiempo indicar el info de ejecucion
          print'(A32,I3)', '<<<Solving System error in time: ', time
        endif
        u_pre = u_fut
        
        !---------- Printing and writing results -----------!
        write(*,*) ' '
        call file_name_inc(solution_file_name)
        print'(A11,I4,A2,F10.3,A3,F10.3,A)',' time step:',time,' = ',nt,' of ',time_fin,' seg'
        call GID_PostProcess(u_pre, solution_file_name, 'res', time)
        
      end do
      
      DEALLOCATE( AK_time, rhs_time, u_pre, u_fut)
    end subroutine BackwardEuler
    
    
    
  !end contains

end module timeInt
