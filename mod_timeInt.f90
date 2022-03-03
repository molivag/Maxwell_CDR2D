module timeInt
  use param
  use library, only: ApplyBVs, PosProcess, MKLsolverResult, MKLfactoResult

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
    
    
    subroutine BackwardEuler(time_ini, time_fin, max_time, u0_cond, nofix, ifpre, presc,&
      &                      S_m, S_n, S_trans, S_nrhs, S_ipiv, S_ldSol , A_K, A_C, A_F)
      
      implicit none
      
      external :: dgbtrf, dgbtrs
      
      double precision, dimension(ldAKban,ntotv), intent(in) :: A_K, A_C
      double precision, dimension(ndofn,nBVs), intent(in)    :: presc
      double precision, dimension(ntotv,1), intent(in)       :: A_F
      integer, dimension(ndofn,nBVs), intent(in)             :: ifpre
      integer, dimension(nBVs), intent(in)                   :: nofix
      character(len=1), intent(in)                           :: S_trans
      integer,intent (in)                                    :: S_m, S_n, S_nrhs, S_ldSol
      real, intent(in)                                       :: time_ini, time_fin, max_time, u0_cond
      ! - - Local Variables - -!
      double precision, allocatable, dimension(:,:) :: AK_time, u_pre, u_fut, rhs_time, u_init
      integer, allocatable, dimension(:)            :: S_ipiv
      double precision :: delta_t
      integer          :: time, info
      real             :: nt 
      
      
      ! do time=I(3)+k:k:I(4)
      !   t=time;     
      !   fnew=(eval(f ))
      !   fnew=h*fnew
      !   fnew=[bc(1); fnew; bc(2)]
      !   b=theta*fnew+(1-theta)*fold+B*uold
      !   y=L\b
      !   u=U\y
      !   uold=u
      ! end
      
      
      allocate( AK_time(ldAKban,ntotv))
      allocate( rhs_time(ntotv,1), u_init(ntotv,1), u_pre(S_ldSol, 1), u_fut(S_ldSol, 1))
      allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))
      
      ! call initial_condition ( node_num, node_xy, u0_cond, u_init ) !could be 0
      delta_t  = ( time_fin - time_ini ) / max_time     !Step size
      u_init = 0.0                                      !Para prueba lo dejo sin subroutina
      u_pre  = u_init                                   !u in present time 
      u_fut  = 0.0                                      !here mkl will rewrite u_fut by the solution vector
      time   = 0                                        !initializing the time
      
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Initial time = ', time_ini
      write ( *, '(a,g14.6)' ) '  Final time =   ', time_fin
      write ( *, '(a,g14.6)' ) '  Step size =    ', delta_t
      write ( *, '(a,g14.6)' ) '  Number of steps = ', max_time
      
      print'(A6,F10.3,A)','time: ',time_ini,'is equal to the initial condiction'
      print*, 'Starting time-stepping. . . . .'
      do nt = time_ini+delta_t,time_fin,delta_t
        
        time = time + 1
        print'(A11,I4,A2,F10.3,A2,F10.3,A)','time step: ',time,'= ',nt,'of ',time_fin,'seg'
        
        !-------- Implicit Scheme --------!
        AK_time  = 1/delta_t*A_C + A_K
        rhs_time = A_F +  (u_pre * 1/delta_t)
        ! rhs_time = A_F +  (A_C * u_pre * 1/elta_t)
        call ApplyBVs(nofix,ifpre,presc,AK_time,rhs_time)
        
        !------- Factorizing Matrix -------!
        call dgbtrf(S_m, S_n, lowban, upban, AK_time, ldAKban, S_ipiv, info)
        call MKLfactoResult('dgbtrf',info)   !Aqui agregar el paso del tiempo para en cada tiempo indicar el info de ejecucion
        !----- Solving System of Eqns -----!
        call dgbtrs( S_trans, S_n, lowban, upban, S_nrhs, AK_time, ldAKban, S_ipiv, u_fut, S_ldSol, info )
        call MKLsolverResult('dgbtrs',info)  !Aqui agregar el paso del tiempo para en cada tiempo indicar el info de ejecucion
        
        u_pre = u_fut 
        !---------- Print and write results -----------!
        call PosProcess(u_pre, File_PostMsh, 'msh') !se debe agregar el nt como dummyvariable
        call PosProcess(u_pre, File_PostRes, 'res')
        write(*,*)
        
      end do
      
      
      DEALLOCATE( AK_time, rhs_time, u_pre, u_fut)
    end subroutine BackwardEuler
    
    
    
  !end contains

end module timeInt
