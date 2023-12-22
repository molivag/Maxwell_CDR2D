module timeInt
  use param
  use library
  use E0field
  use sourceTerm

  contains

    subroutine initialCondition(presc,ifpre, nofix, shapeTime, Uinit) 
      
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
     
      !external                                      :: dgbtrf, dgbtrs, dgbrfs
     
      double precision, dimension(ndofn,nBVs)       ,intent(in) :: presc
      integer, dimension(ndofn,nBVs)                ,intent(in) :: ifpre
      integer, dimension(nBVs)                      ,intent(in) :: nofix
      double precision, allocatable, dimension(:,:)             :: dummy
      double precision, dimension(t_steps)                      :: u
      integer                                                   :: tw, t, itotv
      double precision                                          :: tEz
      double precision, dimension(ntotv)                        :: E0
      double precision, dimension(ntotv,1)         ,intent(out) :: Uinit
      double precision, dimension(t_steps)         ,intent(out) :: shapeTime
      !double precision                             ,intent(out):: delta_t
      
      allocate( dummy(ldAKban,ntotv) )
      
      u  = 1.0
      tw = 1 !time*width how strong the impulse is
      Uinit = 0
      tEz = time_ini
      t = 0
      
      ! call Efield_WholeSpace(t,tEz, E0)
      ! do itotv=1,ntotv
      !   Uinit(itotv,1) = E0(itotv)
      ! enddo
      ! ! call ApplyBVs(nofix,ifpre,presc,dummy,Uinit)
      
      !--Select a signal shape function in time
      select case(signal)
        case(1) !step-on
          do t = 1, t_steps
            if(t.le.tw)then
              u(t) = 0.0
            !elseif(t.gt.10*tw)then!.and.t.le.2*tw)then
            !  u(t) = 0.0
            !elseif(t.gt.2*tw)then
            ! u(t) = 1.0
            endif
          end do
        case(2) !step off
          do t = 1, t_steps
            if(t.le.tw)then
              u(t) = 1.0
            elseif(t.gt.tw.and.t.le.2*tw)then
              u(t) = 0.0
            elseif(t.gt.2*tw)then
              u(t) = 0.0
            endif
          end do
        case(3) !triangular
          do t = 1, t_steps
            if(t.le.tw)then
              u(t) = 1.0/tw * t
            elseif(tw.lt.t.and.t.lt.2*tw)then
              u(t) = 1. - 1./tw * (t-tw)
            !else
            !  u(t) = 0.0
            endif
          end do
        case default
          write(*,*) 'No function in time defined'
          stop
      end select
      
      shapeTime = u 
      
      
      !do t =1,t_steps
      !  print"(I0, 1x, f5.3)", t, shapeTime(t)
      !end do
      
      
      
      deallocate(dummy)
      return
      stop 
    end subroutine initialCondition
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine Timeintegration(k_y, basfun,dN_dxi,dN_deta, hes_xixi, hes_xieta, hes_etaeta, &
      &                        nofix, ifpre, presc, Ex_field)
      
      implicit none
      
      external                                            :: dgbtrf, dgbtrs, dgbrfs
      
      double precision, dimension(nne,TotGp) ,intent(in)  :: basfun, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp) ,intent(in)  :: hes_xixi, hes_xieta, hes_etaeta
      double precision, dimension(ndofn,nBVs),intent(in)  :: presc
      integer         , dimension(ndofn,nBVs),intent(in)  :: ifpre
      integer         , dimension(nBVs)      ,intent(in)  :: nofix

      ! - - Local Variables - -!
      character(len=1)                                    :: S_trans
      double precision, allocatable, dimension(:,:)       :: A_K, A_C, A_F
      double precision, allocatable, dimension(:,:)       :: LHS, lhs_BDF2
      double precision, allocatable, dimension(:,:)       :: u_pre, u_curr, u_fut, Mu_pre
      double precision, allocatable, dimension(:)         :: S_ferr, S_berr, S_work
      double precision             , dimension(ntotv,1)   :: Jsource, Jsource_pre
      double precision             , dimension(ntotv,1)   :: RHS, F_plus_MU, rhs_BDF2, u_init
      double precision             , dimension(t_steps+1) :: shapeTime
      integer         , allocatable, dimension(:)         :: S_ipiv, S_iwork
      double precision                                    :: nt, ttt, k_y
      integer                                             :: time, info, workdim
      integer                                             :: ii, S_m, S_n, S_nrhs, S_ldSol!, t_steps
      double precision, allocatable,dimension(:),intent(out) :: Ex_field
      
      
      allocate( LHS(ldAKban,ntotv), lhs_BDF2(ldAKban,ntotv))
      print'(A28,E11.5)', ' -The current wave number is: ', k_y
      
      call GlobalSystem(k_y, basfun, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, A_C, A_K, A_F)
      !allocate( RHS(ntotv,1), Jsource(ntotv,1), F_plus_MU(ntotv,1), rhs_BDF2(ntotv,1), u_init(ntotv,1) )
      !----- Setting MKL-Solver Parammeters -----!
      S_m     = size(A_K,2)  !antes ntotv
      S_n     = size(A_K,2)  !antes ntotv
      S_ldSol = max(1,S_n)
      S_trans = 'N'
      S_nrhs  = 1
      
      allocate( u_pre(S_ldSol, 1), u_fut(S_ldSol, 1) )
      allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))
      allocate( Ex_field(t_steps+1))
      
      u_curr = 0.0; u_pre  = 0.0; u_fut  = 0.0
      Mu_pre = 0.0; RHS    = 0.0; LHS    = 0.0
      time   = 0  ; nt = time_ini;  ttt    = 0.0
      !delta_t 1e-3!( time_fin - time_ini ) / (t_steps + 1.0)   !Step size
      
      call initialCondition(presc,ifpre, nofix, shapeTime, u_init)
      
      u_pre  = u_init                                   !u in present time 
      !write(*,*) ' '
      call GID_PostProcess(1,u_pre, 'msh'    , time, nt, time_fin, Ex_field)
      call GID_PostProcess(1,u_pre, 'res'    , time, nt, time_fin, Ex_field)
      !call GID_PostProcess(1,u_pre, 'profile', time, nt, time_fin, Ex_field)
      write(*,*) ' '
      print*, 'Starting time integration. . . . .'
      write(*,*) ' '
      ! print 100,' time step:',time,'  = ',time_ini,' is the value of u by the initial condiction'
     
      select case(theta)
        !-------- 1st-order Backward Difference 
        case(2) 
          !Time-stepping
          write(*,*)'        < < < BDF1 Selected > > >'
          !do while(ttt < time_fin)
          do time = 1, t_steps
            nt = nt + delta_t!,time_fin,delta_t
            !time = time+1
            
            call prevTime(k_y,basfun,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,S_ldsol,u_pre,Mu_pre)
            LHS  = (A_C + delta_t*A_K)
            call currDensity(2,time,shapeTime(time),Jsource) 
            RHS = (delta_t*A_F + Mu_pre - delta_t*Jsource)
            ! RHS = ( Mu_pre - delta_t*Jsource )
            call ApplyBVs(nofix,ifpre,presc,LHS,RHS)
            
            !------------- Solver -------------!
            u_fut = RHS   !here mkl will rewrite u_fut by the solution vector
            !--- Factorizing Matrix
            call dgbtrf(S_m, S_n, lowban, upban, LHS, ldAKban, S_ipiv, info)
            call checkMKL('f',time,info)
            !--- Solving System of Eqns
            call dgbtrs(S_trans,S_n,lowban,upban,S_nrhs,LHS,ldAKban,S_ipiv,u_fut,S_ldSol,info )
            call checkMKL('s',time,info)
           
            !---------- Printing and writing results -----------!
            call infoTime(time)
            call GID_PostProcess(1,u_fut, 'res'    , time, nt, time_fin, Ex_field)
            ! call GID_PostProcess(1, u_fut, 'profile', time, nt, time_fin, Ex_field)
            ! call GID_PostProcess(1, u_fut, 'spatial', time, nt, time_fin, Ex_field)
            
            !---------- Updating Variables ---------------------! 
            !Jsource_pre = Jsource
            !ttt = ttt+delta_t
            u_pre = u_fut
            
          end do
        !-------- Crank- Nicholson Scheme 
        case(3)
          write(*,*)'    Crank-Nicholson method Selected'  
          !!Time-stepping
          !do time = 1, t_steps +1
          !  nt = nt + delta_t!,time_fin,delta_t  
          !  call GlobalSystem_Time(basfun,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,S_ldsol,delta_t,u_pre,A_F)
          !  !call GlobalSystem_Time(N, dN_dxi, dN_deta, Hesxieta, S_ldsol, delta_t, u_pre, A_F)
          !  LHS  = (1.0/delta_t)*A_C + 0.5*A_K !A_C + delta_t*A_K
          !  RHS = A_F
          !  !duda en 0.5*(A_Fnext + A_Fcurr) dependencia del tiempo del rhs
          !  call ApplyBVs(nofix,ifpre,presc,LHS,RHS)
          !  
          !  !------------- Solver -------------!
          !  u_fut = RHS   !here mkl will rewrite u_fut by the solution vector
          !  !--- Factorizing Matrix
          !  call dgbtrf(S_m, S_n, lowban, upban, LHS, ldAKban, S_ipiv, info)
          !  if(info.ne.0)then
          !    call MKLfactoResult('dgbtrf',info) 
          !    print'(A32,I3)', '<<<Error in factorization at time: ', time
          !  endif
          !  !--- Solving System of Eqns
          !  call dgbtrs(S_trans, S_n, lowban, upban, S_nrhs, LHS, ldAKban, S_ipiv, u_fut, S_ldSol, info )
          !  if(info.ne.0)then
          !    call MKLsolverResult('dgbtrs',info) 
          !    print'(A32,I3)', '<<<Error in solving system of equation at time: ', time
          !  endif
          !  !---------- Printing and writing results -----------!
          !  print'(A11,I3,A3,F8.3,A5,F8.3,A5)',' time step:',time,' =',nt,'   of',time_fin,' seg'
          !  !if(time.eq.1)print'(A11,I3,A3,F8.3,A5,F8.3,A5)',' time step:',time,' =',nt,'   of',time_fin,' seg'
          !  !if(time.eq.t_steps+1)print'(A11,I3,A3,F8.3,A5,F8.3,A5)',' time step:',time,' =',nt,'   of',time_fin,' seg'
          !  call GID_PostProcess(u_fut,'res', time, nt, time_fin)
          !  u_pre = u_fut
          !  
          !end do      
        case(4)
          write(*,*)'BDF2 Selected'
          
          !!################## Desde aqui 
          !!---1st-Order Backward Euler   
          !LHS  = (1.0/delta_t)*A_C + A_K !A_C + delta_t*A_K
          !call GlobalSystem_Time(basfun,dN_dxi, dN_deta,hes_xixi,hes_xieta,hes_etaeta,S_ldsol, delta_t, u_pre, A_F)
          !!call GlobalSystem_Time(N, dN_dxi, dN_deta, Hesxieta, S_ldsol, delta_t, u_pre, A_F)
          !RHS = A_F 

          !call ApplyBVs(nofix,ifpre,presc,LHS,RHS)

          !!#Aqui entraria LHS, RHS para dar el tiempo current, el prev es el de afuera de loop de tiempo y antes
          !! de BDF1 y luego ya solo se actualiza: ucurr=u_fut  y uprebv=ucurr de la iteracion anterior, la calculada con BDF1
          !
          !!------------- Solver -------------!
          !u_curr = RHS   !here mkl will rewrite u_fut by the solution vector
          !!--- Factorizing Matrix
          !call dgbtrf(S_m, S_n, lowban, upban, LHS, ldAKban, S_ipiv, info)
          !if(info.ne.0)then
          !  call MKLfactoResult('dgbtrf',info) 
          !  print'(A32,I3)', '<<<Error in factorization at time: ', time
          !endif
          !!--- Solving System of Eqns
          !call dgbtrs(S_trans, S_n, lowban, upban, S_nrhs, LHS, ldAKban, S_ipiv, u_curr, S_ldSol, info )
          !if(info.ne.0)then
          !  call MKLsolverResult('dgbtrs',info) 
          !  print'(A32,I3)', '<<<Error in solving system of equation at time: ', time
          !endif
          !
          !!---------- Printing and writing results -----------!
          !print'(A11,I3,A3,F8.3,A5,F8.3,A5)',' time step:',time,' =',nt,'   of',time_fin,' seg'
          !call GID_PostProcess(u_curr, 'res', 1, 0.0, time_fin)
          !!call GID_PostProcess(u_curr, File_PostRes, 'res', 1, 0.0, time_fin)
          !!#################  y hast aqui debe ir afuera del loop de tiempo de BDF2 

          !!Una vez calculado el tiempo presente Ucurr a partir del tiempo pasado Uprev
          !!usando BDF1. Calculamos BDF2
          !
          !do time = 2, t_steps +1 !Empieza en 2?
          !  nt = nt + delta_t!,time_fin,delta_t

          !  ! Uprev(Un-1) --> computed throug the initial cond for time 1   
          !  ! Ucurr(Un) --> computed throug 1st order Euler
          !  !Tghe next subrotuine perform: M*(Un - Un-1)
          !  ! - - - 2nd-Order Backward Euler 
          !  lhs_BDF2  = (3/2*delta_t)*A_C + A_K      !#Este tambien podria ir fuera del loop de tiempo
          !  call TimeLevels(basfun, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, delta_t, u_curr, u_pre, F_plus_MU)
          !  !By the last subroutine we compute the product of right-hand side of 2nd order backward Euler 
          !  !done in elemental form to match the dimension btwn band matrix and global vecotr
          !  rhs_BDF2 = F_plus_MU 

          !  call ApplyBVs(nofix,ifpre,presc, lhs_BDF2, rhs_BDF2)

          !  !------------- Solver for future time -------------!
          !  u_fut = rhs_BDF2   !here mkl will rewrite u_fut by the solution vector  
          !  !--- Factorizing Matrix
          !  call dgbtrf(S_m, S_n, lowban, upban, lhs_BDF2, ldAKban, S_ipiv, info)
          !  if(info.ne.0)then
          !    call MKLfactoResult('dgbtrf',info) 
          !    print'(A32,I3)', '<<<Error in factorization at time: ', time
          !  endif
          !  !--- Solving System of Eqns
          !  call dgbtrs(S_trans, S_n, lowban, upban, S_nrhs, lhs_BDF2, ldAKban, S_ipiv, u_fut, S_ldSol, info )
          !  if(info.ne.0)then
          !    call MKLsolverResult('dgbtrs',info) 
          !    print'(A32,I3)', '<<<Error in solving system of equation at time: ', time
          !  endif
          !  !---------- Printing and writing results -----------!
          !  !print'(A11,I3,A3,F8.3,A5,F8.3,A5)',' time step:',time,' =',nt,'   of',time_fin,' seg'
          !  call GID_PostProcess(u_fut, 'res', time, nt, time_fin)
          !  u_pre=u_curr
          !  u_curr=u_fut
          !
          !end do
          
        case default
           write(*,*) 'Not time integration method definded'
           stop
      end select
      
      
      !print*, 'Shape of Global K: ',shape(LHS)
      !print*, 'Shape of Global F: ',shape(RHS)
      !print*, 'Shape of Solution: ',shape(u_pre)
      !write(*,*)
      DEALLOCATE( LHS, u_pre, u_fut)
      
      
      100 format(A11,I4,1x,A3,e12.5,A) 
      101 format(A11,I4,1x,A3,1x,E12.5,A5,1x, E12.5,A5)
    end subroutine TimeIntegration
    
    
    
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
        !&            RHS, S_ldSol, u_fut, S_ldSol, S_ferr, S_berr, S_work, S_iwork, info )
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
        
