module timeInt
  use param
  use library!, only: ApplyBVs, GlobalSystem_Time, file_name_inc, GID_PostProcess, MKLsolverResult, MKLfactoResult
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
     
      double precision, dimension(ndofn,nBVs)      ,intent(in) :: presc
      integer, dimension(ndofn,nBVs)               ,intent(in) :: ifpre
      integer, dimension(nBVs)                     ,intent(in) :: nofix
      double precision, allocatable, dimension(:,:)            :: dummy
      double precision, dimension(t_steps+1)                  :: u
      integer                                                  :: tw, t
      double precision, dimension(ntotv,1)         ,intent(out):: Uinit
      double precision, dimension(t_steps+1)      ,intent(out):: shapeTime
      !double precision                             ,intent(out):: delta_t
      
      allocate( dummy(ldAKban,ntotv) )
      
      u  = 0.0
      tw = 1 !time*width how strong the impulse is
      !delta_t 1e-3!( time_fin - time_ini ) / (t_steps + 1.0)   !Step size
      
      call ApplyBVs(nofix,ifpre,presc,dummy,Uinit)
      
      
      !--Select a signal shape function in time
      select case(signal)
        case(1) !step-on
          do t = 1, t_steps+1
            if(t.le.tw)then
              u(t) = 0.0
            elseif(t.gt.tw)then!.and.t.le.2*tw)then
              u(t) = 1.0
            !elseif(t.gt.2*tw)then
             ! u(t) = 1.0
            endif
          end do
        case(2) !step off
          do t = 1, t_steps+1
            if(t.le.tw)then
              u(t) = 1.0
            elseif(t.gt.tw.and.t.le.2*tw)then
              u(t) = 0.0
            elseif(t.gt.2*tw)then
              u(t) = 0.0
            endif
          end do
        case(3) !triangular
          do t = 1, t_steps+1
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
      
      !do t =1,t_steps+1
      !  print"(I0, 1x, f5.3)", t, shapeTime(t)
      !end do
      
      
      
      deallocate(dummy)
      return
      stop 
    end subroutine initialCondition
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine Timeintegration(basfun,dN_dxi,dN_deta, hes_xixi,hes_xieta,hes_etaeta, &
      &                        time_ini,time_fin,t_steps, nofix, ifpre, presc,&
      &                        S_m, S_n, S_trans, S_nrhs, S_ipiv, S_ldSol, workdim, Ex_field)
      
      implicit none
      
      external :: dgbtrf, dgbtrs, dgbrfs
      
      double precision, dimension(nne,TotGp) ,intent(in)     :: basfun, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp) ,intent(in)     :: hes_xixi, hes_xieta, hes_etaeta
      double precision, dimension(ndofn,nBVs),intent(in)     :: presc
      integer, dimension(ndofn,nBVs)         ,intent(in)     :: ifpre
      integer, dimension(nBVs)               ,intent(in)     :: nofix
      character(len=1)                       ,intent(in)     :: S_trans
      integer                                ,intent(in)     :: S_m, S_n, S_nrhs, S_ldSol, t_steps
      double precision                       ,intent(in)     :: time_ini, time_fin
      ! - - Local Variables - -!
      double precision, allocatable, dimension(:,:)          :: A_K, A_C, A_F
      double precision, allocatable, dimension(:,:)          :: AK_time, lhs_BDF2


      double precision, dimension(ntotv,1)                   :: Jsource, Jsource_pre
      double precision, dimension(ntotv,1)                   :: rhs_time, F_plus_MU, rhs_BDF2, u_init
      double precision, dimension(t_steps+1)                :: shapeTime
      
      double precision, allocatable, dimension(:,:)          :: u_pre, u_curr, u_fut
      double precision, allocatable, dimension(:)            :: S_ferr, S_berr, S_work
      integer         , allocatable, dimension(:)            :: S_ipiv, S_iwork
      double precision                                       :: nt, ttt
      integer                                                :: time, info, workdim
      integer :: ii
      double precision, allocatable,dimension(:),intent(out) :: Ex_field
      
      
      allocate( AK_time(ldAKban,ntotv), lhs_BDF2(ldAKban,ntotv))
      
      !allocate( rhs_time(ntotv,1), Jsource(ntotv,1), F_plus_MU(ntotv,1), rhs_BDF2(ntotv,1), u_init(ntotv,1) )
      
      allocate( u_pre(S_ldSol, 1))
      allocate( u_fut(S_ldSol, 1))
      allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))
      allocate( S_work(workdim), S_iwork(S_ldSol), S_ferr(S_nrhs), S_berr(S_nrhs) )
      allocate( Ex_field(t_steps+1))
      
      u_curr   = 0.0
      u_pre    = 0.0
      u_fut    = 0.0
      rhs_time = 0.0
      AK_time  = 0.0
      time     = 0                                            !initializing the time
      !delta_t 1e-3!( time_fin - time_ini ) / (t_steps + 1.0)   !Step size
      nt       = time_ini
      ttt      = 0.0
      
      call initialCondition(presc,ifpre, nofix, shapeTime, u_init)
      !call initialCondition(presc,ifpre,nofix,basfun,dN_dxi,dN_deta,S_m,S_n,S_trans,S_nrhs,S_ipiv,S_ldSol,delta_t,Uinit) 
      
      u_pre  = u_init                                   !u in present time 
      Jsource_pre = 0.0 !or equal to initial condition which actually is the value of u at boundary conds?
      
      call GlobalSystem(basfun, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, A_C, A_K, A_F)
      
      !write(*,*) ' '
      print 100,' time step:',time,'  = ',time_ini,' is the value of u by the initial condiction'
      call GID_PostProcess(1,u_pre, 'msh'    , time, nt, time_fin, Ex_field)
      call GID_PostProcess(1,u_pre, 'res'    , time, nt, time_fin, Ex_field)
      call GID_PostProcess(1,u_pre, 'profile', time, nt, time_fin, Ex_field)
      print*, 'Starting time integration. . . . .'
      write(*,*) ' '
     
      select case(theta)
        !-------- 1st-order Backward Difference 
        case(2) 
          !Time-stepping
          write(*,*)'              BDF1 Selected'
          !do while(ttt < time_fin)
           
          do time = 1, t_steps
            u_pre = u_fut
            nt = nt + delta_t!,time_fin,delta_t
            !time = time+1
            
            call GlobalSystem_Time(basfun,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,S_ldsol,u_pre,A_F)
            AK_time  = (A_C/delta_t) + A_K 
            
            call ApplyBVs(nofix,ifpre,presc,AK_time,A_F)
            call currDensity(time,shapeTime(4),Jsource) 
            
            rhs_time =  A_F + Jsource/delta_t
            !rhs_time =  A_F + 1/delta_t*(Jsource + Jsource_pre)
            !print'(f15.5)', 1/delta_t*(Jsource + Jsource_pre) 
            
            !------------- Solver -------------!
            u_fut = rhs_time   !here mkl will rewrite u_fut by the solution vector
            !--- Factorizing Matrix
            call dgbtrf(S_m, S_n, lowban, upban, AK_time, ldAKban, S_ipiv, info)
            if(info.ne.0)then
              call MKLfactoResult('dgbtrf',info) 
              print'(A32,I3)', '<<<Error in factorization at time: ', time
            endif
            !--- Solving System of Eqns
            call dgbtrs(S_trans, S_n, lowban, upban, S_nrhs, AK_time, ldAKban, S_ipiv, u_fut, S_ldSol, info )
            if(info.ne.0)then
              call MKLsolverResult('dgbtrs',info) 
              print'(A32,I3)', '<<<Error in solving system of equation at time: ', time
            endif
           
            !---------- Printing and writing results -----------!
            print 101,' time step:',time,' =',nt,'   of',time_fin,' seg'
            call GID_PostProcess(1,u_fut, 'res'    , time, nt, time_fin, Ex_field)
            call GID_PostProcess(1, u_fut, 'profile', time, nt, time_fin, Ex_field)
            
            !---------- Updating Variables ---------------------! 
            Jsource_pre = Jsource
            
            !ttt = ttt+delta_t
            
          end do
        !-------- Crank- Nicholson Scheme 
        case(3)
          write(*,*)'    Crank-Nicholson method Selected'  
          !!Time-stepping
          !do time = 1, t_steps +1
          !  nt = nt + delta_t!,time_fin,delta_t  
          !  call GlobalSystem_Time(basfun,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,S_ldsol,delta_t,u_pre,A_F)
          !  !call GlobalSystem_Time(N, dN_dxi, dN_deta, Hesxieta, S_ldsol, delta_t, u_pre, A_F)
          !  AK_time  = (1.0/delta_t)*A_C + 0.5*A_K !A_C + delta_t*A_K
          !  rhs_time = A_F
          !  !duda en 0.5*(A_Fnext + A_Fcurr) dependencia del tiempo del rhs
          !  call ApplyBVs(nofix,ifpre,presc,AK_time,rhs_time)
          !  
          !  !------------- Solver -------------!
          !  u_fut = rhs_time   !here mkl will rewrite u_fut by the solution vector
          !  !--- Factorizing Matrix
          !  call dgbtrf(S_m, S_n, lowban, upban, AK_time, ldAKban, S_ipiv, info)
          !  if(info.ne.0)then
          !    call MKLfactoResult('dgbtrf',info) 
          !    print'(A32,I3)', '<<<Error in factorization at time: ', time
          !  endif
          !  !--- Solving System of Eqns
          !  call dgbtrs(S_trans, S_n, lowban, upban, S_nrhs, AK_time, ldAKban, S_ipiv, u_fut, S_ldSol, info )
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
          !AK_time  = (1.0/delta_t)*A_C + A_K !A_C + delta_t*A_K
          !call GlobalSystem_Time(basfun,dN_dxi, dN_deta,hes_xixi,hes_xieta,hes_etaeta,S_ldsol, delta_t, u_pre, A_F)
          !!call GlobalSystem_Time(N, dN_dxi, dN_deta, Hesxieta, S_ldsol, delta_t, u_pre, A_F)
          !rhs_time = A_F 

          !call ApplyBVs(nofix,ifpre,presc,AK_time,rhs_time)

          !!#Aqui entraria AK_time, rhs_time para dar el tiempo current, el prev es el de afuera de loop de tiempo y antes
          !! de BDF1 y luego ya solo se actualiza: ucurr=u_fut  y uprebv=ucurr de la iteracion anterior, la calculada con BDF1
          !
          !!------------- Solver -------------!
          !u_curr = rhs_time   !here mkl will rewrite u_fut by the solution vector
          !!--- Factorizing Matrix
          !call dgbtrf(S_m, S_n, lowban, upban, AK_time, ldAKban, S_ipiv, info)
          !if(info.ne.0)then
          !  call MKLfactoResult('dgbtrf',info) 
          !  print'(A32,I3)', '<<<Error in factorization at time: ', time
          !endif
          !!--- Solving System of Eqns
          !call dgbtrs(S_trans, S_n, lowban, upban, S_nrhs, AK_time, ldAKban, S_ipiv, u_curr, S_ldSol, info )
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
      
      
      !print*, 'Shape of Global K: ',shape(AK_time)
      !print*, 'Shape of Global F: ',shape(rhs_time)
      !print*, 'Shape of Solution: ',shape(u_pre)
      !write(*,*)
      DEALLOCATE( AK_time, u_pre, u_fut)
      
      
      100 format(A11,I4,1x,A3,F8.5,A)  
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
        
