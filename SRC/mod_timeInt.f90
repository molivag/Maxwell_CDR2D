module timeInt
  use param
  use tensor_inputs 
  use library
  ! use E0field
  use sourceTerm

  contains

    ! subroutine initialCondition(presc,ifpre, nofix, shapeTime, Uinit) 
    subroutine Initial_Condition(i_WaveNum, presc,ifpre, nofix,  E0) 
              ! Efield_WholeSpace(t,t0, E0) !---> como t0 se da una sola vez, entonces se quita de variable de entrada
              ! Efield_WholeSpace(t, E0)
              
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
      double precision,              dimension(ndofn,nBVs) ,intent(in) :: presc
      integer,                       dimension(ndofn,nBVs) ,intent(in) :: ifpre
      integer,                       dimension(nBVs)       ,intent(in) :: nofix
      integer,                                              intent(in) :: i_WaveNum
      
      ! double precision, allocatable, dimension(:,:)                    :: dummy
      ! integer                                                   :: 
      double precision                                          :: t0
      ! double precision, dimension(ntotv)                        :: E0

      
      ! character(len=*), parameter           :: fileplace3 = "Exact_Sol_TEM/2D_DoubleLine_WholeSpace/"
      double precision, parameter                                      :: pi = 4*atan(1.d0)
      ! double precision, intent(in)          :: tEz             !correspond to nt 
      ! integer         , intent(in)          :: time            !correspond to time 
      double precision                                                 :: sigma, S, x1,y1,x2,y2,x,y
      double precision                                                 :: Curr_x,Curr_y, angle_loop
      double precision                                                 :: facto, theta2, rho1, rho2, exp1, exp2
      double precision                                                 :: Ez_r_x, Ez_r_y
      integer                                                          :: itest, inode, ii, ipoin, itotv
      double precision                                                 :: Mr,Mx,My, x_profile, mu
      double precision, dimension(1,nnodes)                            :: xcor, ycor
      ! double precision, intent(out)                                    :: Ez_r(ntotv)
      double precision, intent(out)                                    :: E0(ntotv,1)
     
      
      
      t0 = time_ini !este es el tiempo inicial dado en input file
      ! t = 0
      
      xcor  = spread(coord(1,:),dim = 1, ncopies= 1)
      ycor  = spread(coord(2,:),dim = 1, ncopies= 1)
      
      select case(initCond)
        case(0)
          if(i_WaveNum ==0 .or. i_WaveNum==1)then
            ! print*, '* * *'
            print*, '- No initial Condition defined, source term gonna be used'
          endif
          continue
          
        case(1) !Double Line:  eq. (17) of paper by Zhang & Liu (2021)
          if(i_WaveNum ==0 .or. i_WaveNum==1)then
            ! print*, '* * *'
            print*, '- Initial condition for the Double Line problem'
          endif
          
          E0 = 0.0 
          mu = 1.0/lambda
          angle_loop = 90.0
          sigma = 0.01
          x1 = coord(1,Srcloc(1))
          y1 = coord(2,Srcloc(1))
          x2 = coord(1,Srcloc(2))
          y2 = coord(2,Srcloc(2))
          S  = abs(x1*x2) ! not needed, cancels out
          
          angle_loop=angle_loop*pi/180. ! theta=-30 deg
          Mr=Icurr(1)*S
          Mx=Mr*sin(angle_loop)
          My=Mr*cos(angle_loop)
          Curr_x= Mx/S
          Curr_y= My/S
          
          ! loop coord.
          facto =4.0*pi*t0
          theta2=mu*sigma/(4.*t0) ! = theta^2
          
          do inode = 1, nnodes  !simple Function
            x = coord(1,inode)
            y = coord(2,inode)
            
            rho1 = sqrt( (x-x1)**2 + (y-y1)**2 )
            rho2 = sqrt( (x-x2)**2 + (y-y2)**2 )
            exp1 = exp(-rho1*rho1*theta2)
            exp2 = exp(-rho2*rho2*theta2)
            ! Source x
            Ez_r_x = -mu*(Curr_x+Curr_y)/facto*exp1
            Ez_r_y = +mu*(Curr_x+Curr_y)/facto*exp2

            !Electric Field
            ! Ez_r(inode) = Ez_r_x + Ez_r_y
            E0(inode,1) = Ez_r_x + Ez_r_y
            ! if(ndofn.eq.1)then
            !   Ez_r(inode) = Ez_r_x + Ez_r_y
            ! else
            !   Ez_r(ndofn*inode-2) = Ez_r_x + Ez_r_y
            ! endif
           
          end do
        
        case(2) !VMD  eq. (125) of paper by Li & Cheng (2023)
          
          
          
          
          
          
          
          
          
          
          
          
          
          case default
            write(*,*) ' No initial condition defined'
            stop
        end select
      
      ! call Efield_WholeSpace(t,t0, E0)
      ! do itotv=1,ntotv
        ! Uinit(itotv,1) = E0(itotv)
      ! enddo
      ! call ApplyBVs(nofix,ifpre,presc,dummy,Uinit)
      
      
      
      
      
      
      
      
      
      
      
      
      !open(unit=5, file= fileplace3//"E0exact.dat", status='old', action='read')
      !write(5,1)x,Ez_r_x,Ez_r_y,Ez_r
      918 format(I7,3x,E14.6,3x,E14.6) !format for res velocity
      906 format(I7,2(3x,f9.4)) !format for msh
      
      
      
      
      
      return
      stop 
    end subroutine Initial_Condition
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine ShapeTimeSignal(i_WaveNum, shapeTime)
      !si esto es una sola variable de salida, puede quedar como una funcion
      
      implicit none
      integer, intent(in)                               ::  i_WaveNum
      integer                                           :: tw, t
      double precision, dimension(t_steps)              :: u
      double precision, dimension(t_steps), intent(out) :: shapeTime
      
      u  = 0.0
      tw = 1 !time*width how strong the impulse is must be grather or equal to t_steps
      ! t  = 0
      
      select case(signal)
        case(1) !step-on (for test is a continuous turned on signal)
          if(i_WaveNum ==0 .or. i_WaveNum==1)then
            print*, '- Step-on source waveform turned on till end of simulation'
            print*, '* * *'
          endif
          
          ! With the step-on signal, a source current is turned on with a step-function
          ! and left on (source current remains J_p = 1). The initial E fields start from E=0 everywhere.
          ! The fields then converge towards the DC value over time.
          
          do t = 1, t_steps
            if(t.ge.tw)then
              u(t) = 1.0
            !elseif(t.gt.10*tw)then!.and.t.le.2*tw)then
            !  u(t) = 0.0
            !elseif(t.gt.2*tw)then
            ! u(t) = 1.0
            endif
          end do
        case(2) !step off
          if(i_WaveNum ==0 .or. i_WaveNum==1)then
            print*, '- Step-off source waveform selected'
            print*, ' '
          endif
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
          if(i_WaveNum ==0 .or. i_WaveNum==1)then
            print*, '* * *'
            print*, '- Triangular source waveform selected'
            print*, ' '
          endif
          do t = 1, t_steps
            if(t.le.tw)then
              u(t) = 1.0/tw * t
            elseif(tw.lt.t.and.t.lt.2*tw)then
              u(t) = 1. - 1./tw * (t-tw)
            else
             u(t) = 0.0
            endif
          end do
        case default
          write(*,*) 'Signal shape in time not defined'
          stop
      end select
      shapeTime = u 
      
    end subroutine ShapeTimeSignal
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine TimeIntegration(i_WaveNum, basfun,dN_dxi,dN_deta, hes_xixi, hes_xieta, hes_etaeta, &
      &                        nofix, ifpre, presc, Ex_field)
      
      implicit none
      
      external                                            :: dgbtrf, dgbtrs, dgbrfs
      
      double precision, dimension(nne,TotGp) ,intent(in)  :: basfun, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp) ,intent(in)  :: hes_xixi, hes_xieta, hes_etaeta
      double precision, dimension(ndofn,nBVs),intent(in)  :: presc
      integer                                ,intent(in)  :: i_WaveNum
      integer         , dimension(ndofn,nBVs),intent(in)  :: ifpre
      integer         , dimension(nBVs)      ,intent(in)  :: nofix

      ! - - Local Variables - -!
      character(len=1)                                    :: S_trans
      double precision, allocatable, dimension(:,:)       :: A_K, A_C, A_F
      double precision, allocatable, dimension(:,:)       :: LHS!, lhs_BDF2
      double precision, allocatable, dimension(:,:)       :: u_pre, u_curr, u_fut, Mu_pre
      ! double precision, allocatable, dimension(:)         :: S_ferr, S_berr, S_work
      double precision, allocatable, dimension(:,:)       :: Jsource, Jsource_pre, dummy
      double precision             , dimension(ntotv,1)   :: RHS, u_init!, F_plus_MU, rhs_BDF2
      ! double precision             , dimension(t_steps+1) :: shapeTime
      double precision             , dimension(t_steps) :: shapeTime
      double precision, allocatable, dimension(:,:)       :: store_Spec
      integer         , allocatable, dimension(:)         :: S_ipiv!, S_iwork
      double precision                                    :: nt, ttt
      integer                                             :: time, info, workdim
      integer                                             :: ii, S_m, S_n, S_nrhs, S_ldSol!, t_steps
      double precision, allocatable,dimension(:),intent(out) :: Ex_field
      
      
      !* * * Aqui se puede hacer un call a una subroutina que reduzca estas lineas * * *!
      if(TwoHalf == 'Y')then !Just if it is dealing with a 2.5D problem, these variable gonna be updating 
        ky_id = nodal_ky(i_WaveNum)
        k_y = WaveNumbers(i_WaveNum) !Value of wave number for current problem (used in reama)
      else
        continue
      endif
      !* * * Aqui se puede hacer un call a una subroutina que reduzca estas lineas * * *!
      
      allocate( LHS(ldAKban,ntotv))
      ! allocate(lhs_BDF2(ldAKban,ntotv))
      
      if(i_WaveNum.gt.1)allocate(mesh_conductivity(nelem))
      call GlobalSystem(basfun, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, A_C, A_K, A_F)
      !allocate( RHS(ntotv,1), Jsource(ntotv,1), F_plus_MU(ntotv,1), rhs_BDF2(ntotv,1), u_init(ntotv,1) )
      !----- Setting MKL-Solver Parammeters -----!
      S_m     = size(A_K,2)  !antes ntotv
      S_n     = size(A_K,2)  !antes ntotv
      S_ldSol = max(1,S_n)
      S_trans = 'N'
      S_nrhs  = 1
      
      allocate( dummy(ldAKban,ntotv) ) ! esto es Global K para poder llamar a ApplyBVs(nofix,ifpre,presc,dummy,Uinit)
      allocate( u_pre(S_ldSol, 1), u_fut(S_ldSol, 1) )
      allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))
      allocate( Ex_field(t_steps+1))
      
      u_curr = 0.0; u_pre  = 0.0  ; u_fut  = 0.0
      Mu_pre = 0.0; RHS    = 0.0  ; LHS    = 0.0
      time   = 0  ; nt = time_ini ; ttt    = 0.0
      !delta_t 1e-3!( time_fin - time_ini ) / (t_steps + 1.0)   !Step size
      
      call Initial_Condition(i_WaveNum, presc,ifpre, nofix, u_init)
      call ShapeTimeSignal(i_WaveNum,shapeTime)
      !Aqui deberia aplicar las condiciones de frontera para el U-inicial?
      ! call ApplyBVs(nofix,ifpre,presc,dummy,Uinit)
      deallocate(dummy)
      u_pre  = u_init                                          !u in present time 
      
      
      
      
      !* * * Aqui se puede hacer un call a una subroutina que reduzca estas lineas * * *!
      !If it is dealing with a 2.5D problem
      if(TwoHalf == 'Y')then
        deallocate(A_F)
        allocate( store_Spec(S_ldSol,t_steps+1) )
        do ii = 1, S_ldSol
          store_Spec(ii,time+1) = u_pre(ii,1) 
        end do
        call storeSpectrum(i_WaveNum, store_Spec)
      endif
      !* * * Aqui se puede hacer un call a una subroutina que reduzca estas lineas * * *!
      
      !write(*,*) ' '
      
      call GID_PostProcess(i_WaveNum, 1,u_pre, 'msh'    , time, nt, time_fin, Ex_field)
      call GID_PostProcess(i_WaveNum, 1,u_pre, 'res'    , time, nt, time_fin, Ex_field)
      call GID_PostProcess(i_WaveNum, 1,u_pre, 'profile', time, nt, time_fin, Ex_field)
      ! call storeSpectrum('TIME',u_fut, time)
      write(*,*) ' '
      if(i_WaveNum == 0)then
        print*, '- Progress. . . . .'
      else
        if(i_WaveNum ==1) print*, 'Starting time integration for all wavenumbers. . . . .'
      endif
      ! print 100,' time step:',time,'  = ',time_ini,' is the value of u by the initial condiction'
     
      select case(theta)
        !-------- 1st-order Backward Difference 
        case(2) 
          !Time-stepping
          if(i_WaveNum == 0)then
            continue
          else
            if(i_WaveNum ==1) then
              write(*,*)'        < < < BDF1 Selected > > >'
              write(*,*) ' - '
            endif
            print'(A34,I0,A3,E11.4)', " !-------> Time Integration for ky",i_WaveNum,"=",k_y
          endif
          !do while(ttt < time_fin)
          time_stepping: do time = 1, t_steps
            nt = nt + delta_t!,time_fin,delta_t
            !time = time+1
            
            call prevTime(basfun,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,S_ldsol,u_pre,Mu_pre)
            LHS  = (A_C + delta_t*A_K)
            call currDensity(i_WaveNum, Jsource,time,shapeTime(time)) 
            if(oper == 'MAXW')then
              if( (TwoHalf=='Y'.or.TwoHalf=='N') .and. initCond.eq.0)then
                !si se trata de un problema 2.5D sin condicion inicial, entonces usa esto:'
                RHS = (Mu_pre - delta_t*Jsource)
              elseif((TwoHalf=='Y'.or.TwoHalf=='N') .and. initCond.ne.0)then
                !si se trata de un problema 2.5D con condicion inicial, entonces Jsource=0 y se usa esto:'
                RHS = (Mu_pre )
              endif
            else
              RHS = (delta_t*A_F + Mu_pre)
            endif
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
            ! print'(I0, 1x, e15.5)',time, nt
            call GID_PostProcess(i_WaveNum,1,u_fut, 'res'    , time, nt, time_fin, Ex_field)
            call GID_PostProcess(i_WaveNum,1, u_fut, 'profile', time, nt, time_fin, Ex_field)
            ! call GID_PostProcess(1, u_fut, 'spatial', time, nt, time_fin, Ex_field)
            
            !---------- Updating Variables ---------------------! 
            !Jsource_pre = Jsource
            !ttt = ttt+delta_t
            u_pre = u_fut
            if(TwoHalf == 'Y')then
              do ii = 1, S_ldSol
                store_Spec(ii,time+1) = u_pre(ii,1) 
              end do
              call storeSpectrum(i_WaveNum, store_Spec)
            end if
            
          end do time_stepping
          DEALLOCATE( A_K, A_C)
          DEALLOCATE( LHS, u_pre, u_fut)
          
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
      if(i_WaveNum == 0)then
        continue
      else
        print'(A,I0)','-----------------------------------------------------------------> End ky',i_WaveNum
        ! print*,' '
      endif
      
      
      !print*, 'Shape of Global K: ',shape(LHS)
      !print*, 'Shape of Global F: ',shape(RHS)
      !print*, 'Shape of Solution: ',shape(u_pre)
      !write(*,*)
      
      
      
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
        
