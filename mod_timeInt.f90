module timeInt
  use param
  use library!, only: ApplyBVs, GlobalSystem_Time, file_name_inc, GID_PostProcess, MKLsolverResult, MKLfactoResult

  contains


    subroutine initialCondition( presc, ifpre, nofix, N, dN_dxi, dN_deta, S_m, S_n, S_trans, S_nrhs, S_ipiv, S_ldb, delta_t, Uinit) 

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

      external                                      :: dgbtrf, dgbtrs, dgbrfs

      character(len=1), intent(in)                       :: S_trans
      double precision, dimension(nne,TotGp), intent(in) :: N, dN_dxi, dN_deta
      double precision, dimension(ndofn,nBVs), intent(in):: presc
      integer, dimension(ndofn,nBVs), intent(in)         :: ifpre
      integer, dimension(nBVs), intent(in)               :: nofix
      integer, intent(in)                                :: S_m, S_n, S_nrhs, S_ldb

      double precision, allocatable, dimension(:,:) :: A_C0, A_U0, dummy    
      double precision, dimension(DimPr, dimPr)     :: Jaco
      integer, allocatable, dimension(:)            :: S_ipiv
      double precision, dimension(nevab, nevab)     :: c0e, u0e
      double precision, dimension(nne)              :: basis, ux_cond
      integer, dimension(nne)                       :: nodeIDmap
      real, dimension(nne,DimPr)                    :: element_nodes
      integer                                       :: ielem, igaus, ibase, ii
      double precision                              :: dvol, detJ, aa, bb, cc, x0, xcoor
      double precision, intent(out)                 :: delta_t

      double precision, allocatable, dimension(:,:), intent(out) ::  Uinit
      !SE DEBERIA CAMBIAR LAS VARIABLES DE TIEMPO A VARIABLES  LOCALES PONIENDOLAS COMO SALIDA DUMMY VARIABLE
      ! EN EL CALL DE INPUT DATA PUES OSLO SE USAN EN DOS SUBRUTINAS Y NO TIENE CASO QUE ESTEN EN TODO EL CODIGO

      !allocate( Uinit(ntotv,1),  A_U0(ntotv,1), A_C0(ldAKban,ntotv) ) 
      
      allocate( Uinit(ntotv,1), dummy(ldAKban,ntotv))
      
      delta_t  = ( time_fin - time_ini ) / (max_time + 1.0)   !Step size
      write(*,*) ' '
      print*,'!============ TIME DISCRETIZATION ============!'
      write(*,"(A19,8X,I2,1X,A10)")    ' - Method Selected:       ', theta,' '
      write(*,"(A19,4X,F10.3,1X,A10)") ' - Initial time:          ', time_ini,' '
      write(*,"(A19,4X,F10.3,1X,A10)") ' - Final time:          ', time_fin,' '
      write(*,"(A19,8X,I5,1X,A10)")    ' - Number of steps:       ', max_time,' '
      write(*,"(A21,7X,F10.6,1X,A10)") ' - Step size(∆t):         ', delta_t,' '
      write(*,*) ' '

      !aa = 5.0/7.0
      !bb = (7.0* sqrt(2.0) ) / 300.0    
      !x0 = u0cond      

      ! A_C0 = 0.0; A_U0 = 0.0
      ! do ielem = 1, nelem 
      !   c0e=0.0;  u0e=0.0
      !   call SetElementNodes(ielem, element_nodes, nodeIDmap)
        
      !   do ibase = 1, nne
      !     !xcoor = 
      !     !print*, 'xcoor',ibase, ' ', xcoor
      !     !cc    = (  xcoor - x0 )
      !     ux_cond(ibase) = element_nodes(ibase,1)!xcoor!aa* exp( -(cc/bb)**2 )
      !   end do

      !   !do-loop: compute element capacity and stiffness matrix Ke Ce and element vector Fe
      !   do igaus = 1, TotGp
      !     Jaco = J2D(element_nodes, dN_dxi, dN_deta, igaus)
      !     detJ = m22det(Jaco)
      !     !Jinv = inv2x2(Jaco)
      !     dvol = detJ *  weigp(igaus,1)
      !     !hmaxi = elemSize(Jinv)
      !     do ibase = 1, nne
      !       basis(ibase) = N(ibase,igaus)
      !     end do
          
      !     call Galerkin_Init_Cond(dvol, basis, ux_cond, c0e, u0e) 
      !     !call TauMat(hmaxi,tauma)
      !     !call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, c0e, u0e)
      !     !                                             Se debe estabilizar?
      !   end do
        
      !   call Assemb_Glob_Mat(nodeIDmap, C0e, A_C0)   !Assemble Global Capacity Matrix at t=0 
      !   call Assemb_Glob_Vec(nodeIDmap, u0e, A_U0)   !Assemble Global initial condition vector

      ! end do
  
      ! call dgbtrf(S_m, S_n, lowban, upban, A_C0, ldAKban, S_ipiv, S_info)  
      ! if(S_info.ne.0)then
      !   print*, 'In Matrix Factorization of Initial Condition:'
      !   call MKLfactoResult('dgbtrf',S_info) 
      ! endif
      ! Uinit = A_U0
      ! print*, ' '  
      ! call dgbtrs(S_trans, S_n, lowban, upban, S_nrhs, A_C0, ldAKban, S_ipiv, Uinit, S_ldb, S_info )
      ! if(S_info.ne.0)then
      !   print*, 'In System of Equations of Initial Condition :'
      !   call MKLsolverResult( 'dgbtrs', S_info )
      ! end if
      
      do ii = 1 , nnodes
        Uinit(ii,1)  = coord(ii,2)!Uinit             !u in present time 
      end do

      call ApplyBVs(nofix,ifpre,presc,dummy,Uinit)
      deallocate(dummy)
      return

    end subroutine initialCondition
    
    
    subroutine TimeIntegration(N, dN_dxi, dN_deta, Hesxieta, time_ini, time_fin, max_time,&
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
      real, intent(in)                                       :: time_ini, time_fin
      ! - - Local Variables - -!
      double precision, allocatable, dimension(:,:) :: A_K, A_C, A_F, Ftime
      double precision, allocatable, dimension(:,:) :: AK_time, AK_LU, Uprev, Unext, rhs_time, Uinit
      double precision, allocatable, dimension(:)   :: S_ferr, S_berr, S_work
      integer, allocatable, dimension(:)            :: S_ipiv, S_iwork
      double precision :: delta_t
      integer          :: time, info, workdim
      real             :: nt 


      allocate( AK_time(ldAKban,ntotv),  AK_LU(ldAKban,ntotv))
      allocate( rhs_time(ntotv,1), Uinit(ntotv,1) )
      allocate( Uprev(S_ldSol, 1))
      allocate( Unext(S_ldSol, 1))
      allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))

      allocate( S_work(workdim), S_iwork(S_ldSol), S_ferr(S_nrhs), S_berr(S_nrhs) )

      Uprev = 0.0;  Unext = 0.0;  rhs_time = 0.0;  AK_time  = 0.0;  time = 0 

      call initialCondition(presc, ifpre, nofix, N, dN_dxi, dN_deta, S_m, S_n, S_trans, S_nrhs, S_ipiv, S_ldSol, delta_t, Uinit) 
      Uprev = Uinit
      

      write(*,*) ' '
      print'(A11,I3,A3,F8.3,A)',' time step:',time,'  = ',time_ini,' is the value of u by the initial condiction'
      call GID_PostProcess(Uprev, File_PostRes, 'res', time, 0.0, time_fin)
      print*, 'Starting time integration. . . . .'
      write(*,*) ' '

      !the allocate of A_K and A_F are inside of GlobalSystem
      call GlobalSystem(N, dN_dxi, dN_deta, Hesxieta, A_C, A_K, A_F) 
      select case(theta)
        case(2)!--- LHS BDF1
          AK_time  = (1.0/delta_t)*A_C + A_K 
        case(3)!-- CN
          AK_time  = (1.0/delta_t)*A_C + 0.5*A_K !A_C + delta_t*A_K
      end select  
      
      !-- Time Stepping
      nt = 0.0 
      do time = 1, max_time +1
        nt = nt + delta_t!,time_fin,delta_t

        call TimeContribution(N, dN_dxi, dN_deta, Hesxieta, delta_t, Uprev, Ftime)
        !--- LHS
        rhs_time = Ftime  
        call ApplyBVs(nofix,ifpre,presc,AK_time,rhs_time)
        !------------- Solver -------------!
        Unext = rhs_time   !here mkl will rewrite Unext by the solution vector
        AK_LU = AK_time  
        !--- Factorizing Matrix
        call dgbtrf(S_m, S_n, lowban, upban, AK_LU, ldAKban, S_ipiv, info)
        if(info.ne.0)then
          call MKLfactoResult('dgbtrf',info) 
          print'(A32,I3)', '<<<Error in factorization at time: ', time
        endif
        !--- Solving System of Eqns
        call dgbtrs(S_trans, S_n, lowban, upban, S_nrhs, AK_LU, ldAKban, S_ipiv, Unext, S_ldSol, info )
        if(info.ne.0)then
          call MKLsolverResult('dgbtrs',info) 
          print'(A32,I3)', '<<<Error in solving system of equation at time: ', time
        endif
        !--- Update Time
        Uprev = Unext
        !---------- Printing and writing results -----------!
        print'(A11,I3,A3,F8.3,A5,F8.3,A5)',' time step:',time,' =',nt,'   of',time_fin,' seg'
        call GID_PostProcess(Uprev, File_PostRes, 'res', time, nt, time_fin)

      end do

      print*, ' '
      print*, 'Shape of Global K: ',shape(AK_time)
      print*, 'Shape of Global F: ',shape(rhs_time)
      print*, 'Shape of Solution: ',shape(Uprev)
      write(*,*)
      call GID_PostProcess(Uprev, File_PostMsh, 'msh', time, 0.0, time_fin) 
      DEALLOCATE( AK_time, rhs_time, Uprev, Unext)

    end subroutine TimeIntegration
    
    
  !end contains

end module timeInt

!Estas lineas son para revisar la solucion e implementar el refinamiento de la solucion,
!Deben agregarse despues de Uprev = Unext
        ! print*, ' '
        ! print*,"Solution vector"
        ! do i =1, ntotv
        !   print '(F13.5)', Unext(i,1)
        ! end do
        ! !----- Refining Solution -----!
        ! call dgbrfs(S_trans, S_n, lowban, upban, S_nrhs, A_K, ldAKban, AK_LU, ldAKban, S_ipiv,&
        !&            rhs_time, S_ldSol, Unext, S_ldSol, S_ferr, S_berr, S_work, S_iwork, info )
        ! if(info.ne.0)then
        !   call MKLsolverResult('dgbrfs',info) 
        !   print'(A32,I3)', '<<<The parameter ', info, ' had an illegal value in time: ', time
        ! endif
        ! print*, ' '
        ! print*, 'Refined Solution'
        ! do i =1, ntotv
        !   print '(9F13.5)', Unext(i,1)
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
        
