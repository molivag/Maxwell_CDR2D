module timeInt
  use param
  use geometry
  use library!, only: ApplyBVs, GlobalSystem_Time, file_name_inc, GID_Nodal_Vals, MKLsolverResult, MKLfactoResult

  contains
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    !subroutine initialCondition(presc,ifpre, nofix, N, dN_dxi,dN_deta,S_m, S_n, S_trans, S_nrhs, S_ipiv,S_ldb,delta_t, Uinit) 
    subroutine initialCondition(presc,ifpre, nofix, delta_t, Uinit) 
      
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
     
      double precision, dimension(ndofn,nBVs), intent(in):: presc
      integer, dimension(ndofn,nBVs), intent(in)         :: ifpre
      integer, dimension(nBVs), intent(in)               :: nofix
      double precision, allocatable, dimension(:,:) :: dummy
      double precision, allocatable, dimension(:,:), intent(out) ::  Uinit
      double precision, intent(out)                 :: delta_t
      
      
      allocate( Uinit(ntotv,1), dummy(ldAKban,ntotv))
      
      delta_t  = ( time_fin - time_ini ) / (max_time + 1.0)   !Step size
      
      call ApplyBVs(nofix,ifpre,presc,dummy,Uinit)
      deallocate(dummy)
      return
      
    end subroutine initialCondition
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine Time_RHS(N, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, delta_t, ugl_pre, A_F)
      use sourceTerm
      
      implicit none
      
      double precision, allocatable, dimension(:,:), intent(in out) :: ugl_pre
      double precision, dimension(nne,TotGp), intent(in):: N, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp), intent(in):: hes_xixi, hes_xieta, hes_etaeta
      double precision, dimension(ndofn)        :: EMsource
      double precision, dimension(nne)          :: basis, xi_cor, yi_cor
      double precision, dimension(DimPr,nne)    :: dN_dxy
      double precision, dimension(3,nne)        :: HesXY
      double precision, dimension(DimPr, dimPr) :: Jaco, Jinv
      double precision, dimension(nevab, nevab) :: Ke, Ce, rhs_CN
      double precision, dimension(nevab)        :: Fe, Fe_time, ue_pre, time_cont
      double precision, dimension(3,3)          :: tauma
      double precision, dimension(nne,DimPr)    :: element_nodes
      integer, dimension(nne)                   :: nodeIDmap
      double precision                          :: dvol, hmaxi, detJ, delta_t
      integer                                   :: igaus, ibase, ielem
      double precision, allocatable, dimension(:,:), intent(out)  :: A_F
      
      allocate(A_F(ntotv, 1) )
      
      A_F = 0.0
      do ielem = 1, nelem 
        Ke = 0.0; Fe = 0.0; Ce = 0.0
        call SetElementNodes(ielem, element_nodes, nodeIDmap, xi_cor, yi_cor)
        
        call gather(nodeIDmap, ugl_pre, ue_pre)
        time_cont = ue_pre * 1.0/delta_t
        !do-loop: compute element capacity and stiffness matrix Ke Ce and element vector Fe
        do igaus = 1, TotGp
          call Jacobian( element_nodes, dN_dxi, dN_deta, igaus, Jaco, detJ, Jinv)
          call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, dN_dxy, HesXY)
          dvol = detJ *  weigp(igaus,1)
          hmaxi = elemSize(Jinv)
          do ibase = 1, nne
            basis(ibase) = N(ibase,igaus)
          end do
          
          call source_term(ielem, basis, xi_cor, yi_cor, EMsource)
          call Galerkin(hmaxi, dvol, basis, dN_dxy, EMsource, Ke, Ce, Fe) !amate lo llame Ke
          call TauMat(hmaxi,tauma)
          call Stabilization(dvol, basis, dN_dxy, HesXY, EMsource, tauma, Ke, Fe)
          
          select case(theta)
            case(1) !BDF1
              Fe_time = Fe + matmul(Ce,time_cont)
            case(2) !BDF2
              print*, 'BDF2 not available yet '
              stop
            case(3) !Crank-Nicholson
              rhs_CN  = (1.0/delta_t)*Ce - 0.5*Ke
              Fe_time = 0.5*Fe + matmul(rhs_CN,ue_pre)
            case default
              write(*,'(A)') 'Not theta defined'
          end select
          
          !Fe_time = Fe + matmul(Ce,time_cont)
        end do
        
        call Assemb_Glob_Vec(nodeIDmap, Fe_time, A_F) !Assemble Global Source vector F
        
      end do
      
    end subroutine Time_RHS
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine TimeIntegration(basfun, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta,&
    &                    nofix, ifpre, presc, S_m, S_n, S_trans, S_nrhs, S_ipiv, S_ldSol)
     
      implicit none
      
      external :: dgbtrf, dgbtrs, dgbrfs
      
      double precision, dimension(nne,TotGp), intent(in)  :: basfun, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp), intent(in)  :: hes_xixi, hes_xieta, hes_etaeta
      double precision, dimension(ndofn,nBVs),intent(in)  :: presc
      integer, dimension(ndofn,nBVs), intent(in)          :: ifpre
      integer, dimension(nBVs), intent(in)                :: nofix
      character(len=1), intent(in)                        :: S_trans
      integer,intent (in)                                 :: S_m, S_n, S_nrhs, S_ldSol!, max_time
      !real, intent(in)                                    :: time_ini, time_fin
      ! - - Local Variables - -!
      double precision, allocatable, dimension(:,:) :: A_K, A_C, A_F, Ftime
      double precision, allocatable, dimension(:,:) :: AK_time, AK_LU, Uprev, Unext, rhs_time, Uinit
      !double precision, allocatable, dimension(:)   :: S_ferr, S_berr, S_work
      integer, allocatable, dimension(:)            :: S_ipiv!, S_iwork
      double precision :: delta_t
      integer          :: time, info!, workdim
      real             :: nt 
     
     
      allocate( AK_time(ldAKban,ntotv),  AK_LU(ldAKban,ntotv))
      allocate( rhs_time(ntotv,1), Uinit(ntotv,1) )
      allocate( Uprev(S_ldSol, 1))
      allocate( Unext(S_ldSol, 1))
      allocate( S_ipiv(max(1,min(S_m, S_n)) ))  !size (min(m,n))
     
      !allocate( S_work(workdim), S_iwork(S_ldSol), S_ferr(S_nrhs), S_berr(S_nrhs) )
      
      Uprev = 0.0;  Unext = 0.0;  rhs_time = 0.0;  AK_time  = 0.0;  time = 0 
     
      call initialCondition(presc,ifpre, nofix, delta_t, Uinit) 
      !call initialCondition(presc,ifpre,nofix,basfun,dN_dxi,dN_deta,S_m,S_n,S_trans,S_nrhs,S_ipiv,S_ldSol,delta_t,Uinit) 
      Uprev = Uinit
     
      write(*,*) ' '
      print'(A11,I0,A3,F8.3,A)',' time step: ',time,'  = ',time_ini,' is the value of u by the initial condiction'
      call GID_Time_Results(Uprev, 'res', time, 0.0)
      print*, 'Starting time integration. . . . .'
      write(*,*) ' '
     
      !the allocate of A_K and A_F are inside of GlobalSystem
      !call GlobalSystem(basfun, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, A_C, A_K, A_F)
      !select case(theta)
      !  case(1)! BDF1
      !    AK_time  = (1.0/delta_t)*A_C + A_K 
      !  case(2)! BDF2
      !    print*, 'BDF2 not available yet '
      !    stop
      !  case(3)! Crank-Nicholson
      !    AK_time  = (1.0/delta_t)*A_C + 0.5*A_K !A_C + delta_t*A_K
      !end select  
      
      !-- Time Stepping
      nt = 0.0 
      do time = 1, max_time +1
        nt = nt + delta_t!,time_fin,delta_t
       
        call Time_RHS(basfun, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, delta_t, Uprev, Ftime)
        call GlobalSystem(basfun, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, A_C, A_K, A_F)
        deallocate(A_F)
        select case(theta)
          case(1)! BDF1
            AK_time  = (1.0/delta_t)*A_C + A_K 
          case(2)! BDF2
            print*, 'BDF2 not available yet '
            stop
          case(3)! Crank-Nicholson
            AK_time  = (1.0/delta_t)*A_C + 0.5*A_K !A_C + delta_t*A_K
        end select  
        
        !--- LHS
        rhs_time = Ftime  
        call ApplyBVs(nofix,ifpre,presc,AK_time,rhs_time)
        
        !------------- solver -------------!
        Unext = rhs_time   !here mkl will rewrite Unext by the solution vector
        AK_LU = AK_time  
        
        !--- Factorizing Matrix
        call dgbtrf(S_m, S_n, lowban, upban, AK_time, ldAKban, S_ipiv, info)
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
        !----------- end solver ------------!
        
        !--- Update Time
        Uprev = Unext
        AK_time = 0.0
        !---------- Printing and writing results -----------!
        print'(A12,I0,A,F10.5,A5,F10.5,A5)',' time step: ',time,' =',nt,' of',time_fin,' seg'
        call GID_Time_Results(Uprev, 'res', time, nt)
        
        !if(time.eq.2)stop
        
      end do
      
      !print*, ' '
      !print*, 'Shape of Global K: ',shape(AK_time)
      !print*, 'Shape of Global F: ',shape(rhs_time)
      !print*, 'Shape of Solution: ',shape(Uprev)
      !write(*,*)
      call GID_Time_Results(Uprev,'msh', time, 0.0) 
      print*, ' '
      !print"(1x, A26,A30)", File_Nodal_Vals//'.post.res', 'written succesfully in Pos/ '
      DEALLOCATE( AK_time, rhs_time, Uprev, Unext)
    end subroutine TimeIntegration
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine GID_Time_Results(solution, activity, step_value, interval)
      
      implicit none
      
      character(len=*), parameter    :: fileplace = "Pos/"
      double precision, dimension(ntotv, 1), intent(in) :: solution
      character(*), intent(in)                :: activity
      integer     , intent(in)                :: step_value
      real        , intent(in)                :: interval
      character(len=10)                       :: extension1, extension2
      character(len=15)                       :: Elem_Type
      double precision, dimension(1, ntotv)   :: solution_T
      double precision, dimension(1,nnodes)   :: xcor, ycor
      integer                                 :: inode, ipoin, ielem, ii
      
      solution_T = transpose(solution)
      xcor  = spread(coord(1,:),dim = 1, ncopies= 1)
      ycor  = spread(coord(2,:),dim = 1, ncopies= 1)
      
      extension1 = ".post.msh"
      extension2 = ".post.res"
      
      
      if(ElemType.eq.'QUAD')then
        Elem_Type = 'Quadrilateral'
      elseif(ElemType.eq.'TRIA')then
        Elem_Type = 'Triangle'
      endif
      
      !print*, '!====== Output files ======'
      
      if(activity == "msh")then !quitar este if y acomodar el numero de unidad 
        open(unit=100, file= fileplace//File_Nodal_Vals//extension1, ACTION="write", STATUS="replace")
        
        write(100,902) 'MESH', '"Domain"', 'dimension', DimPr, 'ElemType', Elem_Type, 'Nnode', nne
        write(100,"(A)") '#2D Convection-Diffusion-Reaction'
        write(100,900) '#Element tipe: ', ElemType,'/',ElemType
        write(100,"(A)")'Coordinates'
        write(100,"(A)") '#   No        X           Y'
        do ipoin = 1, nnodes
          write(100,906) ipoin, xcor(1,ipoin), ycor(1,ipoin)
        end do
        write(100,"(A)") 'End Coordinates'
        write(100,"(A)") 'Elements'
        do ielem=1,nelem
          write(100,908) ielem,(lnods(ielem,inode),inode=1,nne)
        end do
        write(100,"(A)") 'End Elements'
        close(100)
        print"(A11,A22,A30)", ' Mesh file ',File_Nodal_Vals//extension1, 'written succesfully in Pos/ '
        
      elseif(activity == "res")then
        if(step_value == 0)then
          open(unit=200, file= fileplace//File_Nodal_Vals//extension2, ACTION="write", STATUS="replace")
          write(200,"(A)") 'GiD Post Results File 1.0'
          write(200,"(A)") '#2D Convection-Diffusion-Reaction'
        else
          continue
        endif
        open(unit=200, file= fileplace//File_Nodal_Vals//extension2, ACTION="write", STATUS="old", position="append")
        ! se escribe el res de las componentes del campo electrico
        select case(ndofn)
          case(1)
            write(200,"(A29, I3, A)") 'Result "DoF" "EM field" ', step_value,' Scalar OnNodes'
            write(200,"(A)") 'ComponentNames "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','             ex '
            !  se escribe el res para el caso escalar de un grado de libertad
            write(200,914)
            do ipoin = 1, nnodes
              write(200,916) ipoin, solution(ipoin, 1)
            end do
            write(200,"(A)") 'End Values'
          case(2)
            write(200,"(A29, I3, A)") 'Result "DoF" "EM field" ', step_value,' Vector OnNodes'
            write(200,"(A)") 'ComponentNames "ex" "ey" "ez" "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','             ex ','               ey '
            do ipoin = 1, nnodes
              write(200,918) ipoin, solution_T(1, ndofn*ipoin-1), solution_T(1,ndofn*ipoin)
            end do
            write(200,"(A)") 'End Values'
          case(3)
            write(200,"(A24, I3, A)") 'Result "DoF" "EM field" ', step_value,' Vector OnNodes'
            write(200,"(A)") 'ComponentNames "ex" "ey" "ez" "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','             ex ','               ey'
            do ipoin = 1, nnodes
              write(200,919) ipoin, solution_T(1, ndofn*ipoin-2), solution_T(1,ndofn*ipoin-1)
            end do
            write(200,"(A)") 'End Values'
            write(200,"(A24, I3, A)") 'Result "P" "Multiplier" ', step_value,' Scalar OnNodes'
            write(200,"(A)") 'ComponentNames "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','     p '
            !  se escribe el res para el caso escalar de un grado de libertad
            write(200,914)
            ii=1
            do ipoin = 3, nnodes*3,3
              write(200,916) ii, solution_T(1,ipoin)
              ii=ii+1
            end do
            write(200,"(A)") 'End Values'
        end select
      else
        write(*,"(A)") ' < < Error > > Postprocess activity must be "msh" or "res" non ', activity
        close(200)
        stop
      end if
      !write(200,*) 'GaussPoints', "GP_1", 'ElemType', Elem_Type 
      !write(200,*) 'Number Of Gauss Points:', totGp
      !write(200,*) 'Natural Coordinates:' "Given"
      !do igaus = 1, totGp
      !  write(200,901) (ngaus(igaus,j), j=2,DimPr)
      !end do
      !write(200,*) 'End GaussPoints'
      !
      !write(555,"(A)") 'Result "E-field" "ANALYSIS" 0 Vector OnGaussPoints "GP_1" '
      !write(200,"(A)") 'ComponentNames "Ex" "Ey" '
      !write(200,"(A)") 'Values'
      !do ielem = 1, nelem
      !  do igaus = 1, totGp
      !    write(200,903) ielem, (grad_sol(ielem,igaus,icomp), icomp=1,2)
      !  end do
      !end do
      !write(200,"(A)") 'End Values'
      
      close(200)
      !la siguiente instruccion debe usarse con nt no con time pero solo es para avanzar
      
      if(interval.eq.time_fin) then
        print*, ' '
        write(*,"(A9,A24,A30)") ' Res file ', File_Nodal_Vals//extension2, 'written succesfully in Pos/ '
      endif
      
      
      900 format(A15, A13, A1, A13)
      901 format(f15.5,1x)
      902 format(A4,1x,A8,1X,A9,1X,I1,1X,A8,1X,A13,A6,1X,I1)
      903 format(I5,2(1x,E15.5))
      906 format(I7,2(3x,f9.4)) !format for msh
      908 format(9(2x,I7) )
      914 format('#',3x,'No',     9x, 'Dof')
      916 format(I7,2x,E15.5)  !format for scalar case
      918 format(I7,3x,E15.5,3x,E15.5) !format for res velocity
      919 format(I7,2(4x,E15.5)) !format for res velocity
      
    end subroutine GID_Time_Results
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    
    
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
        
