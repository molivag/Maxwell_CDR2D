module library
  use param
  use geometry
  use biunit

  contains
    
    subroutine GeneralInfo(name_inputFile)
      external :: fdate
      
      character(len=*), parameter :: fileplace = "Res/results/"
      character(len=19)           :: name_inputFile
      character(len=5)            :: file_name 
      character(len=24)           :: date
      character(len=4)            :: aaaa, cccc
      character(len=12)           :: bbbb
      character(len=16)           :: dddd
      double precision            :: delta_t
      integer :: i,j, k, l
      
      if(kstab.eq.0)then
        aaaa = 'NONE'
      elseif(kstab.eq.1)then
        aaaa = 'SUPG'
      elseif(kstab.eq.2)then
        aaaa = 'GLS'
      elseif(kstab.eq.3 .or. kstab.eq.5)then
        aaaa = 'SGS'
      elseif(kstab.eq.4)then
        aaaa = 'CG'
      elseif(kstab.eq.6)then
        aaaa = 'MAVF'
      else
        write(*,'(A)') '> > >Error in stabilization method'
      endif
     
      if(ProbType.ne.'TIME')ProbType = 'STAT'
     
      call fdate(date)
      print*, ' '
      print*, '- - - - 2D Convection-Diffusion-Reaction Simulation - - - - '
      print*, ' '
      print*,' ',date
      print*,'!================= GENERAL INFO ===============!'
      write(*,"(A19,7x,a19,3X,A1)") ' - Input File:             ', name_inputFile,''
      write(*,"(A19,7x,a5,3X,A1)") ' - Element type:           ', InitElemType,''
      write(*,"(A19,7x,a5,3X,A1)")  ' - Problem Type:           ', ProbType,''
      write(*,"(A19,4X,I6,1X,A10)") ' - Problem dimension:      ', DimPr, '  '
      write(*,"(A19,4X,I6,1X,A10)") ' - Elements:               ', initelem,'   '
      write(*,"(A19,4X,I6,1X,A10)") ' - Nodal points:           ', initnodes, ' '
      write(*,"(A19,4X,I6,1X,A10)") ' - DoF per node:           ', ndofn, '  '
      write(*,"(A19,4X,I6,1X,A10)") ' - Nodes per element:      ', nne, '    '
      write(*,"(A19,4X,I6,1X,A10)") ' - Total Gauss points:     ', totGp,'   '
      write(*,"(A19,4X,I6,1X,A10)") ' - Element variabless:     ', initnevab,'   '
      write(*,"(A19,4X,I6,1X,A10)") ' - Total unknowns:         ', initntotv,'   '
      write(*,"(A19,7X,f8.4,1X,A10)") ' - Element size:        ', 2.0**(-i_exp),'   '
      write(*,"(A26,2X,f3.1,1X,A10)") ' - Length ref. element:    ', hnatu,'   '
      
      if(refiType.eq.'NO')then
        write(*,"(A23,2x,a6,3X,A1)") ' - Refinement type:        ', '  NONE',''
      elseif(refiType.eq.'PS'.or.refitype.eq.'CB')then
        if(refitype.eq.'PS')then
          bbbb = 'Powell-Sabin'
        elseif(refiType.eq.'CB')then
          bbbb = 'Cross-Box'
        end if
          print*, ' '
          print*,'!================= REFINMENT INFO ===============!'
          write(*,"(A23,2x,a12,3X,A1)") ' - Refinement type:    ', bbbb,''
          write(*,"(A23,2X,I6,1X,A10)") ' - Final elements:     ', nelem,'   '
          write(*,"(A23,2X,I6,1X,A10)") ' - Final nodes:        ', nnodes, ' '
      else
        write(*,'(A)') '> > >Error in refinment type'
      endif
      
      print*, ' '
      if(kstab.eq.0)then
        print*,'!========== STABILIZATION PARAMETERS ==========!'
        write(*,"(A26,3x,a4,3X,A1)") ' - Stabilization method:   ', aaaa,''
      elseif(kstab.eq.6)then
        print*,'!========== STABILIZATION PARAMETERS ==========!'
        write(*,"(A26,3x,a4,3X,A1)") ' - Stabilization method:   ', aaaa,''
        write(*,"(A30,3X,f5.2,1X,A10)") ' - Reluctivity of medium (λ): ', lambda, '  '
        write(*,"(A30,2X,f5.2,1X,A10)") ' - Algorithmic constant (Cu): ', Cu, ' '
        write(*,"(A30,4X,f5.2,1X,A10)") ' - Constante of length (ℓ):   ', ell, '    '
        write(*,"(A26,5X,e13.5,1X,A10)") ' - Stab. param.1 (Su):       ', Cu*lambda*(helem**2/ell**2),'   '
        write(*,"(A26,4X,e14.5,2X,A10)") ' - Stab. param.2 (Sp):       ', ell**2 / lambda,'   '
      else
        print*,'!========== STABILIZATION PARAMETERS ==========!'
        write(*,"(A26,3x,a4,3X,A1)") ' - Stabilization method:   ', aaaa,''
        write(*,"(A26,3x,I2,3X,A1)")  ' - Type of Tau matrix:    ', ktaum,''
        write(*,"(A26,3X,f3.1,1X,A10)") ' - Param. to obtain TAU:   ', patau, '  '
      endif
      
      print*, ' '
      if(ProbType.eq.'TIME')then
        delta_t  = ( time_fin - time_ini ) / (max_time + 1.0)   !Step size
        print*,'!============ TIME DISCRETIZATION ============!'
        if((theta.eq.2).or.(theta.eq.4))then
          if(theta.eq.2)then
            cccc = 'BDF1'
          else
            cccc = 'BDF2'
          endif
          write(*,"(A19,8X,A,1X,A10)") ' - Method Selected:        ', cccc,' '
        elseif(theta.eq.3)then
          dddd = 'Cranck-Nicholson'
          write(*,"(A23,2x,a16,3X,A1)") ' - Method Selected:        ', dddd,' '
        endif
        write(*,"(A19,4X,F10.3,1X,A10)") ' - Initial time:          ', time_ini,' '
        write(*,"(A19,4X,F10.3,1X,A10)") ' - Final time:            ', time_fin,' '
        write(*,"(A19,5X,I5,1X,A10)")    ' - Number of steps:       ', max_time,' '
        write(*,"(A21,7X,F10.6,1X,A10)") ' - Step size(∆t):         ', delta_t,' '
      else
        continue
      endif
      
      !print*, ' '
      !print*,'!============ TENSOR COEFFICIENTS  ============!'
      !print*, 'Diffusion'
      !do i = 1,dimPr
      !  do j = 1,DimPr
      !    print"(A,2I1)", 'k_',i,j
      !    do k = 1,ndofn
      !      print"(f15.7,1x,f15.7, 1x, f15.7)",( difma(k,l,i,j), l=1,ndofn)
      !    end do
      !    !print*,' '
      !  end do
      !end do
      !print*, ' '  
      !print*, 'Convection'
      !do k = 1, DimPr
      !  print"(A,2I1)",'A_',k
      !  do i = 1, ndofn
      !    write(*, "(f10.5, 1x, f10.5, 1x, f10.5)")( conma(i,j,k) ,j=1, ndofn)
      !  end do
      !end do
      !  print*,' '
      !print*,'Reaction'
      !do i=1,ndofn
      !  write(*,"(f10.5, 1x, f10.5, 1x, f10.5)" )( reama(i,j) ,j=1,ndofn)
      !end do
      !print*, ' '
      !print*, 'External Forces'
      !if(ndofn.eq.1)then
      !  write(*,"(3(f10.3,1x))") force(1)
      !elseif(ndofn.eq.2)then
      !  write(*,"(2(f10.3,1x))") force(1), force(2)
      !else
      !  write(*,"(3(f10.3,1x))") force(1), force(2), force(3)
      !endif
      
      file_name ="test_"
      open(unit=100,file= fileplace//file_name//testID//'.txt', ACTION="write", STATUS="replace")
      
      if(refiType.eq.'NO')then
        bbbb = '    NONE'
      elseif(refiType.eq.'PS')then
        bbbb = 'Powell-Sabin'
      elseif(refiType.eq.'CB')then
        bbbb = 'Cross-Box'
      else
        write(*,'(A)') '> > >Error in refinment type'
      endif
      
      write(100,'(A)')'- - - - 2D Convection-Diffusion-Reaction Simulation - - - - '
      write(100,'(A)')
      write(100,'(A8,1x,A14)') 'test ID: ',testID
      write(100,'(A)') " "
      write(100,'(A)') date
      write(100,'(A)')'!================= GENERAL INFO ===============!'
      write(100,"(A19,7x,a19,3X,A1)") ' - Input File:             ', name_inputFile,''
      write(100,"(A19,7x,a4,3X,A1)") ' - Element type:           ', InitElemType,''
      write(100,"(A19,5X,I6,1X,A10)") ' - Elements:               ', initelem,'   '
      write(100,"(A19,5X,I6,1X,A10)") ' - Nodal points:           ', initnodes, ' '
      write(100,"(A19,4X,I6,1X,A10)") ' - Nodes per element:      ', nne, '    '
      write(100,"(A19,4X,f8.4,1X,A10)") ' - Mesh size 2^-i:        ', 2.0**(-i_exp),'   '
      write(100,"(A26,3X,f3.1,1X,A10)") ' - Length ref. element:     ', hnatu,'   '
      if(refiType.eq.'NO')then
        write(100,"(A19,4x,a7,3X,A1)") ' - Refinement type:        ', '  NONE',''
        write(100,'(A)')
      elseif(refiType.ne.'NO')then
        write(100,'(A)')'!================= REFINMENT INFO ===============!'
        write(100,"(A22,1x,a12,3X,A1)") ' - Refinement type:         ', bbbb,''
        write(100,"(A14,4x,a13,3X,A1)") ' - Final element type:     ', ElemType,''
        write(100,"(A19,6X,I6,1X,A10)") ' - Total Elements:         ', nelem,'   '
        write(100,"(A23,2X,I6,1X,A10)") ' - Total Nodal points:     ', nnodes, ' '
        write(100,'(A)') 
      endif
      if(kstab.eq.0)then
        write(100,'(A)')'!========== STABILIZATION PARAMETERS ==========!'
        write(100,"(A26,3x,a4,3X,A1)") ' - Stabilization method:       ', aaaa,''  
      !write(100,"(A26,3X,f3.1,1X,A10)") ' - Exponent of mesh size:    ', i_exp,'   '
      elseif(kstab.eq.6)then
        write(100,'(A)')'!========== STABILIZATION PARAMETERS ==========!'
        write(100,"(A29,3x,a4,3X,A1)") ' - Stabilization method:        ', aaaa,''
        write(100,"(A29,3X,f8.3,1X,A10)") ' - Algorithmic constant(Cu): ', Cu, ' '
        write(100,"(A31,3X,f8.3,1X,A10)") ' - Constant of length(ℓ):    ', ell, '    '
        write(100,"(A30,2X,f8.3,1X,A10)") ' - Reluctivity of the medium: ', lambda, '  '
        write(100,"(A29,4X,e13.5,1X,A10)")' - Stab. param.1 (Su):       ', Cu*lambda*(helem**2/ell**2),'   '
        write(100,"(A29,3X,e14.5,2X,A10)")' - Stab. param.2 (Sp):       ', ell**2 / lambda,'   '
      else
        write(100,'(A)')'!========== STABILIZATION PARAMETERS ==========!'
        write(100,"(A26,3x,a4,3X,A1)") ' - Stabilization method:       ', aaaa,''
        write(100,"(A26,3x,I2,3X,A1)")  ' - Type of Tau matrix:        ', ktaum,''
        write(100,"(A26,3X,f3.1,1X,A10)") ' - Param. to obtain TAU:    ', patau, '  '
      endif
    
      if(ProbType.eq.'TIME')then
        delta_t  = ( time_fin - time_ini ) / (max_time + 1.0)   !Step size
        write(100,'(A)') 
        write(100,'(A)')'!============ TIME DISCRETIZATION =============!'
        if((theta.eq.2).or.(theta.eq.4))then
          if(theta.eq.2)then
            cccc = 'BDF1'
          else
            cccc = 'BDF2'
          endif
          write(100,"(A19,8X,A,1X,A10)") ' - Method Selected:        ', cccc,' '
        elseif(theta.eq.3)then
          dddd = 'Cranck-Nicholson'
          write(100,"(A23,2x,a16,3X,A1)") ' - Method Selected:        ', dddd,' '
        endif
        write(100,"(A19,4X,F10.3,1X,A10)") ' - Initial time:          ', time_ini,' '
        write(100,"(A19,4X,F10.3,1X,A10)") ' - Final time:            ', time_fin,' '
        write(100,"(A19,5X,I5,1X,A10)")    ' - Number of steps:       ', max_time,' '
        write(100,"(A21,7X,F10.6,1X,A10)") ' - Step size(∆t):         ', delta_t,' '
      else
        continue
      endif
      write(100,'(A)') 
      write(100,'(A)')'!============ TENSOR COEFFICIENTS  ============!'
      write(100,'(A)') 'Diffusion'
      do i = 1,dimPr
        do j = 1,DimPr
          write(100,"(A,2I1)") 'k_',i,j
          do k = 1,ndofn
            write(100,"(f15.7,1x,f15.7, 1x, f15.7)") ( difma(k,l,i,j), l=1,ndofn)
          end do
          !print*,' '
        end do
      end do
      write(100,'(A)')
      write(100,'(A)') 'Convection'
      do k = 1, DimPr
        write(100,"(A,2I1)")'A_',k
        do i = 1, ndofn
          write(100,"(f10.3, 1x, f10.3, 1x, f10.3)") ( conma(i,j,k) ,j=1, ndofn)
        end do
      end do
      write(100,'(A)') 
      write(100,'(A)') 'Reaction'
      do i=1,ndofn
        write(100,"(f10.3, 1x, f10.3, 1x, f10.3)" ) ( reama(i,j) ,j=1,ndofn)
      end do
        write(100,'(A)') 
      write(100,'(A)') 'External Forces'
      if(ndofn.eq.1)then
        write(100,"(1(f10.3,1x))") force(1)
      elseif(ndofn.eq.2)then
        write(100,"(2(f10.3,1x))") force(1), force(2)
      else
        write(100,"(3(f10.3,1x))") force(1), force(2), force(3)
      endif
      
      
      close(100)
      
    endsubroutine GeneralInfo
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    ! 
    subroutine ReadFile(NumRows, NumCols, condition, condition_value)
      
      integer :: i, j, stat1, stat2
      integer, intent(in)            :: NumRows, NumCols
      character(len=*), parameter    :: fileplace = "./"
      character (len=9)              :: FileName1, FileName2
      character(len=180) :: msg
      integer, dimension(NumRows,NumCols-ndofn), intent (out) :: condition
      double precision, dimension(NumRows,ndofn),intent (out) :: condition_value
      
      FileName1 ='ifpre.dat' 
      FileName2 ='BoVal.dat'
      
      open(unit=10,file=fileplace//FileName1, status='old', action='read', iostat=stat1, iomsg=msg)
      open(unit=20,file=fileplace//FileName2, status='old', action='read', iostat=stat2, iomsg=msg)
      
      read(10,*,iostat=stat1,iomsg=msg) ((condition(i,j), j=1,NumCols-ndofn), i=1,NumRows)
      if (stat1.ne.0) then
        print*, ' '
        print*,'!============ STATUS READING BOUND VAL. ============!'
        print "(A38,I2)", "- Status while reading ifPre file is: ", stat1
        print*, ' '
        print'(A8,1x,A180)','iomsg= ',msg
      else
        continue
      end if
      
      read(20,*,iostat=stat2,iomsg=msg) ((condition_value(i,j), j=1,ndofn), i=1,NumRows)
      if (stat2.ne.0) then
        print "(A38,I2)", "- Status while reading BoVal file is: ", stat2
        print*, ' '
        print'(A8,1x,A180)','iomsg= ',msg
      else
        continue
      end if
      
      
      close (10)
      close (20)
      
    end subroutine ReadFile
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    ! 
    subroutine SetElementNodes(elm_num, element_nodes, nodeIDmap, xi_cor, yi_cor)
      
      implicit none
      
      ! number of element for each elem integral in loop of K global
      integer,intent(in)                                   :: elm_num 
      integer                                              :: i,j,global_node_id
      double precision, dimension(DimPr,nne)               :: aaa
      double precision, dimension(nne,DimPr), intent(out)  :: element_nodes
      double precision,dimension(nne), intent(out)         :: xi_cor, yi_cor
      integer,dimension(nne), intent(out)                  :: nodeIDmap
      
      element_nodes = 0.0
      nodeIDmap = 0
      
      do i = 1, nne
        global_node_id = lnods(elm_num,i)
        nodeIDmap(i)   = global_node_id
        do j=1 ,DimPr
          aaa(j,i) = coord(j,global_node_id)
        end do
        xi_cor(i) = aaa(1,i)
        yi_cor(i) = aaa(2,i)
      end do
      
      element_nodes = transpose(aaa)
      
    end subroutine SetElementNodes
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    ! 
    subroutine Jacobian( element_nodes, dN_dxi, dN_deta, Gp, xjacm, djacb, xjaci)
      implicit none
      
      integer, intent(in)      :: Gp ! se usara en el lazo principal con el punto de Gauss
      double precision, dimension(nne,DimPr), intent(in) :: element_nodes
      double precision, dimension(nne,totGp), intent(in) :: dN_dxi, dN_deta
      double precision, dimension(DimPr,nne)  :: derst!, Basis2D, derst
      double precision, dimension(DimPr,nne)  :: elcod
      double precision, dimension(DimPr,DimPr), intent(out) :: xjacm, xjaci
      double precision                        , intent(out) :: djacb
      double precision, parameter                           :: EPS = 1.0E-10
      integer                                               :: idime, jdime, inode
      
      !double precision, dimension(1,nne)      :: Nxi, Neta
      !double precision, dimension(DimPr,DimPr):: J2D
      
      !con estas instrucciones extraigo la columna de Nx como renglon y lo guardo en Nxi, Gp se
      !ira moviendo conforme la funcion J2D sea llamada en el lazo principal para 
      !cada elemento lo mismo para Neta con dN_deta
      
      !Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)
      !Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)
      
      !Las siguientes tres lineas realizan de forma implicita el calculo de las derivadas
      !espaciales es decir dN/dx and dN/dy (eq. 5.114 - 5.117). Las derivadas espaciales
      !no se calcula explicitamente, en su lugar se usa:
      
      !            d/dy = (d/deta)J^-1      (ver eq. 5.77 y 5.109)
      
      !Basis2D(1,:) = Nxi(1,:)
      !Basis2D(2,:) = Neta(1,:)
      !
      !J2D = matmul(Basis2D,element_nodes)
      
      
      do jdime = 1,2
        do inode = 1,nne
          elcod(jdime,inode) = element_nodes(inode, jdime)
        end do
      end do
     
      do inode = 1,nne
        derst(1,inode) = dN_dxi(inode,Gp) 
        derst(2,inode) = dN_deta(inode,Gp) 
      end do
      !print*," " 
      !print"(A8, 9(1x,f9.5))","dN_dxi: ", (derst(1,inode), inode=1,nne)
      !print"(A9, 9(1x,f9.5))","dN_deta: ", (derst(2,inode), inode=1,nne)
     
      ! Jacobian Matrix
      do idime=1,2
        do jdime=1,2
          xjacm(idime,jdime)=0.0
          do inode=1,nne
            !xjacm(idime,jdime) = xjacm(idime,jdime) + derst(idime,inode)*element_nodes(inode,jdime)
            xjacm(idime,jdime) = xjacm(idime,jdime) + derst(idime,inode)*elcod(jdime,inode)
          end do
        end do
      end do
      
      !do idime=1,2
      !  write(*,"(f10.5, 1x, f10.5)" )( xjacm(idime,jdime) ,jdime=1,2)
      !end do
      
      ! Determinant of Jacobian Matrix
      djacb  = xjacm(1,1)*xjacm(2,2) - xjacm(1,2)*xjacm(2,1)
      !print"(A8,f15.5)", 'det Jaco: ',djacb 
      if(djacb.lt.1.e-8)then
        write(*,*)'Element with non-positive Jacobian'
      end if
      
      ! Inverse of Jacobian
      if (abs(djacb) .le. EPS) then
        xjaci = 0.0! 0.0D0
        return
      end if
      xjaci(1,1)= xjacm(2,2)/djacb
      xjaci(2,2)= xjacm(1,1)/djacb
      xjaci(1,2)=-xjacm(1,2)/djacb
      xjaci(2,1)=-xjacm(2,1)/djacb
      
      return
      
    end subroutine Jacobian
    
    !function djacb(xjacm)
    !  
    !  implicit none
    !  
    !  double precision :: djacb
    !  double precision, dimension(2,2), intent(in)  :: xjacm
    !  integer :: idime, jdime 
    !  !m22det =   A(1,1)*A(2,2) - A(1,2)*A(2,1)
    !  do idime=1,2
    !    write(*,"(f10.5, 1x, f10.5)" )( xjacm(idime,jdime) ,jdime=1,2)
    !  end do
    !  
    !  djacb  = xjacm(1,1)*xjacm(2,2) - xjacm(1,2)*xjacm(2,1)
    !  print"(A8,f15.5)", 'det Jaco: ',djacb 
    !  if(djacb.lt.1.e-8)then
    !    write(*,*)'Element with non-positive Jacobian'
    !  end if
    ! 
    !  return
    !  
    !end function djacb
    !
    !function xjaci(detJ,xjacm)
    !  
    !  implicit none
    !  
    !  double precision, parameter :: EPS = 1.0E-10
    !  double precision, dimension(DimPr, DimPr), intent(in) :: xjacm
    !  double precision, intent(in)                          :: detJ
    !  
    !  
    !  if (abs(detJ) .le. EPS) then
    !    xjaci = 0.0! 0.0D0
    !    return
    !  end if
    !  
    !  xjaci(1,1)= xjacm(2,2)/detJ
    !  xjaci(2,2)= xjacm(1,1)/detJ
    !  xjaci(1,2)=-xjacm(1,2)/detJ
    !  xjaci(2,1)=-xjacm(2,1)/detJ
    !  
    !  return
    !  
    !end function xjaci
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    ! 
    subroutine DerivativesXY(Gp,InvJaco,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,dN_dxy, HesXY)
      
      implicit none
      
      double precision, dimension(DimPr,DimPr),intent(in)    :: InvJaco
      double precision, dimension(nne,totGp), intent(in)     :: dN_dxi, dN_deta
      double precision, dimension(nne,totGp), intent(in)     :: Hes_xixi, Hes_xieta, Hes_etaeta 
      integer, intent(in)                                    :: Gp !Variable en el lazo principal
      double precision, dimension(3,nne)                     :: Hesxieta
      double precision, dimension(2,nne)                     :: derst
      double precision, dimension(1,nne)                     :: Nxi, Neta
      double precision, dimension(1,nne)                     :: Hes_1, Hes_2, Hes_3 
      integer                                                :: idime, inode, jdime
      double precision, dimension(DimPr,nne), intent(out)    :: dN_dxy
      double precision, dimension(DimPr+1,nne), intent(out)  :: HesXY
      
      
      Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)
      Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)
      
      Hes_1 = spread(Hes_xixi(:,Gp),dim = 1, ncopies= 1)
      Hes_2 = spread(Hes_xieta(:,Gp),dim = 1, ncopies= 1)
      Hes_3 = spread(Hes_etaeta(:,Gp),dim = 1, ncopies= 1)
     
      do inode = 1, nne
        Hesxieta(1,inode) = Hes_1(1,inode)
        Hesxieta(2,inode) = Hes_2(1,inode)
        Hesxieta(3,inode) = Hes_3(1,inode)
      end do
      
      
      do idime=1,DimPr
        do inode=1,nne
          dN_dxy(idime,inode) = 0.0
          derst(1,inode) = Nxi(1,inode)
          derst(2,inode) = Neta(1,inode)
          do jdime=1,2
            dN_dxy(idime,inode) = dN_dxy(idime,inode) + InvJaco(idime,jdime) * derst(jdime,inode)
          end do
        end do
      end do
      
      
      !The Hessian matrix
      HesXY = 0.0
      do inode=1,nne
        HesXY(1,inode) = InvJaco(1,1)*InvJaco(1,1)*Hesxieta(1,inode)+&
          2.0*InvJaco(1,1)*InvJaco(2,1)*Hesxieta(2,inode)+InvJaco(2,1)*InvJaco(2,1)*Hesxieta(3,inode)
        
        HesXY(2,inode) = InvJaco(1,1)*InvJaco(1,2)*Hesxieta(1,inode) + (InvJaco(1,1)*InvJaco(2,2)+&
          InvJaco(2,1)*InvJaco(1,2)) * Hesxieta(2,inode) + InvJaco(2,1)*InvJaco(2,2)*Hesxieta(3,inode)
        
        HesXY(3,inode) = InvJaco(2,2)*InvJaco(2,2)*Hesxieta(3,inode)+&
          2.0*InvJaco(2,2)*InvJaco(1,2)*Hesxieta(2,inode)+InvJaco(1,2)*InvJaco(1,2)*Hesxieta(1,inode)
      end do
      
    end subroutine DerivativesXY 
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
    subroutine invmtx(a,deter,b)
      
      ! This routine inverts a square matrix A -> Mat(ndofn,ndofn). The
      ! inverse is stored in B. Its determinant is DETER
      
      implicit none
      
      double precision, intent(in) :: a(3,3)
      double precision :: deter, t1, t2, t3, denom
      double precision, intent(out) :: b(3,3)
      
      !nvers of a 1*1 matrix
      
      if(ndofn.eq.1) then
        deter=a(1,1)
        if(deter.eq.0.0) return
        b(1,1) = 1.0/a(1,1)
        return
      endif
      
      !invers of a 2*2 matrix
      
      if(ndofn.eq.2) then
        deter=a(1,1)*a(2,2)-a(2,1)*a(1,2)
        if(deter.eq.0.) return
        denom=1.0/deter
        b(1,1) = a(2,2)*denom
        b(2,2) = a(1,1)*denom
        b(2,1) =-a(2,1)*denom
        b(1,2) =-a(1,2)*denom
        return
      endif
      
      !inverse of a 3*3 matrix
      
      if(ndofn.eq.3) then
        t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
        t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
        t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
        deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
        if(deter.eq.0.0) return
        denom = 1./deter
        b(1,1) = t1*denom
        b(2,1) = t2*denom
        b(3,1) = t3*denom
        b(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
        b(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
        b(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
        b(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
        b(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
        b(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom
        return
      endif
      
    end subroutine invmtx 
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
    subroutine sqrtma(mainp,maout)
      
      ! Square root of matrix. In the case NDOFN = 3, it is assumed that this
      ! matrix has the form diag(A,a), where A is a 2 x 2 matrix.
      
      implicit none
      
      double precision, intent(inout) :: mainp(3,3)
      integer :: i
      double precision, intent(out) :: maout(3,3)
      
      do i = 1,2
        mainp(i,i) = abs(mainp(i,i))
      end do
      
      call sqrtm2(mainp,maout)
      
      if(ndofn.eq.3)then
        maout(3,3) = sqrt(abs(mainp(3,3)))
      endif
      
    end subroutine sqrtma 
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    !
    subroutine sqrtm2(mainp,maout)
      
      ! Square root of a 2 x 2 matrix (from a 3 x 3 matrix)
      
      implicit none
      double precision, intent(in) :: mainp(3,3)
      double precision :: a, b, c, d, aux1, aux2, vap1, vap2, det, sq1, sq2
      double precision, intent(out) :: maout(3,3)
      
      a = mainp(1,1)
      b = mainp(1,2)
      c = mainp(2,1)
      d = mainp(2,2)
      aux1 =  0.5*(a+d)
      aux2 = 0.25*(a-d)*(a-d) + b*c
      if(aux2.lt.0.0) then
        !call runend('SQRTMA: Non real eigenvalue in A')
        print*, 'SQRTMA: Non real eigenvalue in A'
        stop
      else if(aux2.lt.1.0e-10) then                         ! b or c = 0, a = d
        maout(1,1) = sqrt(a)
        maout(1,2) = 0.0
        maout(2,1) = 0.0
        maout(2,2) = sqrt(d)
        if(abs(b).gt.1.0e-10) then
          maout(1,2) = b/(sqrt(a) + sqrt(d))
        else if(abs(c).gt.1.0e-10) then
          maout(2,1) = c/(sqrt(a) + sqrt(d))
        end if
      else
        vap1 = aux1 + sqrt(aux2)                            ! vep1 = ( b ,vap1 -a )
        vap2 = aux1 - sqrt(aux2)                            ! vep2 = ( vap2 -d, c )
        if(abs(b)+abs(vap1-a).lt.1.0e-15) then
          sq1  = vap1
          vap1 = vap2
          vap2 = sq1
        end if
        sq1 = sqrt(vap1)
        sq2 = sqrt(vap2)
        vap1 = vap1 - a
        vap2 = vap2 - d
        det = b*c - vap1*vap2
        maout(1,1) = (b*c*sq1-vap1*vap2*sq2)/det
        maout(1,2) =  b*vap2*(-sq1+sq2)/det
        maout(2,1) =  c*vap1*( sq1-sq2)/det
        maout(2,2) =(-vap1*vap2*sq1+b*c*sq2)/det
      end if
      
    end subroutine sqrtm2
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
    subroutine gather(lnods, vecgl, veclo)
      !        gather(vecgl,veclo,lnods,ndofn,nnode)
      !   call gather(coord,elcod,lnods(1,ielem),2,nnode)
      !*****************************************************************************
      !
      !     Gather operations: Recupera los nodos elementales (veclo) del vector global vecgl
      !
      !*****************************************************************************
      !
      !veclo - vector global
      !vecgl - vector global
      
      implicit none
      
      double precision, dimension(*), intent(in) :: vecgl(*)
      integer, intent(in) ::   lnods(nne)
      integer   inode,idofn,ipoin,ievab,itotv
      double precision, dimension(nevab), intent(out) :: veclo
      
      do inode=1,nne
        ipoin=lnods(inode)
        do idofn=1,ndofn
          ievab=(inode-1)*ndofn+idofn
          itotv=(ipoin-1)*ndofn+idofn
          veclo(ievab)=vecgl(itotv)
        end do
      end do
      
    end subroutine gather 
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
    function elemSize(InvJacobian)
      
      implicit none
     
      double precision, dimension(DimPr,DimPr), intent(in) :: InvJacobian
      double precision :: hx, hy, elemSize
     ! hx    = sqrt(xjaci(1,1)**2+xjaci(2,1)**2)
     ! hy    = sqrt(xjaci(1,2)**2+xjaci(2,2)**2)
      
      hx    = sqrt(InvJacobian(1,1)**2+InvJacobian(2,1)**2)
      hy    = sqrt(InvJacobian(1,2)**2+InvJacobian(2,2)**2)
      
      elemSize = hnatu/(min(hx,hy))     !hnatu = Reference element length en mod_param
      
      return
      
    end function elemsize
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
    subroutine Galerkin(hmaxi, dvol, basis, dNdxy, EMsource, Ke, Ce, Fe)
      
      implicit none
      
      double precision, intent(in) :: basis(nne), EMsource(ndofn), dNdxy(DimPr,nne)
      double precision, intent(in) :: dvol, hmaxi
      integer :: inode, idofn, ievab, jevab, jnode, jdofn, i, j
      double precision ::  diff, convec, reac, cpcty, cte
      double precision, intent(in out) :: Ke(nevab,nevab), Fe(nevab), Ce(nevab,nevab)
      
      ievab=0
      do inode=1,nne
        do idofn=1,ndofn
          ievab=ievab+1
          jevab=0
          do jnode=1,nne
            do jdofn=1,ndofn
              jevab=jevab+1
              diff=0.0
              do i=1,2
                do j=1,2   
                  !if(kstab.eq.6)call param_stab(idofn, jdofn, i, j, hmaxi, cte) !conductivity tensor
                  !if(kstab.eq.6)diff = diff+ dNdxy(i,inode) * cte * difma(idofn,jdofn,i,j)* dNdxy(j,jnode)
                  diff = diff+ dNdxy(i,inode) * difma(idofn,jdofn,i,j)* dNdxy(j,jnode)
                end do
              end do
              convec=0.0
              do i=1,2
                !print*,conma(idofn,jdofn,i)
                convec = convec + basis(inode) * conma(idofn,jdofn,i) * dNdxy(i,jnode)
              end do
              reac = basis(inode) * reama(idofn,jdofn) * basis(jnode)
              cpcty = basis(inode) * basis(jnode)
              Ke(ievab,jevab) = Ke(ievab,jevab) + (diff + convec + reac) * dvol
              Ce(ievab,jevab) = Ce(ievab,jevab) + cpcty * dvol     !element Capacity (Mass) matrix
            end do
          end do
          !Fe(ievab) = Fe(ievab) + basis(inode) * EMsource(idofn) * dvol
          Fe(ievab) = Fe(ievab) + basis(inode) * force(idofn) * dvol
        end do
      end do
      
    end subroutine Galerkin
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine Galerkin_Init_Cond(dvol, basis, u0_cond, C0e, u0e)
     
      !Esta rutina proyecta la condicion inicial al dominio de elementos mediante Galerkin        
      implicit none
      
      double precision, intent(in) :: basis(nne), u0_cond(nne)
      double precision, intent(in) :: dvol
      !real, intent(in)             :: u0_cond
      integer :: inode, idofn, ievab, jevab, jnode, jdofn
      double precision :: cpcty
      
      double precision, intent(out) :: u0e(nevab), C0e(nevab,nevab)
      ievab=0
      do inode=1,nne
        do idofn=1,ndofn
          ievab=ievab+1
          jevab=0
          do jnode=1,nne
            do jdofn=1,ndofn
              jevab=jevab+1
              cpcty = basis(inode) * basis(jnode)
              C0e(ievab,jevab) = C0e(ievab,jevab) + cpcty * dvol       !element Capacity (Mass) matrix
            end do
          end do
          u0e(ievab) = u0e(ievab) + basis(inode) * u0_cond(inode) * dvol
          !u0e(ievab) = u0e(ievab) + basis(inode) * u0_cond * dvol
        end do
      end do
      
    end subroutine Galerkin_Init_Cond 
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
    subroutine pertur( idofn, jdofn, workm, derxy, basis, pertu )
      
      !***************************************************************************
      !
      ! Perturbation of the test function according to the type of method:
      !
      !   SUPG    :                 A_i   V,i           (kstab=1)
      !   GLS     : -(K_ij V,j),i + A_i   V,i + S   V   (kstab=2)
      !   SGS, TG :  (K_ij V,j),i + A^t_i V,i - S^t V   (kstab=3,5)
      !   CG      :            diag(A_i)  V,i           (kstab=4)
      !***************************************************************************
      
      implicit none
      
      double precision, intent(in)     :: workm(2,2),derxy(2),basis
      integer                          :: idofn, jdofn, k, l
      double precision                 :: prod1, prod2, prod3
      double precision, intent(in out) :: pertu
      
      
      ! SUPG
      if(kstab.eq.1) then
        prod1=0.0
        do k=1,2
          prod1=prod1+conma(jdofn,idofn,k)*derxy(k)
        end do
        pertu=prod1
        
        ! galerkin least squares
      else if(kstab.eq.2)then
        prod1=0.0
        do k=1,2
          prod1=prod1+conma(jdofn,idofn,k)*derxy(k)
        end do
        prod2=0.0
        do k=1,2
          do l=1,2
            !call param_stab(jdofn,idofn,k,l,cte1)
            !prod2=prod2+cte1*difma(jdofn,idofn,k,l)*workm(k,l)
            prod2=prod2+difma(jdofn,idofn,k,l)*workm(k,l)
          end do
        end do
        prod3=reama(jdofn,idofn)*basis
        pertu=-prod2+prod1+prod3
        
        ! subgrid scale & taylor galerkin
      else if((kstab.eq.3).or.(kstab.eq.5)) then
        prod1=0.0
        do k=1,2
          prod1=prod1+conma(idofn,jdofn,k)*derxy(k)
        end do
        prod2=0.0
        do k=1,2
          do l=1,2
            !call param_stab(idofn,jdofn,k,l,cte2)
            !prod2=prod2+cte2*difma(idofn,jdofn,k,l)*workm(k,l)
            prod2=prod2+difma(idofn,jdofn,k,l)*workm(k,l)
          end do
        end do
        prod3=reama(idofn,jdofn)*basis
        pertu=prod2+prod1-prod3
        
        ! characteristic galerkin
      else if(kstab.eq.4) then
        prod1=0.0
        if(idofn.eq.jdofn) then
          do k=1,2
            prod1=prod1+conma(jdofn,idofn,k)*derxy(k)
          end do
        end if
        pertu=prod1
      end if
      
    end subroutine pertur
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine TauMat(hmaxi,tauma)
      
      !Matrix of intrinsic time scales, computed as
      ! TAU = PATAU * [ 4 K / h^2 + 2 A / h + S ]^{-1}
     
      implicit none
      
      double precision, intent(in) :: hmaxi
      integer :: i, j, k
      double precision :: chadi(3,3), chaco(3,3), chare(3,3), tauin(3,3)
      double precision :: a, b, c, tau, det                     !hnatu -> declarado en parameters
      double precision, intent(out) :: tauma(3,3)               !ndofn -> en parameters
      
      tauma = 0.0
      if(kstab.eq.0) return
      tauin = 0.0
      chaco = 0.0
      chadi = 0.0
      
      !  Characteristic convection matrix: A = sqrt | A_i A_i |
      do i=1,ndofn
        do j=1,ndofn
          chaco(i,j)=0.0
          do k=1,ndofn
            chaco(i,j) = chaco(i,j) + conma(i,k,1)*conma(k,j,1) + conma(i,k,2)*conma(k,j,2)
          end do
        end do
      end do
      call sqrtma(chaco,chaco)
      !  Characteristic diffusion matrix: K = sqrt( K_ij K_ij )
      do i=1,ndofn
        do j=1,ndofn
          chadi(i,j)=0.0
          do k=1,ndofn
            chadi(i,j) = chadi(i,j) + difma(i,k,1,1) * difma(k,j,1,1) + &
              &difma(i,k,1,2)*difma(k,j,1,2)*2.0 + difma(i,k,2,1)*difma(k,j,2,1)*2.0 + &
              &difma(i,k,2,2)*difma(k,j,2,2)
          end do
        end do
      end do
      
      call sqrtma(chadi,chadi)
      
      !  Characteristic reaction matrix: S = | S |
      do i=1,ndofn
        do j=1,ndofn
          chare(i,j)=0.0
          do k=1,ndofn
            chare(i,j)=chare(i,j) + reama(i,k) * reama(k,j)
          end do
        end do
      end do
      call sqrtma(chare,chare)
      
      ! Invers of the matrix of characteristic times
      do i=1,ndofn
        do j=1,ndofn
          tauin(i,j) = 4.0*chadi(i,j)/(hmaxi*hmaxi) + 2.0*chaco(i,j)/hmaxi + chare(i,j)
        end do
      end do
      
      ! Matrix tau, corresponding to:
      ! KTAUM = 0: T = t I, where t is the minimum of all the admissible tau's
      !       = 1: T = diag(t1,t2,t3), ti is the minimum of the admissible tau's for the i-th row (equation)
      !       = 2: T = [ 4 K / h^2 + 2 A / h + S ]^{-1}
      
      if(ktaum.eq.0) then
        tau = 0.0
        do i=1,ndofn
          do j=1,ndofn
            tau = max(tau,abs(tauin(i,j)))
          end do
        end do
        tau = patau/tau
        do i=1,ndofn
          tauma(i,i) = tau
        end do
      
      else if(ktaum.eq.1) then
        a = 0.0
        b = 0.0
        c = 0.0
        do j=1,ndofn
          a = max(a,abs(tauin(    1,j)))
          b = max(b,abs(tauin(    2,j)))
          c = max(c,abs(tauin(ndofn,j)))
        end do
        a = patau/a
        b = patau/b
        c = patau/c
        tauma(    1,    1) = a
        tauma(    2,    2) = b
        tauma(ndofn,ndofn) = c
        
      else if(ktaum.eq.2) then
        call invmtx(tauin,det,tauma)
        do i = 1,ndofn
          do j = 1,ndofn
            tauma(i,j) = tauma(i,j)*patau
          end do
        end do
        tauma(ndofn,ndofn) = 0.0
       
      else if(ktaum.eq.3) then
        !call param_stab(1,1,1,1,cte)
        !a = 1.0/(patau*cte*difma(1,1,1,1) /(hmaxi*hmaxi) + reama(1,1))
        a = 1.0/(patau*difma(1,1,1,1) /(hmaxi*hmaxi) + reama(1,1))
        tauma(1,1) = a
        tauma(2,2) = a
        a = (hmaxi*hmaxi*hmaxi*hmaxi)/(patau*patau)
        a = a*(patau/(hmaxi*hmaxi*reama(1,1)) + 1.0d0/(difma(1,1,1,1))) !cte*(difma(1,1,1,1)))
        tauma(3,3) = a
      end if
      
    end subroutine TauMat
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    !
    subroutine Stabilization(dvolu, basis, derxy,HesXY,source, tauma, Ke,Fe)
      !subroutine Stabilization(dvolu, basis, derxy,HesXY,tauma,Ke,Fe,pertu,workm,resid)
      
      ! Contribution to the system matrix and RHS from the stabilization term
      
      implicit none
      
      double precision, intent(in)  :: basis(nne), derxy(DimPr,nne), HesXY(3,nne), tauma(3,3)
      double precision, intent(in)  :: dvolu, source(ndofn)
      double precision              :: pertu(nevab,ndofn), workm(2,2),  resid(ndofn,nevab)
      double precision              :: prod1, prod2, prod3
      integer                       :: ievab, inode, idofn, jdofn, jevab, jnode, k, l
      double precision, intent(out) :: Ke(nevab,nevab), Fe(nevab)
      
      ! integer :: nnode,ndofn,nevab,kstab,n_ini
      !difma(3,3,2,2), conma(3,3,2), reama(3,3), force(3)
      !common/proper/difma,conma,reama,force
      
      ievab = 0
      !print*, 'in'
      pertu =  0.0
      
      do inode=1,nne
        workm(1,1)=HesXY(1,inode)
        workm(2,2)=HesXY(3,inode)
        workm(1,2)=HesXY(2,inode)
        workm(2,1)=workm(1,2)
        do idofn=1,ndofn
          
          ievab=ievab+1
          do jdofn=1,ndofn
            prod1 = reama(jdofn,idofn)*basis(inode)
            prod2=0.0
            do k=1,2
              prod2 = prod2 + conma(jdofn,idofn,k)*derxy(k,inode)
            end do
            prod3=0.0
            do k=1,2
              do l=1,2
                !call param_stab(jdofn,idofn,k,l,cte)
                !prod3 = prod3 + cte*difma(jdofn,idofn,k,l)*workm(k,l)
                prod3 = prod3 + difma(jdofn,idofn,k,l)*workm(k,l)
              end do
            end do
            
            resid(jdofn,ievab) = prod1 + prod2 - prod3
            call pertur( idofn, jdofn, workm, derxy(1,inode), basis(inode), pertu(ievab,jdofn) )
          end do
        end do
      end do
      
      ievab=0
      do inode=1,nne
        do idofn=1,ndofn
          ievab=ievab+1
          jevab=0
          do jnode=1,nne
            do jdofn=1,ndofn
              jevab=jevab+1
              prod1=0.0
              do k=1,ndofn
                do l=1,ndofn
                  prod1 = prod1 + pertu(ievab,k)*tauma(k,l)*resid(l,jevab)
                end do
              end do
              Ke(ievab,jevab) = Ke(ievab,jevab) + prod1 * dvolu
            end do
          end do
          
          prod1=0.0
          do k=1,ndofn
            do l=1,ndofn
              prod1 = prod1 + pertu(ievab,k) * tauma(k,l) * force(l) * source(l)
            end do
          end do
          Fe(ievab) = Fe(ievab) + prod1 * dvolu
        end do
      end do
      
    end subroutine Stabilization
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
    subroutine VinculBVs(condition,  BVs, nofix, ifpre, presc )
      
      implicit none
      
      !integer, intent(in)      :: nBvs, nBVscol ya no se ponen estan en el modulo parameter y se comunica el valor
      integer, intent(in) :: condition( nBvs, nBVscol-ndofn)
      double precision, intent(in) :: BVs( nBvs, ndofn)
      integer             :: i, j
      double precision, intent(out) :: presc(ndofn,nBVs)
      integer, intent(out)          :: ifpre(ndofn,nBVs)
      integer, intent(out)          :: nofix(nBVs)
      
      
      select case(ndofn)
        case(1)
          do i =1,ndofn
            do j=1,nBVs
              nofix(j)   = condition(j,1)
              ifpre(i,j) = condition(j,2) !El llenado de ifpre sera por grado de libertad
              presc(i,j) = Bvs(j,1)
            end do
          end do
          
        case(2)
          do i =1,ndofn
            do j=1,nBVs
              nofix(j)   = condition(j,1)
              ifpre(i,j) = condition(j,i+1) !El llenado de ifpre sera por grado de libertad
              presc(i,j) = Bvs(j,i)
            end do
          end do
         
        case(3)
          do i =1,ndofn
            do j=1,nBVs
              nofix(j)   = condition(j,1)
              ifpre(i,j) = condition(j,i+1) !El llenado de ifpre sera por grado de libertad
              presc(i,j) = Bvs(j,i)
            end do
          end do
          
        case DEFAULT
          write(*,*) 'Exceeded DoF'
        end select
    end subroutine VinculBVs
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine BandWidth( )
      !this routine computes the bandwidth 
      
      implicit none
      
      integer :: iband, ielem, inode, ipoin, jnode, jpoin, nband       ! , C, D
      
      iband=0
      do ielem =1, nelem
        do inode = 1, nne
          ipoin = lnods(ielem,inode) !Este +1 es para que comience en los nodos (columna 2) y no del numeor de elemento
          do jnode = 1, nne
            jpoin = lnods(ielem,jnode)
            iband = max(iband,abs(jpoin-ipoin))
          end do
        end do
        !if (iband.gt.nband) nband=iband !nband pasa como variable global por que se usa en  ApplyBVal y otros
        nband = iband
      end do
      nband=(nband+1)*ndofn-1
      
      upban  = nband
      lowban = nband
      totban = lowban + upban + 1
      ldAKban= 2*lowban + upban + 1   
      
      print*, ' '
      print*,'!================ Bandwidth Info ==============!'
      write(*,"(A15,9X,I6,1X,A9)")' - UpBand:      ', upban,'   '
      write(*,"(A15,9X,I6,1X,A9)")' - LowBand:     ', lowban,'  '
      write(*,"(A15,9X,I6,1X,A9)")' - TotBand:     ', totban,'  '
      write(*,"(A15,9X,I6,1X,A9)")' - ledimAK:     ', ldAKban,' '
      
    end subroutine BandWidth
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine param_stab(idofn, jdofn, i, j, hmaxi, coeff)       
      !***********************************************************!
      !                                                           !
      ! Subroutine which check the dofn, x and y position in the  !
      ! diffusion tensor and take the coefficient to multiply     !
      ! the term of the PDE to its corresponding coefficient.     !
      !                                                           !
      ! coeficients:                                              !
      !             Cuλ(h^2/ell^2),  λ  and   ell^2/λ             !
      !                                                           !
      !  λ = 1/µ  = reluctivity of the medium                     !
      !                                                           !
      !  h = 2^-i = hmaxi                                         !
      !                                                           !
      !***********************************************************!
      
      implicit none
      
      integer, intent(in) :: idofn, jdofn, i, j
      double precision, intent(in) :: hmaxi
      double precision, intent(out) :: coeff
      
      
      !print*, 'getting into param_stab'
      
      coeff = 0.0 
      !print*, 'helem^2', hmaxi**2
      if(idofn.eq.1)then
        if(jdofn.eq.1)then                      !difma(idofn,jdofn,i,j)
          if(i==1 .and. j==1)then                      !difma(1,1,1,1)
            coeff = Cu*lambda*(hmaxi**2/ell**2)
            !print*,'Beta', coeff
          end if
         
          if(i==2 .and. j==2)then                      !difma(1,1,2,2)
            coeff = lambda
          endif
          
        elseif(jdofn==2)then                           !difma(1,2,1,2)
          if(i==1 .and. j==2)then
            coeff = Cu*lambda*(hmaxi**2/ell**2)
            !print*,'Beta', coeff
          end if
          
          if(i==2.and.j==1)then                        !difma(1,2,2,1)
            coeff = lambda
          end if
        end if
        
      elseif(idofn==2)then
        if(jdofn.eq.1)then
          if(i==1 .and. j==2)then                      !difma(2,1,1,2)
            coeff = lambda
          end if
          
          if(i==2 .and. j==1)then                      !difma(2,1,2,1)
            coeff = Cu*lambda*(hmaxi**2/ell**2)
            !print*,'Beta', coeff
          endif
         
        elseif(jdofn==2)then
          if(i==1 .and. j==1)then
            coeff = lambda
          end if
          
          if(i==2.and.j==2)then                        !difma(2,2,2,2)
            coeff = Cu*lambda*(hmaxi**2/ell**2)
            !print*,'gamma', coeff
          end if
        end if
        
      elseif(idofn==3 .and. jdofn==3)then              !difma(3,3,1,1) or difma(3,3,2,2)
        if( i==j )then
          coeff = ell**2/lambda
            !print*,'gamma', coeff
        endif
        
      else
        continue
      end if
      !close(10)
      
      15 continue 
      !9 format(A20,A6,I1,A1,I1,A1,I1,A1,I1,A1,I1,A1)
      !Next lines are to taste the 
      !print*, 'hmaxi,', h
      !print*, 'Cu µ h^2/ell^2', Cu*lambda*(h**2/ell**2)
      !print*, 'ell^2/µ', ell**2/lambda
      
    end subroutine param_stab
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine GlobalSystem(N, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, A_C, A_K, A_F)
      
      !use sourceTerm
      
      implicit none
      
      double precision, dimension(nne,TotGp), intent(in) :: N, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp), intent(in) :: hes_xixi, hes_xieta, hes_etaeta
      double precision, dimension(ndofn)        :: EMsource
      double precision, dimension(nne)          :: basis, xi_cor, yi_cor
      double precision, dimension(DimPr,nne)    :: dN_dxy
      double precision, dimension(3,nne)        :: HesXY
      double precision, dimension(DimPr, dimPr) :: Jaco, Jinv
      double precision, dimension(nevab, nevab) :: Ke, Ce
      double precision, dimension(nevab)        :: Fe
      double precision, dimension(3,3)          :: tauma
      double precision, dimension(nne,DimPr)    :: element_nodes
      integer, dimension(nne)                   :: nodeIDmap
      double precision                          :: dvol, hmaxi, detJ!, aaa
      
      integer                                   :: igaus, ibase, ielem
      double precision, allocatable, dimension(:,:), intent(out)  :: A_K, A_C, A_F
      
      allocate(A_K(ldAKban,ntotv), A_C(ldAKban,ntotv), A_F(ntotv, 1))
      
      !duda: Fe se declaró como a(n) y en la rutina assembleF como a(n,1)
      !      pero compila y ejecuta bien. ¿Poooor?
      
      A_K = 0.0
      A_F = 0.0
      hmaxi = 0.0
      
      do ielem = 1, nelem 
        !gather
        Ke = 0.0    !Esto es amate
        Fe = 0.0    !Fe(nevab)
        Ce = 0.0    !elemental capacity matrix (not used in static case)
        call SetElementNodes(ielem, element_nodes, nodeIDmap, xi_cor, yi_cor)
        !do-loop: compute element stiffness matrix Ke
        do igaus = 1, TotGp
        !print*, 'gauss point:', igaus
          call Jacobian( element_nodes, dN_dxi, dN_deta, igaus ,Jaco, detJ, Jinv)
          !Jaco = xjacm(element_nodes, dN_dxi, dN_deta, igaus)
          !detJ = djacb(Jaco)
          !Jinv = xjaci(detJ,Jaco)
          dvol = detJ *  weigp(igaus,1)
          call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, dN_dxy, HesXY)
          hmaxi = elemSize(Jinv)
          do ibase = 1, nne
            basis(ibase) = N(ibase,igaus)
          end do
          call TauMat(hmaxi,tauma)
          
          !call source_term(ielem, basis, xi_cor, yi_cor, EMsource)
          call Galerkin(hmaxi, dvol, basis, dN_dxy, EMsource, Ke, Ce, Fe) !amate lo llame Ke
          !call Galerkin(dvol, basis, dN_dxy, EMsource, Ke, Ce, Fe) !amate lo llame Ke
          !call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, Ke, Fe, pertu,workm,resid)
          if(kstab.ne.6.or.kstab.ne.0)call Stabilization(dvol, basis, dN_dxy, HesXY, EMsource, tauma, Ke, Fe)
        end do
        !stop
        
        call Assemb_Glob_Mat(nodeIDmap, Ce, A_C)     !Assemble Global Conductivity Matrix K
        call Assemb_Glob_Mat(nodeIDmap, Ke, A_K)     !Assemble Global Conductivity Matrix K
        call Assemb_Glob_Vec(nodeIDmap, Fe, A_F)     !Assemble Global Source vector F
        
      end do
      
      !aaa = maxval(coord(1,:))*2**(-i_exp) 
      !if(aaa.ne.hmaxi) write(*,'(A)') '> > >Element size does not match'
      
    end subroutine GlobalSystem
    
    !subroutine AssembleK(K, ke, node_id_map, ndDOF)

    !  implicit none
    !  real(8), dimension(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes),intent(in out)  :: K 
    !  !  Global Stiffnes matrix debe 
    !  !  llevar inout por que entra como variable (IN) 
    !  !  pero en esta funcion se modifica (out)
    !  real(8), dimension(2*nUne, 2*nUne), intent(in)   :: ke
    !  integer, dimension(nUne,1), intent(in)           :: node_id_map
    !  integer, intent(in)                              :: ndDOF 
    !  integer :: i, j, row_node, row, col_node, col !nodal Degrees of Freedom
    !  
    !  do i = 1, nUne
    !    row_node = node_id_map(i,1)
    !    row = ndDOF*row_node - (ndDOF-1)
    !    
    !    do j = 1, nUne
    !      col_node = node_id_map(j,1)
    !      col = ndDOF*col_node - (ndDOF-1)
    !      K(row:row+ndDOF-1, col:col+ndDOF-1) =  K(row:row+ndDOF-1, col:col+ndDOF-1) + &
    !      ke((i-1)*ndDOF+1:i*ndDOF,(j-1)*ndDOF+1:j*ndDOF)
    !    enddo
    !    
    !  enddo
    !  
    !  return
    !  
    !end subroutine AssembleK
   
    !subroutine ApplyBoundCond( NoBV, Fbcsvp, A_K, Sv )
    !  ! - - - - - - - - - - * * * * * * * * * * - - - - - - - 
    !  ! Set velocity (u) and preasure (p) boundary condition by penalty method
    !  ! - - - - - - - - - - * * * * * * * * * * - - - - - - - - - -
    !  implicit none
    !                      !Dof
    !  integer , dimension(NoBV,Dof), intent(in) :: Fbcsvp
    !  double precision, dimension(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes),intent(in out) :: A_K  !Global Stiffnes matrix
    !  double precision, dimension(2*n_nodes+n_pnodes, 1), intent(in out) :: Sv
    !  double precision :: param, coeff
    !  integer          :: preasure_row, NoBV, i, component, node_id, pnode_id
    !  
    !  !Esencialmente la siguiente instruccion hace: A_K(1*2-1,:) = A_K(1,:) Es decir, obtene el valor maximo de
    !  !la primera fila de la matriz global K (A_K). No le veo el caso pero lo dejamos asi.
    !  param = maxval(A_K(int(Fbcsvp(1,1))*2-1,:))
    !  coeff = abs(param) * 1.0E7

    !  print*, 'param', param
    !  print*, 'coeff', coeff
    !  
    !  preasure_row = 2*n_nodes

    !  do i =1, NoBV
    !    !print*, 'iteration', i
    !    node_id   = Fbcsvp(i,1) !se pone este int() pq la 1a y 2a col de Fbcsvp esta leida como integer pero 
    !    !print*, 'node_id', node_id
    !    component = Fbcsvp(i,2)!la matriz completa esta declarada como real en esta subroutina y en el main.
    !    !print*, 'component', component
    !    !print*, shape(Fbcsvp)
    !    !print*, Fbcsvp(i,:)
    !    if ( component .le. 2 ) then
    !      !print*, 'component of Boundary value', component
    !      !print*,'La pausa', 2*node_id-2+component,' ', 2*node_id-2 +component
    !      !read(*,*)
    !      A_K(2*node_id-2+component, 2*node_id-2 +component) = coeff
    !      Sv( 2*node_id-2+component, 1) = Fbcsvp(i,3)*coeff 
    !    else                                                     
    !      pnode_id = pnodes(node_id,2)
    !      !print*, 'pnode_id', pnode_id 
    !      !print*, 'preasure_row', preasure_row
    !      A_K(preasure_row+pnode_id, preasure_row + pnode_id) = coeff
    !      !print*, preasure_row+pnode_id, preasure_row + pnode_id
    !      Sv(preasure_row+pnode_id,1) = Fbcsvp(i,3)*coeff 
    !      !el tres es por que en la columna 3 esta el valor de la condicon de forntera
    !    end if
    !  end do

    !end subroutine ApplyBoundCond
   
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine Assemb_Glob_Mat(lnods,Ke,A_K)
      !subroutine Assemb_Glob_Mat(ielem,lnods,Ke,A_K)
      !*****************************************************************************
      !
      !    Fa l'assembly de les matrius de CDR de cada elemento en la matriu global
      !
      !*****************************************************************************

      implicit none
      !common /contrl/ npoin,nelem,nmats,nvfix,nload,nband,ntotv
      double precision, intent(in) :: Ke(nevab,nevab)
      integer, intent(in) :: lnods(nne)
      integer :: inode, ipoin, idofn, ievab, itotv, jnode, jpoin, jdofn, jevab, jtotv, jband !, i,j,k,l
      double precision, intent(in out) :: A_K(ldAKban,ntotv)
     !                                        ldAKban= 2*lowban + upban + 1
     !                                          totban = lowban + upban + 1
      do inode=1,nne    !nne = number of node in the element
        ipoin=lnods(inode)
        !print*,'ipoin',ipoin
        do idofn=1,ndofn
          ievab=(inode-1)*ndofn+idofn
          itotv=(ipoin-1)*ndofn+idofn
          do jnode=1,nne
            jpoin=lnods(jnode)
            do jdofn=1,ndofn
              jevab=(jnode-1)*ndofn+jdofn
              jtotv=(jpoin-1)*ndofn+jdofn

              !jband=jtotv-itotv+1      ------> Original de retpla
              !A_K(jband,itotv)=A_K(jband,itotv) + ...

              !jband= upban +1 +itotv-jtotv    !Algoritmo de recuperacion para LAPACK
              jband = itotv-jtotv + totban
              if (jband.ge.1)then
                A_K(jband,jtotv)=A_K(jband,jtotv)+Ke(ievab,jevab)
              endif
            end do
          end do
        end do
      end do

      return
    end subroutine Assemb_Glob_Mat
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
    subroutine Assemb_Glob_Vec( nodeIDmap, fe, F_global)
      
      implicit none
      
      double precision, dimension(nevab,1), intent(in)   :: fe  !nevab = Num of element variable
      integer, dimension(nne), intent(in)                :: nodeIDmap
      integer :: i, rowNode, row
      double precision, dimension(ntotv,1),intent(inout):: F_global
      
      do i = 1, nne
        rowNode = nodeIDmap(i)                !global id for row nodes
        row =  ndofn * rowNode - (ndofn-1)    !row number in the global F
        F_global(row:row+ndofn-1,1) =  F_global(row:row+ndofn-1,1) + fe( (i-1)*ndofn+1:i*ndofn, 1)
        
      end do
      
    end subroutine Assemb_Glob_Vec 
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    !
    subroutine ApplyBVs(nofix,ifpre,presc,rigid,gload)
      !                                   ,A_K ,A_F
      !        vincul(rigid,gload,treac,nofix,ifpre,presc)
      !*****************************************************************************
      !
      !   Imposa les condicions de contorn
      !
      !*****************************************************************************
      
      implicit none !real*8(a-h,o-z)
      
      !Agregar un common a BVs para guardar nBVs y nBVscol asi como nband.
      !Ya se hizo y se uso el modulo mod_param para guardarlos ahi y el valor se comparte
      double precision,intent(in) :: presc(ndofn,nBVs)
      integer, intent(in)         :: nofix(nBVs), ifpre(ndofn,nBVs)
      !                                    nvfix             ,nvfix
      !common /contrl/ npoin,nelem,nmats,nvfix,nload,nband,ntotv
      integer :: ivfix, idofn, itotv, jband, nvfix, ktotv, ipoin
      double precision,intent(inout)  :: rigid(ldAKban,ntotv), gload(ntotv)
      
      nvfix = nBVs
      
      !***  Llac sobre els nodes coaccionats (Lazo sobre los nodos preescritos)
      do ivfix=1,nvfix
        ipoin=nofix(ivfix)
        do idofn=1,ndofn !3
          if (ifpre(idofn,ivfix).eq.1) then
            
            !***  Ca de grau de llibertat prescrit
            itotv=(ipoin-1)*ndofn+idofn
            do ktotv=1,ntotv
              jband = ktotv-itotv + totban
              if((lowban +1.le. jband).and.(jband.le.ldAKban))then
                gload(ktotv) = gload(ktotv)-rigid(jband,itotv)*presc(idofn,ivfix)
                rigid(jband,itotv)=0!7.777
              end if
            end do
            do ktotv= 1, ntotv
              jband = itotv-ktotv + totban
              if((lowban +1.le. jband).and.(jband.le.ldAKban))then
                rigid(jband,ktotv)=0!5.555 !columna K
              end if
            end do
            rigid(totban,itotv)=1.0
            gload(itotv)=presc(idofn,ivfix)
          end if
        end do
      end do
      
      return
      
    end subroutine ApplyBVs 
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    
    !subroutine source_term(basis, xi_cor, yi_cor, source)
    !  
    !  implicit none     !interpolating X and Y
    !  
    !  !***********************************************************!
    !  !The source term is given by:                               !
    !  !                                                           !
    !  !              u = (fi*psi',-fi'*psi)                       !
    !  !                                                           !
    !  !   f = Lu       ;   where L is the diferential operator    !
    !  !                                                           !
    !  !***********************************************************!
    !  
    !  double precision,dimension(nne), intent(in) :: basis, xi_cor, yi_cor
    !  double precision :: dey_dydx, dex_dy2, dey_dx2, dex_dxdy
    !  double precision :: x, y
    !  integer :: ii
    !  double precision, dimension(ndofn), intent(out)  :: source
    !  
    !  !En cada llamada de esta funcion entrara basis con diferente gauss point
    !  
    !  x = 0.0
    !  y = 0.0
    !  
    !  !x= matmul(basis,xi_cor)
    !  do ii = 1, nne 
    !    x = x + basis(ii)*xi_cor(ii)
    !    y = y + basis(ii)*yi_cor(ii)
    !    !print"(3(1x,f10.7))",  basis(ii), yi_cor(ii), basis(ii)*yi_cor(ii)
    !  end do
    !  
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N1', basis(1),' xi1', xi_cor(1), 'x1', basis(1)*xi_cor(1)
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N2', basis(2),' xi2', xi_cor(2), 'x2', basis(2)*xi_cor(2)
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N3', basis(3),' xi3', xi_cor(3), 'x3', basis(3)*xi_cor(3)
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N4', basis(4),' xi4', xi_cor(4), 'x4', basis(4)*xi_cor(4)
    !  !print"(A9,f10.5)",'x_gaus = ', x
    !  !
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N1', basis(1),' yi1', yi_cor(1), 'y1', basis(1)*yi_cor(1)
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N2', basis(2),' yi2', yi_cor(2), 'y2', basis(2)*yi_cor(2)
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N3', basis(3),' yi3', yi_cor(3), 'y3', basis(3)*yi_cor(3)
    !  !print"(A4,1x,f9.5,1x,A4,1x,f10.5,1x, A3, f10.5)",  '  N4', basis(4),' yi4', yi_cor(4), 'y4', basis(4)*yi_cor(4)
    !  !print"(A9,f10.5)",'y_gaus = ', y
    !  
    !  !Derivatives in x-direction
    !  dey_dydx = 2.0-(4.0*y)
    !  dex_dy2  = 0.0 
    !  
    !  !Derivatives in y-direction
    !  dey_dx2  = 0.0
    !  dex_dxdy = (4.0*x)-2
    !  
    !  if(ndofn.eq.3)then
    !    source(1) = lambda*(dey_dydx - dex_dy2)
    !    source(2) = lambda*(-dey_dx2  + dex_dxdy)
    !    source(3) = 0.0
    !  elseif(ndofn.eq.2)then
    !    source(2) = lambda*(-dey_dx2  + dex_dxdy)
    !  else
    !    source(1) = lambda*(dey_dydx - dex_dy2)
    !  endif
    !  
    !end subroutine source_term
    
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine Res_Matlab(solution)
      
      implicit none
      
      double precision, parameter :: pi = 4*atan(1.d0)
      !double precision, parameter :: mu=1.25663706d-6 ! 4d0*pi*1d-7
      character(len=*), parameter :: fileplace = "Res/"
      double precision, dimension(ntotv, 1), intent(in) :: solution
      
      double precision, dimension(1, ntotv) :: solution_T
      !double precision, dimension(1,nnodes) :: xcoor, ycoor
      !double precision, dimension(nnodes)   :: x, y
      !double precision, dimension(max_time) :: t
      double precision, dimension(nnodes)   :: exact_y, exact_x, exact_p, FEM_x, FEM_y, FEM_p
      double precision     :: aa, bb, cc, dd, ee, x, y
      double precision     :: sum_error, error_EM, error_p, xmax, x_FEM, y_FEM, p_FEM, uxSol, uySol, multi
      double precision     :: fi, psi, der_fi, der_psi!, errL2_x, errL2_y, errL2_p
      !double precision     :: SrcCurr, z, sigma, ds
      character(len=12)    :: extension, File_Solution
      integer              :: ipoin, ielem, inode 
      
      
      extension ="_matlab.txt"
      File_Solution ="globsolution"
      
      solution_T = transpose(solution)
      
      error_EM = 0.0
      error_p  = 0.0
      sum_error = 0.0
      exact_x = 0.0
      exact_y = 0.0
      exact_p = 0.0
      uxSol = 0.0 
      uySol = 0.0
      multi = 1E-18
      
      select case(simul)
        case(1)
          aa = (2.0/3.0)*n_val
          bb = 0.0
          cc = (n_val/3.0) - 1.0
          !cc = (n_val/3.0) - (1.0/2.0)
          dd = 0.0
          ee = 0.0 
          !write(*,*) '       FEMx','            Ex_x', '            FEMy','           Ex_y'
          do inode = 1, nnodes  
            x = coord(1,inode)
            y = coord(2,inode)
            
           !exact solution
            bb    = ((x**2 + y**2))
            
            
            if(abs(x).ne.0.01.and.abs(x).le.1.0e-4)then
              dd    = y/x
              ee    = atan(dd)
            else
              ee    = pi/2.0
            end if
              
            uxSol = aa * bb**cc * ( x*sin(aa*ee) - y*cos(aa*ee) )
            uySol = aa * bb**cc * ( y*sin(aa*ee) + x*cos(aa*ee) )
            !uxSol = aa * bb**cc * sin(aa*ee)
            !uySol = aa * bb**cc * cos(aa*ee)
            
            !FEM solution
            x_FEM = solution_T(1,ndofn*inode-2)
            y_FEM = solution_T(1,ndofn*inode-1)
            
            !write(*,"(4(f15.5,1x))") x_FEM, uxSol, y_FEM, uySol
            !error = error + (x_FEM - uxSol)**2 + (y_FEM - uySol)**2 
            error_EM = error_EM + ( (uxSol - x_FEM)**2  + (uySol - y_FEM)**2 )
            !print*, 'error', error  
            
            !Write to plotting file 
            exact_x(inode) = uxSol 
            exact_y(inode) = uySol
            
          end do
        case(2)
          
        !  ! Implements eq. (2.50) of Nabighian 1988 EM Methods (EM Theory book) (p. 175)
        !  !Nota: Esta solucion analitica determina el campo electrico en un punto especifico 
        !  !de la malla (x,y) por lo que, para usarse de comparacion se requiere ejecutar el codigo y luego en 
        !  !el post-proceso extraer en un punto determinado ux e uy para todos los tiempos simulados
        !  implicit real*8 (a-h,o-z)
        !  !!  real(8) :: x = 0.17_8
        !  !!  x = erfc(x)
        !  ! Constants
        !  ! Times array
        !  
        !  ! Define variables
        !  ! I*ds = dipole moment, is set to I*ds=1
        !  srcCurr=1d0 ! I
        !  ds=1d0
        !  sigma = 0.01
        !  x = 1.0
        !  y = 1.0
        !  z = 0.0
        !  nt = max_time
        !  do i=1,nt
        !    t(i) = i*0.25
        !  end do
        !  r_vec   = sqrt(x*x+y*y+z*z)
        !  spi = sqrt(pi)
        !  cc   = SrcCurr*ds/(4.*pi*sigma*r_vec**3)
        !  
        !  write(*,*)'      t(i)       ex       ey           ez'
        !  ! r
        !  do i=1,nt
        !    
        !    thta = sqrt(mu*sigma/(4.*t(i)))
        !    aa = 4./spi*theta**3*r_vec**3 + 6./spi*theta*r_vec
        !    arg= -theta**2*r_vec**2
        !    ee = erfc(theta*r_vec)
        !    aa = aa*exp(arg)+3.*ee
        !    bb = 4./spi*theta**3*r_vec**3 + 2./spi*theta*r_vec
        !    bb = bb*exp(arg)+ee
        !    
        !    ! geometry term
        !    ex=cc * (aa*x**2/r_vec**2 - bb)
        !    ey=cc * aa*x*y/r_vec**2
        !    ez=cc * aa*x*z/r_vec**2
        !    
        !    exact_x(i) = ex
        !    exact_y(i) = ey
        !    exact_z(i) = ez
        !    
        !    write(*,'(4(1PE14.6))')t(i), exact_x, exact_y, exact_z
        !  end do
          
        case(3) !new test of polynomial solution
          do inode = 1, nnodes  !simple Function
            x = coord(1,inode)
            y = coord(2,inode)
            
            !exact solution
            uxSol   = x**2 * (1.0-2.0*x + x**2)*(2*y**3 -3*y**2 +y)
            uySol   =-(2.0*x**3 -3.0*x**2 +x)*y**2 *(1.0-2.0*y + y**2)
            
            !FEM solution
            x_FEM = solution_T(1,ndofn*inode-2)
            y_FEM = solution_T(1,ndofn*inode-1)             
            p_FEM = solution_T(1,ndofn*inode)
           
            error_EM = error_EM + ( (uxSol - x_FEM)**2  + (uySol - y_FEM)**2 )
            error_p  = error_p  + (multi - p_FEM )**2
            
            !Write to plotting file 
            exact_x(inode) = uxSol
            exact_y(inode) = uySol
            exact_p(inode) = multi
            FEM_x(inode) = solution_T(1,ndofn*inode-2)
            FEM_y(inode) = solution_T(1,ndofn*inode-1)
            FEM_p(inode) = solution_T(1,ndofn*inode-0) 
          end do
          
        case(4)
          goto 115 
        case(5)
          print*,'No analytic solution'
        case(6)
          print*,'DC simulation No analytic solution'
      end select
      error_EM = sqrt(error_EM/nnodes)
      error_p = sqrt(error_p/nnodes)
     
      !errL2_x=norm2(exact_x - FEM_x)/norm2(exact_x)
      !errL2_y=norm2(exact_y - FEM_y)/norm2(exact_y)
      !errL2_p=norm2(exact_p - FEM_p)/norm2(exact_p)
      
      !write(*,*) '       FEMx','            Ex_x', '            FEMy','           Ex_y'
      open(unit=111, file= fileplace//File_Nodal_Vals//extension, ACTION="write", STATUS="replace")
      !open(unit=112, file= fileplace//File_Solution//extension, ACTION="write", STATUS="replace")
      
      
      if(ndofn.eq.3)then
        do ipoin = 1, nnodes  !   uh_x    uh_y    uex_x   uex_y
         
         ! write(111,906) solution_T(1,ndofn*ipoin-2),solution_T(1,ndofn*ipoin-1),&
         ! &     solution_T(1,ndofn*ipoin), exact_x(ipoin), exact_y(ipoin), exact_p(ipoin)
          
          write(111,906) solution_T(1,ndofn*ipoin-2),solution_T(1,ndofn*ipoin-1),&
          &              exact_x(ipoin), exact_y(ipoin)
          
        end do
        
      elseif(ndofn.eq.1)then 
        do ipoin = 1, nnodes  !   uh_x    uh_y    uex_x   uex_y
          write(111,906) solution(ipoin,1)
          !write(111,906) solution_T(1, ipoin), exact_x(ipoin)
        end do
      else
        print*, 'In Res_Matlab, Problem type not defined'
        stop
      end if
      
      !print*, '!====== Matlab file ======'
      write(*,"(A7,A23,A28)") ' -File ',File_Nodal_Vals//extension, 'written succesfully in Res/'
      
      close(111)
      
      open(unit=444, file= fileplace//coord_name//extension, ACTION="write", STATUS="replace")
      open(unit=333, file= fileplace//conec_name//extension, ACTION="write", STATUS="replace")
      
      do ielem=1,nelem
        write(333,902) (lnods(ielem,inode),inode=1,nne)
      end do
      
      do ipoin = 1, nnodes
        write(444,904) ipoin, coord(1,ipoin), coord(2,ipoin)
      end do
      write(*,"(A7,A12,A5,A12,A29)") ' -File ',coord_name, ' and ', conec_name, 'written succesfully in Res/ '
      print*, ' '
      
      close(333)
      close(444)

      open(unit=777, file= fileplace//error_name//extension, ACTION="write", STATUS="replace")
      write(777,"(1x,E15.5,3x, A)") error_EM,'%error in electric field'
      write(777,"(1x,E15.5,3x, A)") error_p, '%error in multiplier'
      !write(777,"(1x,E15.5,3x, A)") errL2_x,'%L2error in ex'
      !write(777,"(1x,E15.5,3x, A)") errL2_y,'%L2error in ey'
      !write(777,"(1x,E15.5,3x, A)") errL2_p,'%L2error in multiplier'
      close(777)

      xmax = maxval(coord(1,:)) !the greatest number in x column
      print*, '!============== Error Estimation ==============!'
      !write(*,"(A10,f7.5,A25,E13.5)")' -For h = ', xmax*2**(-i_exp), 'the error estimation is '
      write(*,"(A8,E13.5)")' -In u: ', error_EM
      write(*,"(A8,E13.5)")' -in p: ', error_p
      print*, ' '
     
      902 format(1x,i5,10(1x,i5))
      904 format(I7,2(3x,f9.4) ) !format for msh
      906 format(6(E15.5, 3x))
      
      115 continue
    end subroutine Res_Matlab
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine MKLfactoResult( routine_name, num )
      implicit none
      
      character(*), intent(in)  :: routine_name
      integer,      intent(in)  :: num
      character(len=34) :: text
      character(len=44) :: text2
      integer           :: val
      external          :: xerbla
      
      text  = '   *FACTORIZATION DONE WITH STATUS'
      text2 = '   *FACTORIZATION HAS BEEN COMPLETED, BUT U('
      if ( num .eq. 0 ) then
        !print*, ' '
        write(*, 101) text, num, ', THE EXECUTION IS SUCCESSFUL.'
      elseif(num .lt. 0 )then
        val = abs(num)
        print*, ' '
        write(*, 102) '    THE',val,'-TH PARAMETER HAD AN ILLEGAL VALUE.'
        call xerbla( routine_name, num )
      elseif(num .gt. 0 )then
        print*, ' '
        write(*, 103) text2,num,',',num,') = 0.0. Then U IS EXACTLY SINGULAR.'
        print"(A)",'   DIVISION BY 0 WILL OCCUR IF USE THE FACTOR U FOR SOLVING A SYSTEM OF LINEAR EQUATIONS.'
        print*, ' '
        call xerbla( routine_name, num )
        print*, ' ~ ~ ~ '
        print*, ' ~ ~ ~ Stopping the execution'
        print*, ' '
        stop
      endif
      print*, ' '

      101 format (A, 1x, I1, A)
      102 format (A, I4, A)
      103 format (A, I3, A, I3, A)
    end subroutine MKLfactoResult
    
    !subroutine MKLCondNumber(norm, ab, ipiv)
    !  implicit none
    !  
    !  double precision, dimension(ldAKban ,ntotv ), intent(in) :: ab
    !  integer, dimension(ntotv), intent(in) :: ipiv
    !  character(*), intent(in)              :: norm         !must be 1 or I
    !  !double precision          :: work(*)
    !  double precision, allocatable, dimension(:) :: work
    !  integer, dimension(ntotv) :: iwork
    !  double precision          :: rcond, val
    !  integer                   :: orderMatrix
    !  external                  :: xerbla, dgbcon
    !  integer                   :: info
    !  
    !  allocate(work(max(1,3*ntotv)))
    !  orderMatrix = ldAKban*ntotv   
    !  
    !  !val = dlangb(norm, orderMatrix, lowban, upban, ab, ldakban, work )
    !  val = 23550.9023
    !  call dgbcon( norm, orderMatrix, upban, lowban, ab, ldakban, ipiv, val, rcond, work, iwork, info )
    !  
    !  call xerbla('dgbcon',info)

    !  print*,'Reciprocal of Condition Number', rcond 
    !  print*,'Condition Number', 1./rcond 
    !  
    !  
    !  
    !end subroutine MKLCondNumber
    
    subroutine MKLsolverResult(routine_name, num )
      implicit none
      character(*), intent(in)  :: routine_name
      integer,      intent(in)  :: num
      integer :: val
      character(len=30) :: text
      character(len=35) :: text2
      character(len=30) :: text3
      external :: xerbla
      text =  '   *SYSTEM SOLVED WITH STATUS'
      text2 = '-TH PARAMETER HAD AN ILLEGAL VALUE.'
      text3 = '   *THE LEADING MINOR OF ORDER'
      if ( num .eq. 0 ) then
        write(*,101) text, num, ', THE EXECUTION IS SUCCESSFUL.'
      elseif(num .lt. 0 )then
        val = abs(num)
        write(*,102) '    THE',val, text2
        call xerbla( routine_name, num )
      elseif(num .gt. 0 )then
        print*, ' '
        write(*, 103) text3,num,' (THEREFORE THE MATRIX A ITSELF) IS NOT'
        print*,'   POSITIVE-DEFINITE. THE FACTORIZATION COULD NOT BE COMPLETED. '
        print*,'   THIS MAY INDICATE AN ERROR IN FORMING THE MATRIX A.'
        print*, ' '
        call xerbla( routine_name, num )
        print*, ' ~ ~ ~ Stopping the execution'
        print*, ' '
        stop
      endif
      101 format (A, 1x, I1, A)
      102 format (A, I3, A)
      103 format (A30, I3, A)
      !print*,' '
    end subroutine MKLsolverResult
    !
    !--------------------------------------------------------------------------------  
    !
    subroutine writeMatrix(Matrix, name1, Vector, name2)
      implicit none
      
      character(len=*), parameter    :: fileplace = "Res/"
      character(*) :: name1, name2
      integer :: i, j, mrow, ncol
      double precision, dimension(ldAKban ,ntotv ), intent(in) :: Matrix
      !double precision, dimension(ntotv), intent(in) :: Vector
      integer, dimension(ntotv), intent(in) :: Vector
     
     
      mrow = size(Matrix,1)
      ncol = size(Matrix,2)
      open(unit=10, file= fileplace//name1, ACTION="write", STATUS="replace")
     
      do i=1,mrow
        write(10, 100)( Matrix(i,j) ,j=1,ncol)
      end do
      close(10)
     
      open(unit=20, file= fileplace//name2, ACTION="write", STATUS="replace")
      do i=1,ncol
        !write(unit2, 100) Vector(i,1)
        write(20, 200) Vector(i)
      end do
      close(20)
      write(*,*) 'files: ', name1,' and ', name2, ' written succesfully on Res/'
      
      
      100 format (900E15.5)
      200 format (900I7)
      
    end subroutine writeMatrix 
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine PostPro_EMfield(N, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, glob_potential, grad_sol)
      
      implicit none
      
      character(len=*), parameter    :: fileplace = "Res/"
      double precision, dimension(nne,TotGp), intent(in) :: N, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp), intent(in) :: hes_xixi, hes_xieta, hes_etaeta
      double precision, dimension(ntotv), intent(in) :: glob_potential
      double precision, dimension(nevab)        :: elem_potential
      double precision, dimension(DimPr, dimPr) :: Jaco, Jinv
      double precision, dimension(DimPr,nne)    :: dN_dxy
      double precision, dimension(3,nne)        :: HesXY
      double precision, dimension(nne)          :: xi_cor, yi_cor, basis
      double precision, dimension(nne,DimPr)    :: element_nodes
      double precision, dimension(nelem)        :: du_dx, du_dy 
      integer         , dimension(nne)          :: nodeIDmap
      double precision                          :: detJ, x, y
      integer                                   :: ielem, igaus, inode! ,idime, ipoin, ii, jj
      double precision, allocatable, dimension(:,:,:), intent(out) :: grad_sol
      
      allocate(grad_sol(nelem,totGp,DimPr))

      if(postpro.eq.1)then
        write(*,'(A)') 'Running Post Process'
        goto 115 
      elseif(postpro.eq.2)then
        write(*,'(A)') 'None Post Process'
        goto 110
      else
        write(*,'(A)') 'No postrocess option defined'
      end if
     
      115 open(unit=100, file= fileplace//'electric_field'//'.dat', ACTION="write", STATUS="replace")
      
      write(100,'(2x,A5,8x,A5,10x,A2,15x,A2,16x,A,15x,A)') 'ielem','igaus','ex','ey','x', 'y'
      
      do ielem = 1,nelem
        call SetElementNodes(ielem, element_nodes, nodeIDmap, xi_cor, yi_cor)
        call gather(nodeIDmap, glob_potential, elem_potential) 
        
        !grad_sol = 0.0
       do igaus = 1, TotGp
          du_dx = 0.0 
          du_dy = 0.0 
          x = 0.0
          y = 0.0
          
          call Jacobian( element_nodes, dN_dxi, dN_deta, igaus ,Jaco, detJ, Jinv)
          call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, dN_dxy, HesXY)
          
          do inode = 1, nne
            basis(inode) = N(inode,igaus)
          end do
          
          !gradient computation (electric field)
          do inode = 1, nne
            du_dx(ielem) = du_dx(ielem) + dN_dxy(1,inode) * elem_potential(inode)
            du_dy(ielem) = du_dy(ielem) + dN_dxy(2,inode) * elem_potential(inode)
          end do
          
          grad_sol(ielem,igaus,1) = -du_dx(ielem)
          grad_sol(ielem,igaus,2) = -du_dy(ielem) 
          
          
          do inode = 1, nne 
            x = x + basis(inode)*xi_cor(inode)
            y = y + basis(inode)*yi_cor(inode)
            !print"(3(1x,f10.7))",  basis(ibase), yi_cor(ibase), basis(ibase)*yi_cor(ibase)
          end do
          
          write(100,'(2(I5,1x),1x,2(1x,F15.5),3x,2(E15.5))') ielem, igaus, x, y, grad_sol(ielem,igaus,1), grad_sol(ielem,igaus,2)
        end do
        !print*,' '
      end do
      
      close(100)
      110 continue 
      !do ielem =1,nelem
      !  do igaus = 1,TotGp
      !    write(*,'(2(1x,e15.5))') ( grad_sol(ielem, igaus, jj), jj=1,2 )
      !  end do
      !end do
      
    end subroutine PostPro_EMfield
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    !
    subroutine GID_results(solution, grad_sol)
      
      implicit none
      
      character(len=*), parameter    :: fileplace = "Pos/"
      double precision, dimension(ntotv, 1), intent(in) :: solution
      double precision, dimension(nelem,totGp,DimPr), intent(in) :: grad_sol
      character(len=10)                       :: extension1, extension2
      character(len=15)                       :: Elem_Type
      double precision, dimension(1, ntotv)   :: solution_T
      double precision, dimension(1,nnodes)   :: xcor, ycor
      integer                                 :: ipoin, ii, ielem, inode, jj, igaus, icomp
      
      solution_T = transpose(solution)
      xcor  = spread(coord(1,:),dim = 1, ncopies= 1)
      ycor  = spread(coord(2,:),dim = 1, ncopies= 1)
      
     
      print*, ' '
      print*, '!============== Output files ==================!'
      
      extension1 = ".post.msh"
      extension2 = ".post.res"
      if(ElemType.eq.'QUAD')then
        Elem_Type = 'Quadrilateral'
      elseif(ElemType.eq.'TRIA')then
        Elem_Type = 'Triangle'
      endif
      open(unit=555, file= fileplace//File_Nodal_Vals//extension1, ACTION="write", STATUS="replace")
      
      write(555,902) 'MESH', '"Domain"', 'dimension', DimPr, 'ElemType', Elem_Type, 'Nnode', nne
      write(555,"(A)") '#2D Convection-Diffusion-Reaction'
      write(555,900) '#Element tipe: ', ElemType,'/',ElemType
      
      
      write(555,"(A)")'Coordinates'
      write(555,"(A)") '#   No        X           Y'
      do ipoin = 1, nnodes
        write(555,906) ipoin, xcor(1,ipoin), ycor(1,ipoin)
      end do
      write(555,"(A)") 'End Coordinates'
      
      
      
      write(555,"(A)") 'Elements'
      do ielem=1,nelem
        write(555,908) ielem,(lnods(ielem,inode),inode=1,nne)
      end do
      write(555,"(A)") 'End Elements'
      
      
      close(555)
      write(*,"(A7,A21,A28)") ' -File ',File_Nodal_Vals//'.post.msh','written succesfully in Pos/'
      
      
      open(unit=555, file= fileplace//File_Nodal_Vals//extension2, ACTION="write", STATUS="replace")
      write(555,"(A)") 'GiD Post Results File 1.0'
      write(555,"(A)") '#2D Convection-Diffusion-Reaction'
      
      ! se escribe el res de las componentes de la velocidad
      select case(ndofn)
        case(1)
          write(555,"(A)") 'Result "phi" "Electric Potential" 0 Scalar OnNodes'
          write(555,"(A)") 'ComponentNames "" '
          write(555,"(A)") 'Values'
          write(555,*) '#',   'No    ','             ux '
          !  se escribe el res para el caso escalar de un grado de libertad
          write(555,914)
          do ipoin = 1, nnodes
            write(555,916) ipoin, solution(ipoin, 1)
          end do
          write(555,"(A)") 'End Values'
          
        case(2)
          write(555,"(A)") 'Result "EM Field" "EM field" 0 Vectsur OnNodes'
          write(555,"(A)") 'ComponentNames "ex" "ey" "--" "" '
          write(555,"(A)") 'Values'
          write(555,*) '#',   'No    ','             ex ','               ey '
          do ipoin = 1, nnodes
            write(555,918) ipoin, solution_T(1, ndofn*ipoin-1), solution_T(1,ndofn*ipoin)
          end do
          write(555,"(A)") 'End Values'
          
        case(3)
          write(555,"(A)") 'Result "EM field" "EM  field" 0 Vector OnNodes'
          write(555,"(A)") 'ComponentNames "ex" "ey" "--" "" '
          write(555,"(A)") 'Values'
          write(555,*) '#',   'No    ','             ex ','               ey'
          do ipoin = 1, nnodes
            write(555,919) ipoin, solution_T(1, ndofn*ipoin-2), solution_T(1,ndofn*ipoin-1)
          end do
          write(555,"(A)") 'End Values'
          write(555,"(A)") 'Result "Multiplier" "Multiplier" 0 Scalar OnNodes'
          write(555,"(A)") 'ComponentNames "" '
          write(555,"(A)") 'Values'
          write(555,*) '#',   'No    ','             p '
          !  se escribe el res para el caso escalar de un grado de libertad
          write(555,914)
          ii=1
          do ipoin = 3, nnodes*3,3
            write(555,916) ii, solution_T(1,ipoin)
            ii=ii+1
          end do
          write(555,"(A)") 'End Values'
          write(*,"(A7,A21,A28)") ' -File ',File_Nodal_Vals//'.post.res','written succesfully in Pos/'
      end select
      
      if(postpro.eq.1)then
        write(*,'(A)') 'Writing Post Process.....'
        goto 115 
      elseif(postpro.eq.2)then
        !write(*,'(A)') 'None Post Process'
        goto 110
      else
        write(*,'(A)') 'No postrocess option defined'
      end if

      115 write(555,"(A,A)") 'GaussPoints "GP_1" ElemType ', Elem_Type 
      write(555,"(A24,I1)") 'Number Of Gauss Points: ', totGp
      write(555,"(A,A)") 'Natural Coordinates: ', "Given"
      do igaus = 1, totGp
        write(555,901) (ngaus(igaus,jj), jj=1,DimPr)
      end do
      write(555,"(A)") 'End GaussPoints'
      
      write(555,"(A)") 'Result "Electric field" "ANALYSIS" 0 Vector OnGaussPoints "GP_1" '
      write(555,"(A)") 'ComponentNames "Ex" "Ey" '
      write(555,"(A)") 'Values'
      do ielem = 1, nelem
        write(555,'(I0)',Advance='NO') ielem
        do igaus = 1, totGp
          write(555,903) (grad_sol(ielem,igaus,icomp), icomp=1,2)
        end do
      end do
      write(555,"(A)") 'End Values'
      
      110 close(555)
      
      900 format(A15, A13, A1, A13)
      901 format(2(f10.5))
      902 format(A4,1x,A8,1X,A9,1X,I1,1X,A8,1X,A13,A6,1X,I1)
      903 format(2(E15.5))
      906 format(I7,4x,2(f15.5,3x)) !format for msh
      908 format(10(2x,I7) )
      914 format('#',3x,'No',     9x, 'Dof')
      916 format(I7,2x,E15.5)  !format for scalar case
      918 format(I7,3x,E12.5,3x,E12.5) !format for res velocity
      919 format(I7,2(4x,E15.5)) !format for res velocity
      
    end subroutine GID_results 
   
    
    
    
    
    !- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- - -!- -
    ! Agregamos las rutinas de tiempo del commit 625fbee: May2 2022
    
    
    subroutine GID_PostProcess(solution, activity, step_value, interval, time_final)
      
      implicit none
      
      character(len=*), parameter    :: fileplace = "Pos/"
      double precision, dimension(ntotv, 1), intent(in) :: solution
      character(*), intent(in)                :: activity
      integer, intent(in)                     :: step_value
      real, intent(in)                        :: interval, time_final
      character(len=10)                       :: extension1, extension2
      character(len=15)                       :: Elem_Type
      double precision, dimension(1, ntotv)   :: solution_T
      double precision, dimension(1,nnodes)   :: xcor, ycor
      integer                                 :: ipoin, ii, ielem, inode
      
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
      !open(unit=555, file= fileplace//File_Nodal_Vals//extension1, ACTION="write", STATUS="replace")
      
      
      
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
        print"(A11,A19,A30)", ' Mesh file ',File_Nodal_Vals//'.post.msh', 'written succesfully in Pos/ '
      ! if(status.eq.0)then continue 
      elseif(activity == "res")then
       
        if(step_value == 0)then
          open(unit=200, file= fileplace//File_Nodal_Vals//extension2, ACTION="write", STATUS="replace")
          write(200,"(A)") 'GiD Post Results File 1.0'
          write(200,"(A)") '#2D Convection-Diffusion-Reaction'
        else
          continue
        endif
        open(unit=200, file= fileplace//File_Nodal_Vals//extension2, ACTION="write", STATUS="old", position="append")
        
        ! se escribe el res de las componentes de la velocidad
        select case(ndofn)
          case(1)
            write(200,"(A29, I3, A)") 'Result "DoF" "Concentration" ', step_value,' Scalar OnNodes'
            write(200,"(A)") 'ComponentNames "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','             ux '
            !  se escribe el res para el caso escalar de un grado de libertad
            write(200,914)
            do ipoin = 1, nnodes
              write(200,916) ipoin, solution(ipoin, 1)
            end do
            !An alternative work around is to explicitly designate the elements to be read using an io-implied-do.
            !Something like
            !read (unit=10, fmt=*, iostat=iostat) (mat(pcnt,i),i=1,m)
            write(200,"(A)") 'End Values'
          case(2)
            write(200,"(A29, I3, A)") 'Result "DoF" "Concentration" ', step_value,' Vector OnNodes'
            write(200,"(A)") 'ComponentNames "u" "v" "--" "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','             ux ','               uy '
            do ipoin = 1, nnodes
              write(200,918) ipoin, solution_T(1, ndofn*ipoin-1), solution_T(1,ndofn*ipoin)
            end do
            write(200,"(A)") 'End Values'
          case(3)
            write(200,"(A29, I3, A)") 'Result "DoF" "Concentration" ', step_value,' Vector OnNodes'
            write(200,"(A)") 'ComponentNames "u" "v" "w" "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','             ux ','               uy'
           ! do ipoin = 1, nnodes
           !   write(200,919) ipoin, solution_T(1, ndofn*ipoin-2), solution_T(1,ndofn*ipoin-1), solution_T(1,ndofn*ipoin)
           ! end do
            do ipoin = 1, nnodes
              write(200,919) ipoin, solution_T(1, ndofn*ipoin-2), solution_T(1,ndofn*ipoin-1)
            end do
            write(200,"(A)") 'End Values'
            write(200,"(A22, I3, A)") 'Result "P" "Preassure"', step_value,' Scalar OnNodes'
            write(200,"(A)") 'ComponentNames "" '
            write(200,"(A)") 'Values'
            write(200,*) '#',   'No    ','     P '
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
      close(200)
      !la siguiente instruccion debe usarse con nt no con time pero solo es para avanzar
      if(interval == time_final) then
        print*, ' '
        print"(1x, A19,A30)", File_Nodal_Vals//'.post.res', 'written succesfully in Pos/ '
      endif
      
      !if(simul.eq.6.and.ndofn.1)then
      !  do ipoin = 1, nnodes
      !    write(555,906) ipoin, xcor(1,ipoin), ycor(1,ipoin)
      !  end do
      !else
      !  write(*,'A')
      !
      
      
      
      
      
      
      
      
      
      900 format(A15, A13, A1, A13)
      902 format(A4,1x,A8,1X,A9,1X,I1,1X,A8,1X,A13,A6,1X,I1)
      906 format(I7,2(3x,f9.4)) !format for msh
      908 format(9(2x,I7) )
      914 format('#',3x,'No',     9x, 'Dof')
      916 format(I7,2x,E12.5)  !format for scalar case
      918 format(I7,3x,E15.5,3x,E15.5) !format for res velocity
      919 format(I7,3(3x,E15.5)) !format for res velocity
      
    end subroutine GID_PostProcess
   



    
    
    
    
    
    
    subroutine GlobalSystem_Time(N,dN_dxi,dN_deta,hes_xixi,hes_xieta,hes_etaeta,S_ldSol,delta_t, ugl_pre , A_F)
      
      use sourceTerm
      
      implicit none
      
      double precision, allocatable, dimension(:,:), intent(in out) :: ugl_pre
      double precision, dimension(nne,TotGp), intent(in) :: N, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp), intent(in) :: hes_xixi, hes_xieta, hes_etaeta
      !double precision, dimension(3,nne), intent(in)     :: Hesxieta
      double precision, dimension(ndofn)        :: EMsource
      integer, intent(in)                                :: S_ldSol
      double precision, dimension(nne)          :: basis, xi_cor, yi_cor
      double precision, dimension(DimPr,nne)    :: dN_dxy
      double precision, dimension(3,nne)        :: HesXY
      double precision, dimension(DimPr, dimPr) :: Jaco, Jinv
      double precision, dimension(nevab, nevab) :: Ke, Ce, rhs_CN
      double precision, dimension(nevab)        :: Fe, Fe_time, ue_pre, time_cont
      double precision, dimension(3,3)          :: tauma
      double precision, dimension(nne,DimPr)                :: element_nodes
      integer, dimension(nne)                   :: nodeIDmap
      double precision                          :: dvol, hmaxi, detJ, delta_t
      integer                                   :: igaus, ibase, ielem
      double precision, allocatable, dimension(:,:), intent(out)  :: A_F
      
      allocate( A_F(ntotv, 1) )
      !allocate(ugl_pre(S_ldSol,1) )
      
                                                                                                
      A_F = 0.0
      do ielem = 1, nelem 
        !gather
        Ke = 0.0
        Fe = 0.0    !Fe(nevab)
        Ce = 0.0
        call SetElementNodes(ielem, element_nodes, nodeIDmap, xi_cor, yi_cor)
        !call SetElementNodes(ielem, element_nodes, nodeIDmap)
        call gather(nodeIDmap, ugl_pre, ue_pre)
        time_cont = ue_pre * 1.0/delta_t
        
        !do-loop: compute element capacity and stiffness matrix Ke Ce and element vector Fe
        do igaus = 1, TotGp
          call Jacobian( element_nodes, dN_dxi, dN_deta, igaus ,Jaco, detJ, Jinv)
          !Jaco = J2D(element_nodes, dN_dxi, dN_deta, igaus)
          !detJ = m22det(Jaco)
          !Jinv = inv2x2(Jaco)
          dvol = detJ *  weigp(igaus,1)
          !call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, Hesxieta, dN_dxy, HesXY)
          call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, dN_dxy, HesXY)
          hmaxi = elemSize(Jinv)
          do ibase = 1, nne
            basis(ibase) = N(ibase,igaus)
          end do
          
          call TauMat(hmaxi,tauma)
          !call source_term(ielem, basis, xi_cor, yi_cor, EMsource)
          call Galerkin(hmaxi, dvol, basis, dN_dxy, EMsource, Ke, Ce, Fe) !amate lo llame Ke
          !call Galerkin(dvol, basis, dN_dxy, Ke, Ce, Fe) 
          !!call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, Ke, Fe, pertu,workm,resid)
          if(kstab.ne.6.or.kstab.ne.0)call Stabilization(dvol, basis, dN_dxy, HesXY, EMsource, tauma, Ke, Fe)
          
          
          select case(theta)
          case(2)
          Fe_time = Fe + matmul(Ce,time_cont)
          case(3)
          rhs_CN  = (1.0/delta_t)*Ce - 0.5*Ke
          Fe_time = 0.5*Fe + matmul(rhs_CN,ue_pre)
          endselect
        
                                                                                                
                                                                                                
        end do
        
        
        !call Assemb_Glob_Mat(nodeIDmap, Ke, A_K)      !Assemble Global Conductivity Matrix K
        !call Assemb_Glob_Mat(nodeIDmap, Ce, A_C)      !Assemble Global Capacity Matrix C 
        call Assemb_Glob_Vec(nodeIDmap, Fe_time, A_F) !Assemble Global Source vector F
        
      end do
      
      
    end subroutine GlobalSystem_Time




    subroutine Galerkin_prevTime(dvol, basis, Ce, Fe)
      
      implicit none
      
      double precision, intent(in) :: basis(nne)
      double precision, intent(in) :: dvol
      integer :: inode, idofn, ievab, jevab, jnode, jdofn
      double precision ::  cpcty
      double precision, intent(out) :: Fe(nevab), Ce(nevab,nevab)
      ievab=0
      do inode=1,nne
        do idofn=1,ndofn
          ievab=ievab+1
          jevab=0
          do jnode=1,nne
            do jdofn=1,ndofn
              jevab=jevab+1
              cpcty = basis(inode) * basis(jnode)
              Ce(ievab,jevab) = Ce(ievab,jevab) + cpcty * dvol                                 !element Capacity (Mass) matrix
            end do
          end do
          Fe(ievab) = Fe(ievab) + basis(inode) * force(idofn) * dvol
        end do
      end do
    
    end subroutine Galerkin_prevTime
   
    
    
    
    subroutine TimeLevels(N, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, delta_t, Ucurr, Uprev, F_plus_MU)
      
      use sourceTerm
      
      
      implicit none
      
      double precision, allocatable, dimension(:,:), intent(in out) :: Uprev, Ucurr
      double precision, dimension(nne,TotGp), intent(in):: N, dN_dxi, dN_deta
      double precision, dimension(nne,TotGp), intent(in) :: hes_xixi, hes_xieta, hes_etaeta
      !double precision, dimension(3,nne), intent(in)    :: Hesxieta
      double precision, dimension(ndofn)        :: EMsource
      double precision, dimension(nne)          :: basis, xi_cor, yi_cor
      double precision, dimension(DimPr,nne)    :: dN_dxy
      double precision, dimension(3,nne)        :: HesXY
      double precision, dimension(DimPr, dimPr) :: Jaco, Jinv
      double precision, dimension(nevab, nevab) :: Ke, Ce
      double precision, dimension(nevab)        :: Fe, Fe_time, ue_prev, ue_curr, AvrgeTime
      double precision, dimension(3,3)          :: tauma
      double precision, dimension(nne,DimPr)                :: element_nodes
      integer, dimension(nne)                   :: nodeIDmap
      double precision                          :: dvol, hmaxi, detJ, delta_t
      integer                                   :: igaus, ibase, ielem
      double precision, dimension(ntotv, 1), intent(out)  :: F_plus_MU
      
      !allocate( F_plus_MU )
      F_plus_MU = 0.0
    
      do ielem = 1, nelem 
        Ke = 0.0; Fe = 0.0; Ce = 0.0
        call SetElementNodes(ielem, element_nodes, nodeIDmap, xi_cor, yi_cor)
        call gather(nodeIDmap, Uprev, ue_prev)
        call gather(nodeIDmap, Ucurr, ue_curr)
        AvrgeTime =  (2*ue_curr + 0.5*ue_prev) / delta_t 
        
        !do-loop: compute element capacity and stiffness matrix Ke Ce and element vector Fe
        do igaus = 1, TotGp
          call Jacobian( element_nodes, dN_dxi, dN_deta, igaus ,Jaco, detJ, Jinv)
          !Jaco = J2D(element_nodes, dN_dxi, dN_deta, igaus)
          !detJ = m22det(Jaco)
          !Jinv = inv2x2(Jaco)
          !dvol = detJ *  weigp(igaus,1)
          !call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, Hesxieta, dN_dxy, HesXY)
          call DerivativesXY(igaus, Jinv, dN_dxi, dN_deta, hes_xixi, hes_xieta, hes_etaeta, dN_dxy, HesXY)
          hmaxi = elemSize(Jinv)
          do ibase = 1, nne
            basis(ibase) = N(ibase,igaus)
          end do
          
          
          !call Galerkin(dvol, basis, dN_dxy, Ke, Ce, Fe) 
          call source_term(ielem, basis, xi_cor, yi_cor, EMsource)
          call Galerkin_prevTime(dvol, basis, Ce, Fe)
          call TauMat(hmaxi,tauma)
          !!call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, Ke, Fe, pertu,workm,resid)
          if(kstab.ne.6.or.kstab.ne.0)call Stabilization(dvol, basis, dN_dxy, HesXY, EMsource, tauma, Ke, Fe)
          !call Stabilization(dvol, basis, dN_dxy, HesXY, tauma, Ke, Fe)
          Fe_time = Fe - matmul(Ce,AvrgeTime)
        end do
    
        call Assemb_Glob_Vec(nodeIDmap, Fe_time, F_plus_MU) !Assemble Global Source vector F
        
      end do
      
    end subroutine TimeLevels



    
    
    
    
    
    
    
    
  !Fin de contains


end module library
