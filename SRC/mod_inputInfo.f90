module inputInfo
  use param
  use geometry

  contains
    !
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
    ! 
    subroutine GeneralInfo(name_inputFile, geometry_File)
      external :: fdate
      
      character(len=*), parameter :: fileplace = "Info/"
      character(len=19)           :: name_inputFile
      character(len=5)            :: file_name 
      character(len=24)           :: date
      character(len=4)            :: aaaa, cccc
      character(len=9)            :: Prob_Type
      character(len=12)           :: bbbb, geometry_File
      character(len=16)           :: dddd, OrderElemType
      !double precision            :: delta_t
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
     
      if(ProbType.ne.'TIME')then
        Prob_Type = 'STATIC'
      else
        Prob_Type = 'TRANSIENT'
      endif
      
      if((ElemType.eq.'QUAD').and.(nne.eq.4))then
        OrderElemType = 'QUADRILATERAL Q1'
      elseif((ElemType.eq.'QUAD').and.(nne.eq.9))then
        OrderElemType = 'QUADRILATERAL Q2'
      elseif((ElemType.eq.'TRIA').and.(nne.eq.3))then
        OrderElemType = 'TRIANGULAR P1'
      elseif((ElemType.eq.'TRIA').and.(nne.eq.6))then
        OrderElemType = 'TRIANGULAR P2'
      end if
      
      
     
      call fdate(date)
      print*, ' '
      print*, '- - - - 2D Convection-Diffusion-Reaction Simulation - - - - '
      print*, ' '
      print*,' ',date
      print*,'!================= GENERAL INFO ===============!'
      write(*,"(A30,2x,a19  ,3X,A1 )") ' - Input File               : ', name_inputFile,''
      write(*,"(A30,2x,a13  ,3X,A1 )") ' - Mesh File                : ', geometry_File,''
      write(*,"(A30,2x,a16  ,3X,A1 )") ' - Element type             : ', OrderElemType,''
      write(*,"(A30,2x,a9   ,3X,A1 )") ' - Problem Type             : ', Prob_Type,''
      write(*,"(A30,2X,I6   ,1X,A10)") ' - Problem dimension        : ', DimPr, '  '
      write(*,"(A30,2X,I6   ,1X,A10)") ' - Elements                 : ', initelem,'   '
      write(*,"(A30,2X,I6   ,1X,A10)") ' - Nodal points             : ', initnodes, ' '
      write(*,"(A30,2X,I6   ,1X,A10)") ' - DoF per node             : ', ndofn, '  '
      write(*,"(A30,2X,I6   ,1X,A10)") ' - Nodes per element        : ', initnne, '    '
      write(*,"(A30,2X,I6   ,1X,A10)") ' - Total Gauss points       : ', totGp,'   '
      write(*,"(A30,2X,I6   ,1X,A10)") ' - Element variabless       : ', initnevab    ,'  '
      write(*,"(A30,2X,I6   ,1X,A10)") ' - Total unknowns           : ', initntotv    ,'  '
      write(*,"(A29,3X,f11.4,1X,A10)") ' - Characteristic mesh size : ', helem,' '
      write(*,"(A30,2X,f9.2 ,1X,A10)") ' - Length reference element : ', hnatu        ,' '
      write(*,"(A30,6X,f9.5 ,1X,A10)") ' - Model condutivity (σ)    : ', sigma        ,' ' 
      
      if(refiType.eq.'NO')then
        write(*,"(A30,2x,A7,3X,A10)") ' - Refinement type          : ','  NONE',' '
      elseif(refiType.eq.'PS'.or.refitype.eq.'CC')then
        if(refitype.eq.'PS')then
          bbbb = 'Powell-Sabin'
        elseif(refiType.eq.'CC')then
          bbbb = 'Criss-Cross'
        end if
          print*, ' '
          print*,'!================= REFINMENT INFO ===============!'
          write(*,"(A30,2x,a12,3X,A1)") ' - Refinement type           : ', bbbb,''
          write(*,"(A30,2X,I6,1X,A10)") ' - Elements after refinement : ', nelem,'   '
          write(*,"(A30,2X,I6,1X,A10)") ' - Nodes after refinement    : ', nnodes, ' '
          write(*,"(A30,2X,I6,1X,A10)") ' - Elm. unk. after refinement: ', nevab, ' '
          write(*,"(A30,2X,I6,1X,A10)") ' - Unknowns after refinement : ', ntotv, ' '
      else
        write(*,'(A)') '> > >Error in refinment type'
      endif
       
        
      print*,' '
      print*,'!============ FOURIER PARAMNETERS =============!'
      if(TwoHalf == 'Y')then
        write(*,"(A30,2X,A8   ,1X,A10)") ' - A 2.5D modeling?         : ', 'YES', '    '
        write(*,"(A30,6X,e11.4 ,1X,A10)") ' - Min wavenumber           : ', ky_min,'   '
        write(*,"(A30,6X,e11.4 ,1X,A10)") ' - Max wavenumber           : ', ky_max    ,'  '
        write(*,"(A30,5X,I3   ,1X,A10)") ' - Total wavenumbers        : ', tot_ky   ,'  '
        write(*,"(A30,2X,f9.2 ,1X,A10)") ' - Location of receiver     : ', y_iFT        ,' '
      else
        write(*,"(A30,2X,A7   ,1X,A10)") ' - A 2.5D modeling?         : ', 'NONE', '    '
      endif
      
      
      print*, ' '
      if(kstab.eq.0)then
        print*,'!========== STABILIZATION PARAMETERS ==========!'
        write(*,"(A26,5x,A8,3X,A1)") ' - Stabilization method:   ', aaaa,''
      elseif(kstab.eq.6)then
        print*,'!========== STABILIZATION PARAMETERS ==========!'
        write(*,"(A30,2x, A10 ,1X,A10)") ' - Stabilization method     : ', aaaa   ,' '
        write(*,"(A30,2X,f13.5,1X,A10)") ' - Reluctivity of medium (λ): ', lambda ,' '
        write(*,"(A30,1X,f13.5,1X,A10)") ' - Algorithmic constant (Cu): ', Cu     ,' '
        write(*,"(A31,2X,f13.5,1X,A10)") ' - Constante of length (ℓ)  : ', ell    ,' '
        write(*,"(A30,2X,e16.5,1X,A10)") ' - Stab. param.1 (Su)       : ', Cu*lambda*(helem**2/ell**2),'   '
        write(*,"(A30,2X,e16.5,2X,A10)") ' - Stab. param.2 (Sp)       : ', ell**2 / lambda,'   '
      else
        print*,'!========== STABILIZATION PARAMETERS ==========!'
        write(*,"(A30,2x,a4  ,3X,A1 )")  ' - Stabilization method     : ', aaaa   ,' '
        write(*,"(A30,2x,I2  ,3X,A1 )")  ' - Type of Tau matrix       : ', ktaum  ,' '
        write(*,"(A30,2X,f3.1,1X,A10)")  ' - Param. to obtain TAU     : ', patau  ,' '
      endif
      
      print*, ' '
      if(ProbType.eq.'TIME')then
        !!delta_t ( time_fin - time_ini ) / (t_steps + 1.0)   !Step size
        print*,'!============ TIME DISCRETIZATION =============!'
        if((theta.eq.2).or.(theta.eq.4))then
          if(theta.eq.2)then
            cccc = 'BDF1'
          else
            cccc = 'BDF2'
          endif
          write(*,"(A29,2X,A8,1X,A10)") ' - Method Selected          : ', cccc,' '
        elseif(theta.eq.3)then
          dddd = 'Cranck-Nicholson'
          write(*,"(A29,3x,A16,1X,A10)") ' - Method Selected          : ', dddd,' '
        endif
        write(*,"(A29,6X,E13.5,1X,A11)") ' - Begining time            : ', time_ini,' '
        write(*,"(A29,6X,E13.5,1X,A11)") ' - End time                 : ', time_fin,' '
        write(*,"(A29,2X,I7   ,1X,A11)") ' - Number of steps          : ', t_steps,' '
        write(*,"(A31,6X,E13.5,1X,A11)") ' - Step size (∆t)           : ', delta_t ,' '
      else
        continue
      endif
      print*, ' '
      print*,'!============ Source Parameters =============!'
      if(srcRHS == 0) then 
        write(*,'(A)') ' -Location'
        if(nodalSrc.eq.1)then
          write(*,'(A)')' -Source point' 
          write(*,'(A,F10.3,A,F10.3,A)') '(',coord(1,Srcloc(1)),',',coord(2,Srcloc(1)),' ) '
          print*, ' '
          write(*,'(A24,99(I0,4x))')' -Source point at node: ', (Srcloc(i), i=1,nodalSrc) 
          print*, ' '
          write(*,'(A)') ' -Dipole lenght: Single source point'
        elseif(nodalSrc.eq.2)then
          write(*,*)'              Begining                End' 
          write(*,'(A,F8.3,A,F8.3,A,A,F8.3,A,F8.3,A)') &
            &'        (',coord(1,Srcloc(1)),',',coord(2,Srcloc(1)),' ) ',&
            &       ' (',coord(1,Srcloc(2)),',',coord(2,Srcloc(2)),' )'
          print*, ' '
          write(*,'(A25,99(I0,4x))')' -Nodes involves source: ', (Srcloc(i), i=1,nodalSrc) 
          print*, ' '
          write(*,'(A17,f5.2)') ' -Dipole lenght: ', abs(coord(1,Srcloc(2)) - coord(1,Srcloc(1))) 
        else
          write(*,*)'              Begining                End' 
          write(*,'(A,F8.3,A,F8.3,A,A,F8.3,A,F8.3,A)') &
            &'        (',coord(1,Srcloc(1)),',',coord(2,Srcloc(1)),' ) ',&
            &       ' (',coord(1,Srcloc(nodalSrc)),',',coord(2,Srcloc(nodalSrc)),' )'
          print*, ' '
          write(*,'(A25,99(I0,4x))')' -Nodes involves source: ', (Srcloc(i), i=1,nodalSrc) 
          !write(*,'(I6)') nodalSrc 
          print*, ' '
          write(*,'(A17,f5.2)') ' -Dipole lenght: ', abs(coord(1,Srcloc(nodalSrc)) - coord(1,Srcloc(1))) 
          !write(*,'(f8.2)')  abs(coord(1,Srcloc(2)) - coord(1,Srcloc(1))) 
        end if
        print*, ' '
        write(*,'(A)') ' -Intensity current'
        if(ndofn.eq.1)then
          write(*,"(1(f10.3,1x))") Icurr(1)
        elseif(ndofn.eq.3)then
          write(*,"(3(f10.3,1x))") Icurr(1), Icurr(2), Icurr(3)
        endif
      else
        write(*,'(A)') ' - Not geophysical source'
      end if
      
      write(*,'(A)') 
      print*,'!============ TENSOR COEFFICIENTS  ============!'
      print*, 'Diffusion'
      do i = 1,dimPr
        do j = 1,DimPr
          print"(A,2I1)", 'k_',i,j
          do k = 1,ndofn
            print"(e15.7,1x,e15.7, 1x, e15.7)",( difma(k,l,i,j), l=1,ndofn)
          end do
          !print*,' '
        end do
      end do
      print*, ' '  
      print*, 'Convection'
      do k = 1, DimPr
        print"(A,2I1)",'A_',k
        do i = 1, ndofn
          write(*, "(f10.5, 1x, f10.5, 1x, f15.5)")( conma(i,j,k) ,j=1, ndofn)
        end do
      end do
        print*,' '
      print*,'Reaction'
      do i=1,ndofn
        write(*,"(f10.5, 1x, f10.5, 1x, f15.5)" )( reama(i,j) ,j=1,ndofn)
      end do
      print*, ' '
      print*, 'External Forces'
      if(ndofn.eq.1)then
        write(*,"(3(f10.5,1x))") force(1)
      elseif(ndofn.eq.2)then
        write(*,"(2(f10.5,1x))") force(1), force(2)
      else
        write(*,"(3(f10.5,1x))") force(1), force(2), force(3)
      endif
      write(*,'(A)') 
      
      file_name ="test_"
      open(unit=100,file= fileplace//file_name//testID//'.txt', ACTION="write", STATUS="replace")
      
      if(refiType.eq.'NO')then
        bbbb = '    NONE'
      elseif(refiType.eq.'PS')then
        bbbb = 'Powell-Sabin'
      elseif(refiType.eq.'CC')then
        bbbb = 'Criss-Cross'
      else
        write(*,'(A)') '> > >Error in refinment type'
      endif
      
      write(100,'(A)')'- - - - 2D Convection-Diffusion-Reaction Simulation - - - - '
      write(100,'(A)')
      write(100,'(A8,1x,A14)') 'test ID: ',testID
      write(100,'(A)') " "
      write(100,'(A)') date
      write(100,'(A)')'!================= GENERAL INFO ===============!'
      write(100,"(A30,2x,a19  ,3X,A1 )") ' - Input File               : ', name_inputFile,''
      write(100,"(A30,2x,a12  ,3X,A1 )") ' - Mesh File                : ', geometry_File,''
      write(100,"(A30,2x,a16  ,3X,A1 )") ' - Element type             : ', OrderElemType,''
      write(100,"(A30,2x,a9   ,3X,A1 )") ' - Problem Type             : ', Prob_Type,''
      write(100,"(A30,2X,I6   ,1X,A10)") ' - Elements                 : ', initelem,'   '
      write(100,"(A30,2X,I6   ,1X,A10)") ' - Nodal points             : ', initnodes, ' '
      write(100,"(A30,2X,I6   ,1X,A10)") ' - DoF per node             : ', ndofn, '  '
      write(100,"(A30,2X,I6   ,1X,A10)") ' - Nodes per element        : ', initnne, '    '
      write(100,"(A30,2X,I6   ,1X,A10)") ' - Total Gauss points       : ', totGp,'   '
      write(100,"(A30,2X,I6   ,1X,A10)") ' - Total unknowns           : ', initntotv    ,'  '
      write(100,"(A29,3X,f11.4,1X,A10)") ' - Element size             : ', helem,' '
      write(100,"(A30,6X,f9.5 ,1X,A10)") ' - Model condutivity (σ)    : ', sigma        ,' ' 
      if(refiType.eq.'NO')then
        write(100,"(A30,2x,a7 ,3X,A10)") ' - Refinement type          : ', '  NONE',' '
        write(100,'(A)')
      elseif(refiType.ne.'NO')then
        write(100,'(A)') 
        write(100,'(A)')'!================= REFINMENT INFO ===============!'
        write(100,"(A30,2x,a12,3X,A1 )")' - Refinement type           : ', bbbb    ,' '
        write(100,"(A30,2X,I6,1X,A10)") ' - Elements after refinement : ', nelem,'   '
        write(100,"(A30,2X,I6,1X,A10)") ' - Nodes after refinement    : ', nnodes, ' '
        write(100,"(A30,2X,I6,1X,A10)") ' - Elm unkns after refinement: ', nevab, ' '
        write(100,"(A30,2X,I6,1X,A10)") ' - Glb unkns after refinement: ', ntotv, ' '
        write(100,'(A)')
      endif 
      write(100,'(A)')'!============ FOURIER PARAMNETERS =============!'
      if(TwoHalf == 'Y')then
        write(100,"(A30,2X,A8   ,1X,A10)") ' - A 2.5D modeling?         : ', 'YES', '    '
        write(100,"(A30,6X,e11.4 ,1X,A10)") ' - Min wavenumber           : ', ky_min,'   '
        write(100,"(A30,6X,e11.4 ,1X,A10)") ' - Max wavenumber           : ', ky_max    ,'  '
        write(100,"(A30,5X,I3   ,1X,A10)") ' - Total wavenumbers        : ', tot_ky   ,'  '
        write(100,"(A30,2X,f9.2 ,1X,A10)") ' - Location of receiver     : ', y_iFT        ,' '
      else
        write(100,"(A30,2X,A7   ,1X,A10)") ' - A 2.5D modeling?         : ', 'NONE', '    '
      endif
      if(kstab.eq.0)then
        write(100,'(A)') 
        write(100,'(A)')'!========== STABILIZATION PARAMETERS ==========!'
        write(100,"(A26,6x,a7,3X,A1)") ' - Stabilization method:       ', aaaa,''  
      !write(100,"(A26,3X,f3.1,1X,A10)") ' - Exponent of mesh size:    ', i_exp,'   '
      elseif(kstab.eq.6)then
        write(100,'(A)') 
        write(100,'(A)')'!========== STABILIZATION PARAMETERS ==========!'
        
        write(100,"(A30,2x, A10 ,1X,A10)") ' - Stabilization method     : ', aaaa   ,' '
        write(100,"(A30,2X,f13.5,1X,A10)") ' - Reluctivity of medium (λ): ', lambda ,' '
        write(100,"(A30,1X,f13.5,1X,A10)") ' - Algorithmic constant (Cu): ', Cu     ,' '
        write(100,"(A31,2X,f13.5,1X,A10)") ' - Constante of length (ℓ)  : ', ell    ,' '
        write(100,"(A30,2X,e16.5,1X,A10)") ' - Stab. param.1 (Su)       : ', Cu*lambda*(helem**2/ell**2),'   '
        write(100,"(A30,2X,e16.5,2X,A10)") ' - Stab. param.2 (Sp)       : ', ell**2 / lambda,'   '
      else
        write(100,'(A)') 
        write(100,'(A)')'!========== STABILIZATION PARAMETERS ==========!'
        write(100,"(A30,2X,A10  ,2X,A10)") ' - Stabilization method     : ', aaaa   ,' '
        write(100,"(A30,2X,I2   ,2X,A10)") ' - Type of Tau matrix       : ', ktaum  ,' '
        write(100,"(A30,2X,f10.3,2X,A10)") ' - Param. to obtain TAU     : ', patau  ,' '
      endif
      if(ProbType.eq.'TIME')then
        !delta_t 1e-3! ( time_fin - time_ini ) / (t_steps + 1.0)   !Step size
        write(100,'(A)') 
        write(100,'(A)')'!============ TIME DISCRETIZATION =============!'
        if((theta.eq.2).or.(theta.eq.4))then
          if(theta.eq.2)then
            cccc = 'BDF1'
          else
            cccc = 'BDF2'
          endif
          write(100,"(A29,6X,A8  ,1X,A10)") ' - Method Selected          : ', cccc,' '
        elseif(theta.eq.3)then
          dddd = 'Cranck-Nicholson'
          write(100,"(A29,3x,a16,3X,A11)") ' - Method Selected          : ', dddd,' '
        endif
        write(100,"(A29,6X,E13.5,1X,A11)") ' - Begining time            : ', time_ini,' '
        write(100,"(A29,6X,E13.5,1X,A11)") ' - Time simulated           : ', time_fin,' '
        write(100,"(A29,2X,I7   ,1X,A11)") ' - Number of steps          : ', t_steps,' '
        write(100,"(A31,6X,E13.5,1X,A11)") ' - Step size (∆t)           : ', delta_t ,' '
      else
        continue
      endif
      write(100,'(A)') 
      write(100,'(A)')'!============ Source Parameters =============!'
      if(srcRHS == 0)then
        write(100,'(A)') ' -Location'
        if(nodalSrc.eq.1)then
          write(100,'(A)')' -Source point' 
          write(100,'(A,F8.3,A,F8.3,A)') '(',coord(1,Srcloc(1)),',',coord(2,Srcloc(1)),') '
          write(100,'(A)') 
          write(100,'(A25,2(I0,4x))')' -Nodes involves source: ', (Srcloc(i), i=1,nodalSrc) 
          write(100,'(A)') 
          write(100,'(A)') ' -Dipole lenght: Single source point'
        elseif(nodalSrc.eq.2)then
          write(100,'(A)')'              Begining                End' 
          write(100,'(A,F8.3,A,F8.3,A,A,F8.3,A,F8.3,A)') &
            &'        (',coord(1,Srcloc(1)),',',coord(2,Srcloc(1)),' ) ',' (',coord(1,Srcloc(2)),',',coord(2,Srcloc(2)),' )'
          write(100,'(A)') 
          write(100,'(A25,99(I0,4x))')' -Nodes involves source: ', (Srcloc(i), i=1,nodalSrc) 
          write(100,'(A)') 
          write(100,'(A17,f5.2)') ' -Dipole lenght: ', abs(coord(1,Srcloc(2)) - coord(1,Srcloc(1))) 
        else
          write(100,'(A)')'              Begining                End' 
          write(100,'(A,F8.3,A,F8.3,A,A,F8.3,A,F8.3,A)') &
            &'        (',coord(1,Srcloc(1)),',',coord(2,Srcloc(1)),' ) ',&
            &       ' (',coord(1,Srcloc(nodalSrc)),',',coord(2,Srcloc(nodalSrc)),' )'
          write(100,'(A)') 
          write(100,'(A25,99(I0,4x))')' -Nodes involves source: ', (Srcloc(i), i=1,nodalSrc) 
          !write(100,'(I6)') nodalSrc 
          write(100,'(A)') 
          write(100,'(A17,f5.2)') ' -Dipole lenght: ', abs(coord(1,Srcloc(nodalSrc)) - coord(1,Srcloc(1))) 
          !write(100,'(f8.2)')  abs(coord(1,Srcloc(2)) - coord(1,Srcloc(1))) 
        end if
      else
        write(100,'(A)') ' - Not geophysical source'
      endif
      
      
      
      
      write(100,'(A)') 
      write(100,'(A)') ' -Intensity current'
      if(ndofn.eq.1)then
        write(100,"(1(f10.3,1x))") Icurr(1)
      elseif(ndofn.eq.3)then
        write(100,"(3(f10.3,1x))") Icurr(1), Icurr(2), Icurr(3)
      endif
      write(100,'(A)') 
      write(100,'(A)') 
      write(100,'(A)')'!============ TENSOR COEFFICIENTS  ============!'
      write(100,'(A)') 'Diffusion'
      write(100,'(A)') 
      do i = 1,dimPr
        do j = 1,DimPr
          write(100,"(A,2I1)") 'k_',i,j
          do k = 1,ndofn
            write(100,"(e15.7,1x,e15.7, 1x, e15.7)") ( difma(k,l,i,j), l=1,ndofn)
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
end module inputInfo
