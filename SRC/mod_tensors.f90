module tensor_inputs
  use param
  implicit none
  
  double precision, allocatable, dimension(:,:,:,:) :: difma
  double precision, allocatable, dimension(:,:,:)   :: conma
  double precision, allocatable, dimension(:,:)     :: reama 
  double precision, allocatable, dimension(:)       :: force
  contains
    subroutine PDE_Coefficients
      
      implicit none
      !En primera instancia, esta rutina no tendria salida
      !Pues las variables de salida, los tensores mismos, son gliobales y deben
      !Estar disponibles para todo el codigo
      double precision, allocatable, dimension(:)   :: difma_elements, conma_elements, reama_elements, force_elements
      
      allocate( difma(ndofn,ndofn,DimPr,DimPr) )
      allocate( conma(ndofn,ndofn,DimPr) )
      allocate( reama(ndofn,ndofn) )
      allocate( force(ndofn) )
      
      allocate( difma_elements(ndofn*ndofn*DimPr*DimPr) )
      allocate( conma_elements(ndofn*ndofn*DimPr) )
      allocate( reama_elements(ndofn*ndofn) )
      allocate( force_elements(ndofn) )
      
      
      !print*, 'Entra a PDE'
      !print*, ' '
      !print*, PHYSICAL_PROBLEM 
      !print*, ''
      !
      difma = 0.0 ; conma = 0.0
      reama = 0.0 ; force = 0.0
      
      select case(PHYSICAL_PROBLEM)
        case('Maxwell_In_Non_Convex_Domain')
          difma_elements = (/ &
          ! Capa 1 (dimension 1:1)
          1.0 , 0.0, 0.0, &
          0.0 , 1.0, 0.0, &
          0.0 , 0.0, 1.0, &
          ! Capa 2 (dimension 1:2)
          0.0 , 1.0, 0.0, &
          -1.0, 0.0, 0.0, &
          0.0 , 0.0, 0.0, &
          ! Capa 3 (dimension 2:1)
          0.0 ,-1.0, 0.0, &
          1.0 , 0.0, 0.0, &
          0.0 , 0.0, 0.0, &
          ! Capa 4 (dimension 2:2)
          1.0 , 0.0, 0.0, &
          0.0 , 1.0, 0.0, &
          0.0 , 0.0, 1.0 /)
          
          conma_elements = (/ &
          ! Capa 1 (dimension 1)
          0.0, 0.0, 1.0, &
          0.0, 0.0, 0.0, &
          1.0, 0.0, 0.0, &
          ! Capa 2 (dimension 2)
          0.0, 0.0, 0.0, &
          0.0, 0.0, 1.0, &
          0.0, 1.0, 0.0 /)
         
          reama_elements = (/ &
          ! Capa 1 (dimension 1:1)
          0.0, 0.0, 0.0, &
          0.0, 0.0, 0.0, &
          0.0, 0.0, 0.0 /)
          
          force_elements = (/ 1.0, 1.0, 0.0 /)
         
          
        case('Cavity_Driven_Flow')
          difma_elements = (/ &
          ! Capa 1 (dimension 1:1)
          1.0, 0.0, 0.0, &
          0.0, 1.0, 0.0, &
          0.0, 0.0, 0.0, &
          ! Capa 2 (dimension 1:2)
          0.0, 0.0, 0.0, &
          0.0, 0.0, 0.0, &
          0.0, 0.0, 0.0, &
          ! Capa 3 (dimension 2:1)
          0.0, 0.0, 0.0, &
          0.0, 0.0, 0.0, &
          0.0, 0.0, 0.0, &
          ! Capa 4 (dimension 2:2)
          1.0, 0.0, 0.0, &
          0.0, 1.0, 0.0, &
          0.0, 0.0, 0.0 /)
          
          conma_elements = (/ &
          ! Capa 1 (dimension 1)
          0.0, 0.0, 1.0, &
          0.0, 0.0, 0.0, &
          1.0, 0.0, 0.0, &
          ! Capa 2 (dimension 2)
          0.0, 0.0, 0.0, &
          0.0, 0.0, 1.0, &
          0.0, 1.0, 0.0 /)
        
          reama_elements = (/ &
          ! Capa 1 (dimension 1:1)
          0.0, 0.0, 0.0, &
          0.0, 0.0, 0.0, &
          0.0, 0.0, 0.0 /)
          
          force_elements = (/ 1.0, 1.0, 0.0 /)
         
        case('Direct_Current_Electrical_Resistivity_in_2.5-D')
          difma_elements = (/ &
          ! Capa 1 (dimension 1:1)
          1.0, &
          ! Capa 2 (dimension 1:2)
          0.0, &
          ! Capa 3 (dimension 2:1)
          0.0, &
          ! Capa 4 (dimension 2:2)
          1.0 /)
          
          conma_elements = (/ &
          ! Capa 1 (dimension 1)
          0.0, &
          ! Capa 2 (dimension 2)
          0.0 /)
          
          reama_elements = (/ &
          ! Capa 1 (dimension 1:1)
          1.0 /)
          
          force_elements = (/ 0.0 /)
          
        case('Electric_Field_Excited_By_A_Double_Line_Source')
          difma_elements = (/ &
          ! Capa 1 (dimension 1:1)
          1.0, &
          ! Capa 2 (dimension 1:2)
          0.0, &
          ! Capa 3 (dimension 2:1)
          0.0, &
          ! Capa 4 (dimension 2:2)
          1.0 /)
          
          conma_elements = (/ &
          ! Capa 1 (dimension 1)
          0.0, &
          ! Capa 2 (dimension 2)
          0.0 /)
          
          reama_elements = (/ &
          ! Capa 1 (dimension 1:1)
          0.0 /)
          
          force_elements = (/ 0.0 /)
          
        case('Horizontal_Electric_Dipole_in_3-D_A_Wrong_capture_of_solution')
          difma_elements = (/ &
          ! Capa 1 (dimension 1:1)
          1.0 ,0.0, 0.0, &
          0.0 ,1.0, 0.0, &
          0.0 ,0.0, 1.0, &
          ! Capa 2 (dimension 1:2)
          0.0 ,1.0, 0.0, &
          -1.0,0.0, 0.0, &
          0.0 ,0.0, 0.0, &
          ! Capa 3 (dimension 2:1)
          0.0,-1.0, 0.0, &
          1.0, 0.0, 0.0, &
          0.0, 0.0, 0.0, &
          ! Capa 4 (dimension 2:2)
          1.0, 0.0, 0.0, &
          0.0, 1.0, 0.0, &
          0.0, 0.0, 1.0 /)
          
          conma_elements = (/ &
          ! Capa 1 (dimension 1)
          0.0, 0.0, 1.0, &
          0.0, 0.0, 0.0, &
          1.0, 0.0, 0.0, &
          ! Capa 2 (dimension 2)
          0.0, 0.0, 0.0, &
          0.0, 0.0, 1.0, &
          0.0, 1.0, 0.0 /)
         
          reama_elements = (/ &
          ! Capa 1 (dimension 1:1)
          0.0, 0.0, 0.0, &
          0.0, 0.0, 0.0, &
          0.0, 0.0, 0.0 /)
          
          force_elements = (/ 0.0, 0.0, 0.0 /)
          
        case('Transient_Electromagnetic_in_2.5-D')
          difma_elements = (/ &
          ! #DIFMA_xx
          1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 ,-1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 ,-1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,-1.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,-1.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0, &     
          ! #DIFMA_xz
          0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          ! #DIFMA_zx
          0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          ! #DIFMA_zz
          -1.0, 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , -1.0, 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , -1.0, 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,-1.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 /)
          
          conma_elements = (/ &
          ! #CONMAT_x
          0.0 , 0.0 , 0.0 , 1.0 , 0.0 ,-1.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 ,-1.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0, &
          1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0, &
          ! #CONMAT_z
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,-1.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 1.0 , 0.0 ,-1.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0, &
          0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 /)
          
          reama_elements = (/ &
          ! #REAMA
          1.0 ,  0.0 ,  0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 , -1.0 ,  0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,-1.0, &
          0.0 ,  0.0 ,  1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, &
          0.0 ,  0.0 ,  0.0 ,-1.0 , 0.0 ,-1.0 , 0.0 , 0.0, &
          0.0 ,  0.0 ,  0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0, &
          0.0 ,  0.0 ,  0.0 ,-1.0 , 0.0 ,-1.0 , 0.0 , 0.0, &
          0.0 ,  0.0 ,  0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0, &
          0.0 ,  1.0 ,  0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,-1.0 /)
          
          ! #FORCE
          force_elements = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
          
        case default
          write(*,*) 'At Tensors Module, No PDE Structure defined'
          write(*,*) 'Program being stopped'
          ! stop
          
      end select
      
      !Utilizar la función reshape para rellenar las matrices
      difma = reshape(difma_elements, shape(difma))
      conma = reshape(conma_elements, shape(conma))
      reama = reshape(reama_elements, shape(reama))
      force = reshape(force_elements, shape(force))
      
     
    end subroutine PDE_Coefficients
    
  !end contains


end module tensor_inputs
