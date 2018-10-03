!     .................................................
!             ____  _       _   ____  _____   _        
!            |  _ \| |     |_| |  _ \|  ___| |_|       
!            | |_) | |___   _  | |_) | |___   _        
!            |  _ /|  _  | | | |  _ /|___  | | |       
!            | |   | | | | | | | |    ___| | | |       
!            |_|   |_| |_| |_| |_|   |_____| |_|       
!     .................................................
!     PhiPsi:     a general-purpose computational      
!                 mechanics program written in Fortran.
!     Website:    http://phipsi.top                    
!     Author:     Shi Fang from Huaiyin Institute of   
!                 Technology, HuaiAn, JiangSu, China   
!     Contact me: shifang@hyit.edu.cn                  
!     ------------------------------------------------ 
!     Please cite the following papers:                
!     (1)Shi F, Wang X L, Liu C, Liu H, Wu H A. An     
!        XFEM-based method with reduction technique    
!        for modeling hydraulic fracture propagation   
!        in formations containing frictional natural   
!        fractures. Engineering Fracture Mechanics,    
!        2017, 173: 64-90.                             
!     (2)Shi F, Wang X L, Liu C, Liu H, Wu H A. A      
!        coupled extended finite element approach      
!        for modeling hydraulic fracturing in          
!        consideration of proppant. Journal of         
!        Natural Gas Science and Engineering, 2016,    
!        33: 885-897.                                  
!     (3)Shi F, Wang X L, Liu C, Liu H, Wu H A. An     
!        XFEM-based numerical model to calculate       
!        conductivity of propped fracture considering  
!        proppant transport, embedment and crushing.   
!        Journal of Petroleum Science and Engineering, 
!        2018, 167: 615-626..                          
 
      SUBROUTINE Get_Material_Matrix
c     Get the material matrix D

      use Global_Common
      use Global_Model
      use Global_Material
     
      implicit none
      integer i_mat
      double precision D1,D2,D3,c_E,c_v   
      
      print *,' '
      print *,'    Calculating material matrix D......'
      
      ALLOCATE(D(num_of_Material,3,3))
      ALLOCATE(S(num_of_Material,3,3))
      ALLOCATE(St(num_of_Material,2))
      ALLOCATE(Sc(num_of_Material,2))
      ALLOCATE(KIc(num_of_Material,2))
      ALLOCATE(E(num_of_Material,3))
      ALLOCATE(v(num_of_Material,3))
      ALLOCATE(thick(num_of_Material))
      ALLOCATE(density(num_of_Material))
      
      do i_mat =1,num_of_Material
          if (Material_Type(i_mat). eq. 1  )then      ! ISO material
              c_E  = Material_Para(i_mat,1)                                            ! Young's modulus
		      c_v  = Material_Para(i_mat,2)                                            ! Poisson's ratio
		      density(i_mat)   = Material_Para(i_mat,3)          
		      St(i_mat,1)      = Material_Para(i_mat,5)
		      Sc(i_mat,1)      = Material_Para(i_mat,7)
		      KIc(i_mat,1)     = Material_Para(i_mat,6)  
              
		      if(Key_Problem .eq. 1  )then                                              ! Plane stress
                  D1 = c_E/(1-c_v**2)                                                   ! Constant for elastic constant matrix
			      D2 = c_E*c_v/(1-c_v**2)                                               ! Constant for elastic constant matrix
			      D3 = c_E/2/(1+c_v)                                                    ! Constant for elastic constant matrix 
			      D(i_mat,1,:)  = [D1,   D2,   0.0D0]
                  D(i_mat,2,:)  = [D2,   D1,   0.0D0]
                  D(i_mat,3,:)  = [0.0D0,0.0D0,   D3]    
			      thick(i_mat)  = Material_Para(i_mat,4)                                  ! Plane stress thickness
		    elseif (Key_Problem .eq.2 )then                                               ! Plane strain
			      D1 = c_E*(1-c_v)/(1+c_v)/(1-2*c_v)                                      ! Constant for elastic constant matrix
			      D2 = c_E*c_v/(1+c_v)/(1-2*c_v)                                          ! Constant for elastic constant matrix
			      D3 = c_E/2/(1+c_v)                                                      ! Constant for elastic constant matrix
			      D(i_mat,1,:)  = [D1,   D2,   0.0D0]
                  D(i_mat,2,:)  = [D2,   D1,   0.0D0]
                  D(i_mat,3,:)  = [0.0D0,0.0D0,   D3]    
			      thick(i_mat)  = 1.0D0
		    end if                                  
            !S matrix       
            call Matrix_Inverse(D(i_mat,:,:),S(i_mat,:,:),3)                     
            
		    E(i_mat,1) = c_E
		    v(i_mat,1) = c_v
          end if      
      end do
      
      RETURN
      END SUBROUTINE Get_Material_Matrix
