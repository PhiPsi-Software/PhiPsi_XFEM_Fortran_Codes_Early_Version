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
 
      SUBROUTINE Assemble_Stiffness_Matrix_XFEM(isub,globalK)
c     Assemble the stiffness matrix.

      use Global_Crack
      use Global_Model
      use Global_Filename
      use Global_Common
      use Global_Material
      
      implicit none
      integer,intent(in)::isub
      double precision globalK(Total_Freedom,Total_Freedom)
      integer i_C,i_E,i_G
      double precision c_thick,c_D(3,3)
      double precision c_X_NODES(4),c_Y_NODES(4)
      integer c_NN(4)
      double precision kesi_Enr(Num_Gauss_Points),         !增强单元的高斯点数据，积分方案2
     &                 yita_Enr(Num_Gauss_Points),
     &                 weight_Enr(Num_Gauss_Points)    
      double precision kesi_N_Enr(4),                      !非增强单元的高斯点数据，积分方案2
     &                 yita_N_Enr(4),
     &                 weight_N_Enr(4)    
      double precision,ALLOCATABLE:: kesi(:),yita(:),weight(:)         
      
      integer:: Location_ESM(200)    !40 *5, 每个单元最多40个增强自由度，每个单元假设最多跟5条裂纹有关
      integer::Location_ESM_C_Crack(40)
      double precision B(3,40),tem_B(3,40)
      integer num_B,num_tem_B
      integer num_Loc_ESM
      integer num_Loc_ESM_C_Crack
      integer c_Num_Gauss_Point
      double precision detJ
      
      print *,'    Assembling the global stiffness matrix......'     
      
      globalK(1:Total_Freedom,1:Total_Freedom) = 0.0D0

      call Cal_Gauss_Points_QUAD(Num_Gauss_Points,
     &                         kesi_Enr,
     &                         yita_Enr,
     &                         weight_Enr)
      call Cal_Gauss_Points_QUAD(4,
     &                         kesi_N_Enr,
     &                         yita_N_Enr,
     &                         weight_N_Enr)
      do i_E = 1,Num_Elem
          c_thick = thick(Elem_Mat(i_E))
          c_D     = D(Elem_Mat(i_E),:,:)     
          c_NN    = G_NN(i_E,:)
          c_X_NODES = G_X_NODES(i_E,:)
          c_Y_NODES = G_Y_NODES(i_E,:)
          !if the current element are enriched element, then 8x8 gauss points is suggested:
          if (maxval(Enriched_Node_Type(c_NN,:)).ne.0)then
              ALLOCATE(kesi(Num_Gauss_Points))
              ALLOCATE(yita(Num_Gauss_Points))
              ALLOCATE(weight(Num_Gauss_Points))
              kesi   = kesi_Enr
              yita   = yita_Enr
              weight = weight_Enr
              c_Num_Gauss_Point = Num_Gauss_Points
          !if the current element are not enriched element, then 2x2 gauss points:
          else 
              ALLOCATE(kesi(4))
              ALLOCATE(yita(4))
              ALLOCATE(weight(4))
              kesi   = kesi_N_Enr
              yita   = yita_N_Enr
              weight = weight_N_Enr
              c_Num_Gauss_Point = 4
          end if 
          
          !Decide the location of each element stiffness matrix in the global stiffness matrix.              
          Location_ESM(1:200)  = 0           
          num_Loc_ESM = 0                  
          do i_C =1,num_Crack 
              call Location_Element_Stiff_Matrix(i_E,i_C,
     &                                     c_POS(:,i_C),
     &                                      Location_ESM_C_Crack,
     &                                      num_Loc_ESM_C_Crack)
              
              Location_ESM(num_Loc_ESM+1:
     &                   num_Loc_ESM+num_Loc_ESM_C_Crack) = 
     &                   Location_ESM_C_Crack(1:num_Loc_ESM_C_Crack)
              num_Loc_ESM  =  num_Loc_ESM + num_Loc_ESM_C_Crack                     
          end do
          
          !loop through gauss points.
          do i_G = 1,c_Num_Gauss_Point
              B(1:3,1:40) = 0.0D0
              num_B = 0
              call Cal_detJ(kesi(i_G),yita(i_G),
     &                                     c_X_NODES,c_Y_NODES,
     &                                     detJ)    
              !Calculate the B Matrix, Loop through each crack.
              do i_C =1,num_Crack 
                  call Cal_B_Matrix(kesi(i_G),yita(i_G),
     &                                i_C,i_E,i_G,
     &                                c_NN,c_X_NODES,c_Y_NODES,
     &                                tem_B,num_tem_B)
                  
                  B(1:3,num_B+1:num_B+num_tem_B) =  
     &                                      tem_B(1:3,1:num_tem_B)
                  num_B = num_B + num_tem_B                        
              end do
              !Assemble the global stiffness matrix.    
              globalK(Location_ESM(1:num_Loc_ESM),
     &                   Location_ESM(1:num_Loc_ESM)) = 
     &        globalK(Location_ESM(1:num_Loc_ESM),
     &                   Location_ESM(1:num_Loc_ESM)) + 
     &                     c_thick*weight(i_G)*detJ*
     &                       MATMUL(MATMUL(transpose
     &                           (B(1:3,1:num_B)),c_D),
     &                            B(1:3,1:num_B))   
          end do                 
          DEALLOCATE(kesi)
          DEALLOCATE(yita)
          DEALLOCATE(weight)           
      end do          
      
      RETURN
      END SUBROUTINE Assemble_Stiffness_Matrix_XFEM
