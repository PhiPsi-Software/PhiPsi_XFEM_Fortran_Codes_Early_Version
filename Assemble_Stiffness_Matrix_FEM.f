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
 
      SUBROUTINE Assemble_Stiffness_Matrix_FEM(isub,globalK)
c     Assemble the stiffness matrix.

      use Global_Model
      use Global_Filename
      use Global_Common
      use Global_Material
      
      implicit none
      integer,intent(in)::isub
      double precision globalK(Total_Freedom,Total_Freedom)
      integer i_E
      double precision c_thick,c_D(3,3)
      double precision c_X_NODES(4),c_Y_NODES(4)
      integer c_NN(4) 
      double precision kesi(4),yita(4),weight(4)        
      double precision localK(8,8)
      integer local(8),i_row,i_col,nIndex
      
      
      print *,'    Assembling the global stiffness matrix......'     
      
      globalK(1:Total_Freedom,1:Total_Freedom) = 0.0D0

      call Cal_Gauss_Points_QUAD(4,
     &                         kesi,
     &                         yita,
     &                         weight)
      nIndex = 0
      do i_E = 1,Num_Elem
          c_thick = thick(Elem_Mat(i_E))
          c_D     = D(Elem_Mat(i_E),:,:)     
          c_NN    = G_NN(i_E,:)
          c_X_NODES = G_X_NODES(i_E,:)
          c_Y_NODES = G_Y_NODES(i_E,:)    
	      !Traditional index locations
          local=[c_NN(1)*2-1,c_NN(1)*2,c_NN(2)*2-1,c_NN(2)*2,
     &           c_NN(3)*2-1,c_NN(3)*2,c_NN(4)*2-1,c_NN(4)*2]
          !Get the element stiffness matrix of the current element	
	      call Cal_Ele_Stiffness_Matrix_N4(c_X_NODES,c_Y_NODES,
     &                           c_thick,c_D,kesi,yita,weight,
     &                           localK)
          do i_row = 1,8
              do i_col = 1,8
                  globalK(local(i_row),local(i_col)) = 
     &                   globalK(local(i_row),local(i_col)) +
     &                   localK(i_row,i_col)
              end do
          end do
      end do
      
      RETURN
      END SUBROUTINE Assemble_Stiffness_Matrix_FEM
