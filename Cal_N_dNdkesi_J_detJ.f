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
 
      subroutine Cal_N_dNdkesi_J_detJ(kesi,yita,
     &                                X_NODES,Y_NODES,
     &                                detJ,J,N,dNdkesi)
     
c     This function calculates N, dNdkesi, J and the determinant of Jacobian matrix.
      
      implicit none
      
      double precision,intent(in)::kesi,yita,X_NODES(4),Y_NODES(4)
      double precision,intent(out)::detJ,dNdkesi(4,2)
      double precision,intent(out):: N(2,8)
      double precision N1,N2,N3,N4
      double precision o     !zero
      double precision one   !1
      double precision fo    !4
      
      double precision temp(2,4),Coor(2,4),J(2,2) 
      
      o   = 0.0D0
      one = 1.0D0
      fo  = 4.0D0
      
	  ! Coordinates of the element.
	  Coor(1,:) = [X_NODES(1),X_NODES(2),X_NODES(3),X_NODES(4)]
	  Coor(2,:) = [Y_NODES(1),Y_NODES(2),Y_NODES(3),Y_NODES(4)]
			
	  ! Calculate N.
	  N1 = (one-kesi)*(one-yita)/fo
	  N2 = (one+kesi)*(one-yita)/fo
	  N3 = (one+kesi)*(one+yita)/fo
	  N4 = (one-kesi)*(one+yita)/fo
	  N(1,:) = [N1,o,N2,o,N3,o,N4,o]
	  N(2,:) = [o,N1,o,N2,o,N3,o,N4]

	! Calculate dNdkesi.
	  temp(1,:) =[(yita-one)/fo,(-yita+one)/fo,
     &                    (yita+one)/fo,(-yita-one)/fo]
      temp(2,:) =[(kesi-one)/fo,(-kesi-one)/fo,
     &                    (kesi+one)/fo,(-kesi+one)/fo]
              
	  dNdkesi = transpose(temp)

	  ! Calculate the Jacobian matrix.
	  J = MATMUL(Coor,dNdkesi)
      
	  ! Calculate the determinant of Jacobian matrix.
      call Matrix_Det(2,J,detJ)  
      
      return 
      end SUBROUTINE Cal_N_dNdkesi_J_detJ                     
