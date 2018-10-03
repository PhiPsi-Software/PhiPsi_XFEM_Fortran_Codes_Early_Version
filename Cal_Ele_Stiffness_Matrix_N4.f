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
 
      subroutine Cal_Ele_Stiffness_Matrix_N4(X_NODES,Y_NODES,
     &                           thick,c_D,kesi,yita,weight,
     &                           localK)
      !Calculate element stiffness matrix.
     
      implicit none
      
      double precision,intent(in)::thick,c_D(3,3),
     &                     kesi(4),yita(4),weight(4)
      double precision,intent(in):: X_NODES(4),Y_NODES(4)
      double precision,intent(out):: localK(8,8)
      
      integer i
      double precision X(4),Y(4),JM(4,4),a,b,c,d,
     &                Nkesi(4),Nyita(4),detJ,one,fo,
     &                B1(3,2),B2(3,2),B3(3,2),B4(3,2),ToTal_B(3,8)
      
      one = 1.0D0
      fo  = 4.0D0
      localK(1:8,1:8)=0.0D0

      X=[X_NODES(1),X_NODES(2),X_NODES(3),X_NODES(4)]
      Y=[Y_NODES(1),Y_NODES(2),Y_NODES(3),Y_NODES(4)]

      !Loop through each gauss point
      do i=1,4
	      JM(1,:)=[0.0D0,one-yita(i),yita(i)-kesi(i),  kesi(i)-one]
		  JM(2,:)=[yita(i)-one,0.0D0,one+kesi(i), -kesi(i)-yita(i)]
		  JM(3,:)=[kesi(i)-yita(i), -kesi(i)-one,0.0D0,yita(i)+one]
		  JM(4,:)=[one-kesi(i), yita(i)+kesi(i),-yita(i)-one,0.0D0]
	
	      Nkesi=[(yita(i)-one)/fo,(-yita(i)+one)/fo,
     &           (yita(i)+one)/fo,(-yita(i)-one)/fo]
	      Nyita=[(kesi(i)-one)/fo,(-kesi(i)-one)/fo,
     &           (kesi(i)+one)/fo,(-kesi(i)+one)/fo]
	
	      a=(Y(1)*(kesi(i)-one)+Y(2)*(-one-kesi(i))+
     *                     Y(3)*(one+kesi(i))+Y(4)*( one-kesi(i)))/fo
	      b=(Y(1)*(yita(i)-one)+Y(2)*(one-yita(i))+ 
     *                     Y(3)*(one+yita(i))+Y(4)*(-one-yita(i)))/fo
	      c=(X(1)*(yita(i)-one)+X(2)*(one-yita(i))+ 
     *                     X(3)*(one+yita(i))+X(4)*(-one-yita(i)))/fo
	      d=(X(1)*(kesi(i)-one)+X(2)*(-one-kesi(i))+
     *                     X(3)*(one+kesi(i))+X(4)*( one-kesi(i)))/fo

	      detJ  =   dot_product(MATMUL(X,JM),Y)/8.0D0

          B1(1,:)=[a*Nkesi(1)-b*Nyita(1),0.0D0 ]
	      B1(2,:)=[0.0D0, c*Nyita(1)-d*Nkesi(1)]
	      B1(3,:)=[c*Nyita(1)-d*Nkesi(1), a*Nkesi(1)-b*Nyita(1)]
	      B2(1,:)=[a*Nkesi(2)-b*Nyita(2),0.0D0 ]
	      B2(2,:)=[0.0D0 ,c*Nyita(2)-d*Nkesi(2)]
	      B2(3,:)=[c*Nyita(2)-d*Nkesi(2), a*Nkesi(2)-b*Nyita(2)]
	      B3(1,:)=[a*Nkesi(3)-b*Nyita(3),0.0D0]
	      B3(2,:)=[0.0D0 , c*Nyita(3)-d*Nkesi(3)]
	      B3(3,:)=[c*Nyita(3)-d*Nkesi(3), a*Nkesi(3)-b*Nyita(3)]
	      B4(1,:)=[a*Nkesi(4)-b*Nyita(4),0.0D0]
	      B4(2,:)=[0.0D0 ,c*Nyita(4)-d*Nkesi(4)]
	      B4(3,:)=[c*Nyita(4)-d*Nkesi(4),a*Nkesi(4)-b*Nyita(4)]
 	   
 	      ToTal_B(1:3,1:2)= B1
          ToTal_B(1:3,3:4)= B2
          ToTal_B(1:3,5:6)= B3
          ToTal_B(1:3,7:8)= B4
          ToTal_B=ToTal_B/detJ
	
 	      localK = localK + thick*weight(i)*detJ*
     &             MATMUL(MATMUL(transpose(ToTal_B),c_D),ToTal_B)
      end do
      
      return 
      end SUBROUTINE Cal_Ele_Stiffness_Matrix_N4               
