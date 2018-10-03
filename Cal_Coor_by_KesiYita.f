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
 
      subroutine Cal_Coor_by_KesiYita(kesi,yita,X_NODES,Y_NODES,
     &                                Out_x,Out_y)
C     Calculate global coordinates using local coordinates kesi and yita.
      
      implicit none
      double precision,intent(in)::kesi,yita,X_NODES(4),Y_NODES(4)
      double precision,intent(out)::Out_x,Out_y
      double precision N(4)
      
      N  = 0.25D0 * [(1-kesi)*(1-yita),(1+kesi)*(1-yita),
     &               (1+kesi)*(1+yita),(1-kesi)*(1+yita)]
      
      Out_x  = N(1)*X_NODES(1)+N(2)*X_NODES(2)+
     &         N(3)*X_NODES(3)+N(4)*X_NODES(4)
      Out_y  = N(1)*Y_NODES(1)+N(2)*Y_NODES(2)+
     &         N(3)*Y_NODES(3)+N(4)*Y_NODES(4)
     
      return 
      end SUBROUTINE Cal_Coor_by_KesiYita                          
