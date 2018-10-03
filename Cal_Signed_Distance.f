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
 
      subroutine Cal_Signed_Distance(Line_AB,Point_C,S_Distance)
c     This function calculates the signed distance from the Point_C to the Line_AB.
      
      implicit none
      double precision,intent(in)::Line_AB(2,2),Point_C(2)
      double precision,intent(out)::S_Distance
      
      double precision tem_1(2,2), tem_2(2)
      double precision tem_Det,tem_Norm
      
      tem_1(1,:) = Line_AB(2,:)-  Line_AB(1,:)
      tem_1(2,:) = Point_C     -  Line_AB(1,:)
      tem_2(:)   = Line_AB(2,:)-Line_AB(1,:)
      tem_Det = tem_1(1,1)*tem_1(2,2) - tem_1(2,1)*tem_1(1,2)
      
  
      call Vector_Norm2(2,tem_2,tem_Norm) 
      
      S_Distance = tem_Det / tem_Norm
      
      return 
      end SUBROUTINE Cal_Signed_Distance                          
