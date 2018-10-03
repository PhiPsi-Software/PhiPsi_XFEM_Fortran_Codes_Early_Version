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
 
      subroutine Cal_Shorten_or_Extend_Line(Line_AB,delta_L,
     &                                      Point_String,
     &                                      new_Line_AB,new_Point)
      ! Shorten or extend line_AB at point a or b by the increment of offset_L.
      ! delta_L can be negative.
      ! Point_String ='A' or 'B'.
      !                A                       B
      !
      !                ¡ñ----------------------¡ñ
      !
      ! Point_String ='A',delta_L < 0:
      !
      !                |<---delta_L--->|
      !                                ¡ñ------¡ñ
      !              
      !              
      ! Point_String ='B',delta_L > 0:
      !
      !                                        |<---delta_L--->|
      !                ¡ñ--------------------------------------¡ñ

      
      implicit none
      character*1,intent(in) :: Point_String
      double precision,intent(in) ::Line_AB(2,2),delta_L
      double precision,intent(out)::new_Line_AB(2,2),new_Point(2)
      
      double precision a_x,a_y,b_x,b_y,theta
      
      a_x = Line_AB(1,1)
      a_y = Line_AB(1,2)
      b_x = Line_AB(2,1)
      b_y = Line_AB(2,2)
			

      theta = atan2(b_y-a_y,b_x-a_x)
      
	  !Move C by offset_delta along the line_AB from B to A or from A to B.
      
      select case(Point_String)
      !Case 1: Point_String ='A'    
      case('A')
          new_Line_AB(1,1) = Line_AB(1,1)-delta_L*cos(theta)
          new_Line_AB(1,2) = Line_AB(1,2)-delta_L*sin(theta)
          new_Line_AB(2,:) = Line_AB(2,:)
          new_Point = [new_Line_AB(1,1),new_Line_AB(1,2)]
      !Case 2: Point_String ='B'
      case('B')  
          new_Line_AB(2,1) = Line_AB(2,1)+delta_L*cos(theta)
          new_Line_AB(2,2) = Line_AB(2,2)+delta_L*sin(theta)
          new_Line_AB(1,:) = Line_AB(1,:)
          new_Point = [new_Line_AB(2,1),new_Line_AB(2,2)]
      end select
      
      return 
      end SUBROUTINE Cal_Shorten_or_Extend_Line                          
