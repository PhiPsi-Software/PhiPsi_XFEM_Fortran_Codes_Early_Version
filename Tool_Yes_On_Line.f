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
 
      subroutine Tool_Yes_On_Line(x,y,A,B,Yes_ON)
C     判断某点P(x,y)是否在线段AB上

      use Global_Elem_Area
      
      implicit none
      double precision,intent(in):: x,y,A(2),B(2)
      logical,intent(out):: Yes_ON

      double precision L_AB,L_AP,L_BP
      
      
      Yes_ON = .False.

c     The length of AB.
      L_AB = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)

c     The length of AP.
      L_AP = sqrt((A(1)-x)**2+(A(2)-y)**2)

c     The length of BP.
      L_BP = sqrt((B(1)-x)**2+(B(2)-y)**2)


      if (((L_AP + L_BP)-L_AB) .le. 1.0D-6*sqrt(Ave_Elem_Area)) then
          Yes_ON = .True. 
      end if
      
C     print *,x,y
C     print *,A
C     print *,B
C     print *,Yes_ON
      return 
      end SUBROUTINE Tool_Yes_On_Line                          
