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
 
      subroutine Tool_Cal_Intersection(A,B,C,D,X,Y,Yes_Cross)
C     get intersection point of two line segments

c     A: Coordinates of point A of line segment AB;
c     B: Coordinates of point B of line segment AB;
c     C: Coordinates of point C of line segment CD;
c     D: Coordinates of point D of line segment CD.

      use Global_Elem_Area
      
      implicit none
      double precision,intent(in)::A(2),B(2),C(2),D(2)
      double precision,intent(out)::X,Y
      logical,intent(out)::Yes_Cross
      
      double precision L_AB,L_CD,tem,k2,k1,b1,b2
      double precision L_AX,L_BX,L_CX,L_DX
      logical Yes_in_AB,Yes_in_CD
      double precision minX_AB,maxX_AB,minY_AB,maxY_AB,
     &                 minX_CD,maxX_CD,minY_CD,maxY_CD
     
c     The length of AB.
      L_AB = sqrt((A(1)-B(1))**2+(A(2)-B(2))**2)
c     The length of CD.
      L_CD = sqrt((C(1)-D(1))**2+(C(2)-D(2))**2)
      
      if (A(1).eq.B(1)) then
          X=A(1)
          tem = D(1)-C(1)
          k2=(D(2)-C(2))/tem
          b2=C(2)-k2*C(1) 
          Y=k2*X+b2;
      end if

      if (C(1).eq.D(1)) then
          X=C(1)
          tem =B(1)-A(1)
          k1=(B(2)-A(2))/tem
          b1=A(2)-k1*A(1)
          Y=k1*X+b1
      end if

      if ((A(1).ne.B(1)) .and. (C(1).ne.D(1))) then
          k1=(B(2)-A(2))/(B(1)-A(1))
          k2=(D(2)-C(2))/(D(1)-C(1))
          b1=A(2)-k1*A(1)
          b2=C(2)-k2*C(1)
          if (k1.eq.k2) then
              X=1.0D9
              Y=1.0D9
          else
              X=(b2-b1)/(k1-k2)
              Y=k1*X+b1
          end if
      end if

      minX_AB = min(A(1),B(1))
      maxX_AB = max(A(1),B(1))
      minY_AB = min(A(2),B(2))
      maxY_AB = max(A(2),B(2))
      minX_CD = min(C(1),D(1))
      maxX_CD = max(C(1),D(1))
      minY_CD = min(C(2),D(2))
      maxY_CD = max(C(2),D(2))

      Yes_Cross = .False.
      Yes_in_AB = .False.
      Yes_in_CD = .False.

c     Check if (X,Y) is in AB.
c     The length of AX.
      L_AX = sqrt((A(1)-X)**2+(A(2)-Y)**2)
c     The length of BX.
      L_BX = sqrt((B(1)-X)**2+(B(2)-Y)**2)
      if (((L_AX + L_BX)-L_AB) .le. 1.0D-6*sqrt(Ave_Elem_Area)) then
          Yes_in_AB = .True.
      end if


c     Check if (X,Y) is in CD.
c     The length of CX.
      L_CX = sqrt((C(1)-X)**2+(C(2)-Y)**2)
c     The length of DX.
      L_DX = sqrt((D(1)-X)**2+(D(2)-Y)**2)
      if (((L_CX + L_DX) - L_CD) .le. 1.0D-6*sqrt(Ave_Elem_Area)) then
          Yes_in_CD = .True.
      end if

c     if (X,Y) is in both AB and CD, then YesCross = 1.
      if ((Yes_in_AB .eqv..True.).and.(Yes_in_CD.eqv..True.)) then 
          Yes_Cross = .True.
      else
          Yes_Cross = .False.
      end if
      
      return 
      end SUBROUTINE Tool_Cal_Intersection                          
