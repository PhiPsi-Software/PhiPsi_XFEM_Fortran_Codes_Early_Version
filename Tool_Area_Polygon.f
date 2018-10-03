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
 
      subroutine Tool_Area_Polygon(x,y,nb,area)
	  !get area of a polygon
      
      ! Code converted using TO_F90 by Alan Miller
      ! Date: 2000-07-04  Time: 12:24:06

      IMPLICIT NONE
      integer nb
      double precision x(nb)
      double precision y(nb)
      INTEGER i, n, nm1
      double precision a,area

      !*****************************************************************

      ! GIVEN A SEQUENCE OF NB POINTS (X(I),Y(I)),  polyarea COMPUTES THE AREA
      ! BOUNDED BY THE CLOSED POLYGONAL CURVE WHICH PASSES THROUGH THE POINTS IN
      ! THE ORDER THAT THEY ARE INDEXED.  THE FINAL POINT OF THE CURVE IS ASSUMED
      ! TO BE THE FIRST POINT GIVEN.  THEREFORE, IT NEED NOT BE LISTED AT THE END
      ! OF X AND Y.  THE CURVE IS NOT REQUIRED TO BE SIMPLE.  e.g. It may cross over
      ! itself.

      !*****************************************************************


      n = nb
      IF ((x(1) .eq. x(n)) .AND. (y(1) .eq. y(n))) n = n - 1

      SELECT CASE (n)
      CASE (:2)
          area = 0.0D0

      CASE (3)
          area=0.5D0*((x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1)))

      CASE DEFAULT
          nm1 = n - 1
          a = x(1)*(y(2) - y(n)) + x(n)*(y(1) - y(nm1))

          DO  i = 2, nm1
              a = a + x(i)*(y(i+1) - y(i-1))
          END DO
          
          area = 0.5D0*a
      END SELECT

      RETURN
      END subroutine Tool_Area_Polygon
 