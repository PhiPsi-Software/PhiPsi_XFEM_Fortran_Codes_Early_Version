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
 
      subroutine Tool_Yes_In_Poly(x,y,xpol,ypol,npol,Yes_INOUT)
C     check if a point inside a polygon
      integer i,j,npol
      double precision x,y,xpol(npol),ypol(npol)
      logical Yes_INOUT

      Yes_INOUT = .False.
      if((x.gt.maxval(xpol)) .or.
     &   (x.lt.minval(xpol)) .or.
     &   (y.gt.maxval(ypol)) .or.
     &   (y.lt.minval(ypol))) then
             Yes_INOUT = .False.
         goto 100
         
      else
          j = npol-1 
          do i=1,npol-1
              if ( ((ypol(i).gt.y).neqv. (ypol(j).gt.y)) .and.
     &               (x .lt. (xpol(j)-xpol(i)) * (y-ypol(i)) 
     &                / (ypol(j)-ypol(i)) + xpol(i)) ) then
                  Yes_INOUT = .NOT. Yes_INOUT
              end if
              j = i
          end do
      end if
      
  100 continue
  
      return 
      end SUBROUTINE Tool_Yes_In_Poly                          
