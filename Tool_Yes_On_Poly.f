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
 
      subroutine Tool_Yes_On_Poly(x,y,xpol,ypol,npol,Yes_ONOUT)
C     check if a point on the boundary of a polygon

      integer i,npol
      double precision x,y,xpol(npol),ypol(npol)
      logical,intent(out)::Yes_ONOUT

      logical Yes_ON
      

      Yes_ONOUT = .False.
      
      do i=1,npol-1 
          call Tool_Yes_On_Line(x,y,
     &                  [xpol(i),ypol(i)],
     &                  [xpol(i+1),ypol(i+1)],Yes_ON)
          if (Yes_ON.eqv..True.) then
              Yes_ONOUT = .True.
              exit
          endif

      end do 
      
      return 
      end SUBROUTINE Tool_Yes_On_Poly                          
