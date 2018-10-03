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
 
      subroutine Cal_Ele_Num_by_Coors(x,y,OUT_Elem)
C     get element number using coordinates.

      use Global_Model
      
      implicit none
      
      integer i
      double precision x,y
      integer,intent(out)::OUT_Elem
      double precision c_x_max,c_x_min,c_y_max,c_y_min
      
      integer c_count,Potent_Elem(5)
      logical Yes_In, Yes_On
      double precision xpol(5),ypol(5)
      
      c_count = 0
      
      !get the potential elements number.
      do i=1,Num_Elem
          c_x_max = maxval(G_X_NODES(i,:))
          c_x_min = minval(G_X_NODES(i,:))
          c_y_max = maxval(G_Y_NODES(i,:))
          c_y_min = minval(G_Y_NODES(i,:))
          if   ((x.le.c_x_max).and.(x.ge.c_x_min).
     &     and. (y.le.c_y_max).and.(y.ge.c_y_min))then
              c_count = c_count +1
              Potent_Elem(c_count) = i   
          end if
      end do 
      
      !Looking for the correct element.
      do i =1,c_count
          xpol(1:4) = G_X_NODES(Potent_Elem(i),1:4)
          ypol(1:4) = G_Y_NODES(Potent_Elem(i),1:4)
          xpol(5) = G_X_NODES(Potent_Elem(i),1)
          ypol(5) = G_Y_NODES(Potent_Elem(i),1)
          call Tool_Yes_In_Poly(x,y,xpol,ypol,5,Yes_In)
          call Tool_Yes_On_Poly(x,y,xpol,ypol,5,Yes_On)
          if ((Yes_In.eqv..True.).or.(Yes_On.eqv..True.))then
              OUT_Elem = Potent_Elem(i)
              exit
          end if
      end do
      
      return 
      end SUBROUTINE Cal_Ele_Num_by_Coors                          
