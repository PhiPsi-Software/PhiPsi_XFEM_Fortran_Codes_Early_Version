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
 
      subroutine Cal_Gauss_Points_QUAD(Num_GP,kesi,yita,weight)
C     ¸ßË¹µã
      
      implicit none
      integer,intent(in)::Num_GP
      double precision,intent(out)::kesi(Num_GP),
     &                              yita(Num_GP),
     &                              weight(Num_GP)
      integer i,j
      double precision np(8),nw(8)
      
      select case(Num_GP)
          case(4)
              np(1) = -0.577350269189626D0
	          np(2) =  0.577350269189626D0
	          nw(1) =  1.0D0
	          nw(2) =  1.0D0
	          do i = 1,2
	              do j = 1,2
	                  kesi((i-1)*2+j)   = np(i)
		              yita((i-1)*2+j)   = np(j)
                      weight((i-1)*2+j) = nw(i)*nw(j)
	              end do
              end do          
          case(9)
	          np(1) = -0.774596669241483D0
	          np(2) =  0.0D0
              np(3) =  0.774596669241483D0
	          nw(1) =  0.555555555555556D0
	          nw(2) =  0.888888888888889D0
	          nw(3) =  0.555555555555556D0
 	          do i = 1,3
	              do j = 1,3
	                  kesi((i-1)*3+j)   = np(i)
		              yita((i-1)*3+j)   = np(j)
                      weight((i-1)*3+j) = nw(i)*nw(j)
	              end do
              end do          
          case(16)
    	     np(1) = -0.861136311594053D0
    	     np(2) = -0.339981043584856D0
    	     np(3) =  0.339981043584856D0
    	     np(4) =  0.861136311594053D0
    	     nw(1) =  0.347854845137454D0
    	     nw(2) =  0.652145154862546D0
    	     nw(3) =  0.652145154862546D0
    	     nw(4) =  0.347854845137454D0
 	          do i = 1,4
	              do j = 1,4
	                  kesi((i-1)*4+j)   = np(i)
		              yita((i-1)*4+j)   = np(j)
                      weight((i-1)*4+j) = nw(i)*nw(j)
	              end do
              end do           
          case(25)
    	     np(1) = -0.906179845938664D0
    	     np(2) = -0.538469310105683D0
    	     np(3) =  0.000000000000000D0
    	     np(4) =  0.538469310105683D0
    	     np(5) =  0.906179845938664D0
    	     nw(1) =  0.236926885056189D0
    	     nw(2) =  0.478628670499366D0
    	     nw(3) =  0.568888888888889D0
    	     nw(4) =  0.478628670499366D0
    	     nw(5) =  0.236926885056189D0
 	          do i = 1,5
	              do j = 1,5
	                  kesi((i-1)*5+j)   = np(i)
		              yita((i-1)*5+j)   = np(j)
                      weight((i-1)*5+j) = nw(i)*nw(j)
	              end do
              end do           
          case(36)
             np(1) = -0.932469514203152D0
    	     np(2) = -0.661209386466265D0
    	     np(3) = -0.238619186083197D0
    	     np(4) =  0.238619186083197D0
    	     np(5) =  0.661209386466265D0
    	     np(6) =  0.932469514203152D0
    	     nw(1) =  0.171324492379170D0
    	     nw(2) =  0.360761573048139D0
    	     nw(3) =  0.467913934572691D0
    	     nw(4) =  0.467913934572691D0
    	     nw(5) =  0.360761573048139D0
    	     nw(6) =  0.171324492379170D0
 	          do i = 1,6
	              do j = 1,6
	                  kesi((i-1)*6+j)   = np(i)
		              yita((i-1)*6+j)   = np(j)
                      weight((i-1)*6+j) = nw(i)*nw(j)
	              end do
              end do 
          case(49)
    	     np(1) = -0.949107912342759D0
    	     np(2) = -0.741531185599394D0
    	     np(3) = -0.405845151377397D0
    	     np(4) =  0.0D0
    	     np(5) =  0.405845151377397D0
    	     np(6) =  0.741531185599394D0
    	     np(7) =  0.949109712342759D0
             
    	     nw(1) =  0.129484966168870D0
    	     nw(2) =  0.279705391489277D0
    	     nw(3) =  0.381830050505119D0
    	     nw(4) =  0.417959183673469D0
    	     nw(5) =  0.381830050505119D0
    	     nw(6) =  0.279705391489277D0
    	     nw(7) =  0.129484966168870D0
 	          do i = 1,7
	              do j = 1,7
	                  kesi((i-1)*7+j)   = np(i)
		              yita((i-1)*7+j)   = np(j)
                      weight((i-1)*7+j) = nw(i)*nw(j)
	              end do
              end do 
          case(64)
    	     np(1) = -0.960289856497536D0
    	     np(2) = -0.796666477413627D0
    	     np(3) = -0.525532409916329D0
    	     np(4) = -0.183434642495650D0
    	     np(5) =  0.183434642495650D0
    	     np(6) =  0.525532409916329D0
    	     np(7) =  0.796666477413627D0
    	     np(8) =  0.960289856497536D0
             
    	     nw(1) =  0.101228536290376D0
    	     nw(2) =  0.222381034453374D0
    	     nw(3) =  0.313706645877887D0
    	     nw(4) =  0.362683783378362D0
    	     nw(5) =  0.362683783378362D0
    	     nw(6) =  0.313706645877887D0
    	     nw(7) =  0.222381034453374D0
    	     nw(8) =  0.101228536290376D0
 	          do i = 1,8
	              do j = 1,8
	                  kesi((i-1)*8+j)   = np(i)
		              yita((i-1)*8+j)   = np(j)
                      weight((i-1)*8+j) = nw(i)*nw(j)
	              end do
              end do 
      end select
      
      return 
      end SUBROUTINE Cal_Gauss_Points_QUAD                          
