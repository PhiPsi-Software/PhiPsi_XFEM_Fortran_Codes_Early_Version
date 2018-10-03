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
 
      SUBROUTINE Matrix_Solve_LSOE(c_Key_SLOE,K,F,D,n)
c     solve D=F/K.

      use Global_Common
      use Global_Model
       
      implicit none

      integer,intent(in)::n,c_Key_SLOE
      double precision,intent(in)::K(n,n),F(n)
      double precision,intent(out)::D(n)
      integer i,j,NZ_NUM,J_Sparse_K(n*n)
      double precision Sparse_K(n*n),D_Estimate(n)
      integer I_Map_Sparse_K(n+1),c_count,IWKSP(3*n), 
     &        IPARM(12),ITMAX,NW,IERR
      parameter(ITMAX=50000)      
      double precision RPARM(12)
      double precision WKSP(4*n + 4*ITMAX)  
      
      select case(c_Key_SLOE)
      case(1)

      case(2)

      case(3)

      case(4)
          !http://rene.ma.utexas.edu/CNA/ITPACK/
          print *, "          ****************************************" 
          print *, "               < LSOE Iterative Solver ITPACK >" 
          print *, "          ----------------------------------------" 
          print *, "          Preparing the sparse format matrix......" 
          NZ_NUM = 0       
          do i=1,n
              c_count = 0
              do j=1,n
                  if(K(i,j).ne.0.0D0) then
                      c_count = c_count +1
                      NZ_NUM             = NZ_NUM +1
                      Sparse_K(NZ_NUM)   = K(i,j)
                      J_Sparse_K(NZ_NUM) = j
                      if (c_count.eq.1) then
                          I_Map_Sparse_K(i)  = NZ_NUM
                      end if
                  end if
              end do
              I_Map_Sparse_K(i+1) = NZ_NUM+1
          end do
          
          print *, "          Iterative solving......"     
          CALL DFAULT (IPARM,RPARM) 
          IPARM(1)  = ITMAX
          IPARM(2)  = 1        
          IPARM(5)  = 1     
          IPARM(9)  = -1     
          IPARM(12) = 0    
          
          RPARM(1)  = 1.0D-9*max(Max_X_Coor-Min_X_Coor,
     &                           Max_Y_Coor-Min_Y_Coor)   
          
          NW = 4*n + 4*ITMAX   
          D_Estimate(1:n) = 0.0D+00    
          call JSI(n,I_Map_Sparse_K,J_Sparse_K(1:NZ_NUM),
     &               Sparse_K(1:NZ_NUM),F,
     &               D_Estimate,IWKSP,
     &               NW,WKSP,IPARM,RPARM,IERR) 
          D = D_Estimate       
          print *, "          ----------------------------------------" 
          print *, "                  < Exit Solver ITPACK >" 
          print *, "          ****************************************"          
      end select
      
      
      RETURN
      END SUBROUTINE Matrix_Solve_LSOE