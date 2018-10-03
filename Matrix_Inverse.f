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
 
      SUBROUTINE Matrix_Inverse(a,c,n)   
C     get inverse matrix
       
      !============================================================
      ! Inverse matrix
      ! Method: Based on Doolittle LU factorization for Ax=b
      ! Alex G. December 2009
      !-----------------------------------------------------------
      ! input ...
      ! a(n,n) - array of coefficients for matrix A
      ! n      - dimension
      ! output ...
      ! c(n,n) - inverse matrix of A
      ! comments ...
      ! the original matrix a(n,n) will be destroyed 
      ! during the calculation
      !===========================================================
      implicit none 
      integer,intent(in)::n
      double precision,intent(in):: a(n,n)
      double precision,intent(out):: c(n,n)
      double precision L(n,n), U(n,n), b(n), d(n), x(n)
      double precision coeff,tem_a(n,n)
      integer i, j, k
      
      tem_a = a
      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L=0.0D0
      U=0.0D0
      b=0.0D0

      ! step 1: forward elimination
      do k=1, n-1
          do i=k+1,n
              coeff=tem_a(i,k)/tem_a(k,k)
              L(i,k) = coeff
              do j=k+1,n
                  tem_a(i,j) = tem_a(i,j)-coeff*tem_a(k,j)          
              end do
          end do
      end do

      ! Step 2: prepare L and U matrices 
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
		  L(i,i) = 1.0D0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
          do i=1,j
              U(i,j) = tem_a(i,j)
          end do
      end do

      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
          b(k)=1.0D0
          d(1) = b(1)
		! Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
              d(i)=b(i)
              do j=1,i-1
                 d(i) = d(i) - L(i,j)*d(j)
              end do
          end do
		! Step 3b: Solve Ux=d using the back substitution
		  x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
              x(i) = d(i)
              do j=n,i+1,-1
			      x(i)=x(i)-U(i,j)*x(j)
              end do
          x(i) = x(i)/u(i,i)
          end do
		! Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
              c(i,k) = x(i)
          end do
		  b(k)=0.0D0
      end do
		  
      return
      END SUBROUTINE Matrix_Inverse
    


