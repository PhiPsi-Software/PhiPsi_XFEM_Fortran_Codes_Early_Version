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
 
      SUBROUTINE Matrix_Det(n,Matrix,det)   
C     get det of matrix
       
      implicit none
      integer n
      double precision,intent(in)::Matrix(n,n)
      double precision,intent(out)::det
      
      double precision m,temp
      double precision temp_Matrix(n,n)
      
      INTEGER :: i, j, k, l

      !Flag to know if the matrix is singular
      LOGICAL :: DetExists = .TRUE.
      !l stores the sign change of the matrix in case there are row exchanges
      l = 1
		  
      !====================================================================  
      !                  Convert to upper triangular form
      !====================================================================  
      temp_Matrix = Matrix
      
      !For each pivot row
      DO k = 1, n-1
          !If the pivot element is 0 then the matrix is probably singular 
          !We'll need to exchange this row with another one which
          !has a nonzero pivot.
          IF (temp_Matrix(k,k) .eq. 0.0D0) THEN 
			  DetExists = .FALSE.            
			  !This DO loop searches for the nearest row with a non zero pivot. A better
			  !way would be to look for the largest non-zero pivot.
              DO i = k+1, n 
				  !Check if the pivot is nonzero
				  !if it is, perform the row exchange
                  IF (temp_Matrix(i,k) /= 0.0D0) THEN
                      !This loop does the row exchange
                      DO j = 1, n 
                          temp = temp_Matrix(i,j)
                          temp_Matrix(i,j)= temp_Matrix(k,j)
                          temp_Matrix(k,j) = temp
                      END DO
                      !Because we were able to find a non-zero pivot, the matrix
                      !may be non-singular so set the flag
                      DetExists = .TRUE.
                      !We'll need to change the sign because of the row exchange 
                      !#Possible Bug#: I think the sign change should be outside 
                      !the if
                      l=-l
                      !Exit the outer Do loop because we found a row we can exchange with
                      EXIT
                  END IF
              END DO
			  !If no row was found, the matrix is singular and we can stop here.
              IF (DetExists .EQV. .FALSE.) THEN
				det = 0
				return
              END IF
          END IF
          !This Do loop performs the row operations as per the Gauss elimination algorithm.
          DO j = k+1, n
			  !Find the correct multiple for the current row
			  m = temp_Matrix(j,k)/temp_Matrix(k,k)
			  !Subtract the pivot row from the current row.
              DO i = k+1, n
                   temp_Matrix(j,i) = temp_Matrix(j,i) 
     &                             - m*temp_Matrix(k,i)
              END DO
          END DO
      END DO  
		  
      !====================================================================  
      !     Calculate determinant by finding product of diagonal elements
      !====================================================================  
      det = l
      DO i = 1, n
          det = det * temp_Matrix(i,i)
      END DO
		  
      return
      END SUBROUTINE Matrix_Det
    


