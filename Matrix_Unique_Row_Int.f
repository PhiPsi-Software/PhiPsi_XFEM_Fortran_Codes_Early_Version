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
 
      SUBROUTINE Matrix_Unique_Row_Int(m,n,m_Finish,Matrix,
     &                             Uniqued_Matrix,Uniqued_m,
     &                             Uni_Mat_Count)   
c     delete reduplicated rows of a matrix

      implicit none
      integer,intent(in):: m,n,m_Finish
      integer,intent(in):: Matrix(m,n)
      integer,intent(out):: Uniqued_m
      integer,intent(out):: Uniqued_Matrix(m,n)
      integer,intent(out):: Uni_Mat_Count(m)
      
      integer tem(n)
      integer i,j,k
      
      Uniqued_Matrix(1:m,1:n) = 0
      Uni_Mat_Count(1:m)=1
      
      k = 1
      Uniqued_Matrix(1,:) = Matrix(1,:)
      
      outer: do i=2,m_Finish
          do j=1,k
              tem = Uniqued_Matrix(j,:) - Matrix(i,:)
              if ((maxval(tem).eq. 0).and.
     &            (minval(tem).eq. 0)) then
                  !Found a match so start looking again
                  Uni_Mat_Count(j) = Uni_Mat_Count(j)+1
                  cycle outer
              end if
          end do
          !No match found so add it to the output
          k = k + 1
          Uniqued_Matrix(k,:) = Matrix(i,:)
      end do outer
                  
      Uniqued_m = k            

      return
      END SUBROUTINE Matrix_Unique_Row_Int
    


