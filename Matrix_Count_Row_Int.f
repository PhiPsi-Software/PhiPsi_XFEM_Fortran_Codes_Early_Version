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
 
      SUBROUTINE Matrix_Count_Row_Int(m,n,Matrix,Count_Row,Num_Once)   
C     count number of each row 
      
      implicit none
      integer i,j,m,n                
      integer Num_Once               
      integer Matrix(m,n),Count_Row(m)
      integer Vector_1(n),Vector_2(n)
      
      do i=1,m
          Count_Row(i)=1
          Vector_1 = Matrix(i,:)
          do j=1,m
              if(j.ne.i)then
                  Vector_2 = Matrix(j,:)
                  if ((MaxVal(Vector_1-Vector_2).eq.0).and.
     &                (MinVal(Vector_1-Vector_2).eq.0)) then
                      Count_Row(i) = Count_Row(i) + 1
                  end if
              end if
          end do
      end do
      
      Num_Once =0
      do i=1,m
          if (Count_Row(i).eq.1)then
              Num_Once = Num_Once +1
          end if
      end do   
      
      return
      END SUBROUTINE Matrix_Count_Row_Int
    


