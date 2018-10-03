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
 
      SUBROUTINE Matrix_Sort_Row_Int(m,n,Start_m,Finish_m,Matrix)   
C     sorting matrix                 
      implicit none
      integer,intent(in)::m,n,Start_m,Finish_m
      integer,intent(inout) :: Matrix(m,n)
      
      integer i,j,iii
      integer Tem_Vector(n),a
      
      do iii=Start_m,Finish_m
          Tem_Vector = Matrix(iii,:)
          do j=2, n
              a=Tem_Vector(j)
              do i=j-1,1,-1
                  if (Tem_Vector(i).le.a) goto 10
                  Tem_Vector(i+1)=Tem_Vector(i)
              end do
              i=0
   10         Tem_Vector(i+1) = a
          end do
          Matrix(iii,:) = Tem_Vector
      end do
      
      return
      END SUBROUTINE Matrix_Sort_Row_Int
    


