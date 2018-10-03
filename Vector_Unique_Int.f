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
 
      SUBROUTINE Vector_Unique_Int(n,n_Op_F,Vector,
     &                             Uniqued_Vec,Uniqued_n)   
c     delete reduplicated element of a vector

      implicit none
      integer,intent(in):: n,n_Op_F
      integer,intent(in):: Vector(n)
      integer,intent(out)::Uniqued_Vec(n)
      integer,intent(out)::Uniqued_n
      
      integer i,j,k
      
      k = 1
      Uniqued_Vec(1) = Vector(1)
          
      outer: do i=2,n_Op_F
          do j=1,k
              if (Uniqued_Vec(j) .eq. Vector(i)) then
                  ! Found a match so start looking again
                  cycle outer
              end if
          end do
          ! No match found so add it to the output
          k = k + 1
          Uniqued_Vec(k) = Vector(i)
      end do outer
             
      Uniqued_n   = k

      return
      END SUBROUTINE Vector_Unique_Int
    


