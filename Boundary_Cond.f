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
 
      SUBROUTINE Boundary_Cond(Total_FD,isub,freeDOF,num_freeDOF)
c     Deal with the boudary condition.

      use Global_Common
      use Global_Model
      use Global_Elem_Area
      use Global_Crack
      
      implicit none
      integer,intent(in)::isub
      integer,intent(in)::Total_FD
      integer,intent(out)::freeDOF(Total_FD),num_freeDOF
      integer i,cur_Node,i_fd
      integer c_index,fixedDOF(Total_FD)
      
      print *,'    Dealing with the boundary condition......'   
      
      !Initialize vector of fixed DOFs.
      freeDOF(1:Total_FD) = 0
      c_index = 0
      
      !Fix displacement in x-direction.
      do i = 1,Num_Bou_x
          cur_Node = Bou_x(i)
          c_index = c_index+1     
	      fixedDOF(c_index) = 2*cur_Node-1         
      end do   
      
      !Fix displacement in y-direction.
      do i = 1,Num_Bou_y
          cur_Node = Bou_y(i)
          c_index = c_index+1  
	      fixedDOF(c_index) = 2*cur_Node        
      end do  
      num_freeDOF = 0              
      do i_fd =1,Total_FD
          if ( .not.(ANY( fixedDOF(1:c_index) .eq. i_fd)) ) then
              num_freeDOF = num_freeDOF +1
              freeDOF(num_freeDOF) = i_fd
          end if
      end do
      
      RETURN
      END SUBROUTINE Boundary_Cond
