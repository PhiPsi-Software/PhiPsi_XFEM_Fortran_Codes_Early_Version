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
 
      SUBROUTINE Force_Vector(Total_FD,isub,Lambda,globalF)
c     Force vector.


      use Global_Common
      use Global_Model
      use Global_Elem_Area 
      implicit none

      integer,intent(in)::isub,Total_FD
      double precision,intent(in)::Lambda
      double precision,intent(out)::globalF(Total_FD)
      
      double precision h,density
      integer i,i_E,i_N,cur_Node,mat_num,c_mat_type,c_node
      
      print *,'    Constructing global force vector......'   
      
      !Thickness.
      if (Key_Problem ==1) then               ! Plane stress.
          if (Material_Type(1) ==1) then      ! ISO material.
              h = Material_Para(1,4)          ! Plane stress thickness.
	      elseif (Material_Type(1) ==2) then  ! Orthotropic material.
	          h = Material_Para(1,6)          ! Plane stress thickness.
	      end if
      elseif (Key_Problem ==2) then           ! Plane strain.
          h =1.0D0
      end if
      
      !Initialize the global force vector. 
      globalF(1:Total_FD) = 0.0D0
      
      do i = 1,Num_Foc_x
	      cur_Node = int(Foc_x(i,1))
          globalF(2*cur_Node-1) = Lambda*Foc_x(i,2)                   
      end do

      do i = 1,Num_Foc_y
	      cur_Node = int(Foc_y(i,1))
	      globalF(2*cur_Node) =   Lambda*Foc_y(i,2)                  
      end do
      
      !Gravity.
      if(Key_Gravity==1) then
          do i_E = 1,Num_Elem
              !Material number of the current element.
	          mat_num    = Elem_Mat(i_E)
	          !Material type of the current element.
	          c_mat_type = Material_Type(mat_num)
              if (c_mat_type ==1) then       ! ISO material.
                  density = Material_Para(mat_num,3)
              elseif (c_mat_type ==2 .or. c_mat_type ==3)  then! Orthotropic material.
			      density = Material_Para(mat_num,5); 
              end if
      
	          do i_N = 1,4
                  c_node = G_NN(i_E,i_N) 
                  globalF(2*c_node)  = globalF(2*c_node)-
     &                          9.8D0*density*h*Elem_Area(i_E)/4
              end do
	      end do
      end if

      RETURN
      END SUBROUTINE Force_Vector
