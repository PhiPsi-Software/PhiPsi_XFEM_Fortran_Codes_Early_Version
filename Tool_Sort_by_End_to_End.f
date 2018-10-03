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
 
      subroutine Tool_Sort_by_End_to_End(m,m_OP,Input_Outline,
     &                                      Output_Outline,cou)
C     sort by end to end

      implicit none
      integer i,j,m,m_OP,c_Outline,c_node_num_2,cou
      logical tem_Yes
      integer Input_Outline(m,2),Output_Outline(m,2)
      integer Location
      
      Output_Outline(1,:)  = Input_Outline(1,:)
      cou = 1
      c_Outline    = 1
      
      do i = 1,m_OP
          if (i.eq.1)then
              c_Outline = 1
	          c_node_num_2 = Input_Outline(1,2)
	      end if
          do j = 1,m
              if(j .ne. c_Outline) then
                if (c_node_num_2 .eq. Input_Outline(j,1)) then
                   call Vector_belongs_Matrix_Is_Int
     &                (cou,2,Output_Outline,
     &                [Input_Outline(j,1),Input_Outline(j,2)],
     &                Location,tem_Yes)
                   if(tem_Yes.eqv..False.) then
                       cou = cou + 1
					   Output_Outline(cou,:) = 
     &                          [Input_Outline(j,1),Input_Outline(j,2)]
					   c_Outline = j
					   c_node_num_2 = Input_Outline(j,2)

                       exit
                   end if
     
                else if(c_node_num_2 .eq. Input_Outline(j,2)) then
                    call Vector_belongs_Matrix_Is_Int
     &                (cou,2,Output_Outline,
     &                [Input_Outline(j,2),Input_Outline(j,1)],
     &                Location,tem_Yes)
                    if(tem_Yes.eqv..False.)then
                       cou = cou + 1
                       Output_Outline(cou,:)=
     &                          [Input_Outline(j,2),Input_Outline(j,1)]
					   c_Outline = j
					   c_node_num_2 = Input_Outline(j,1)
                       exit
	                end if
	           end if
	         end if
	     end do
      end do

      
      RETURN
      END subroutine Tool_Sort_by_End_to_End