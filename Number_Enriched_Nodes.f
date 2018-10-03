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
 
      SUBROUTINE Number_Enriched_Nodes(isub)
c     Number the enriched nodes.

      use Global_Crack
      use Global_Model
      use Global_Filename
      
      implicit none
      
      integer,intent(in)::isub
      integer i_C,i_N
      integer i,j
      character(200) c_File_name_1
      character(5) temp 
      
      print *,'    Numbering enriched nodes......'
      
      n_h_Node =0
      n_t_Node =0
      Total_Freedom = 0
      
      c_POS(1:Num_Node,1:num_Crack) = 0
      
      !Loop through each crack.
      do i_C = 1,num_Crack
          !Loop through each node.
          do i_N = 1,Num_Node
	          if (Enriched_Node_Type(i_N,i_C).eq.2) then       !Heaviside node
                  c_POS(i_N,i_C) = (Num_Node + n_h_Node +
     &                           n_t_Node*4 + n_j_Node) + 1
                  n_h_Node = n_h_Node + 1
		      elseif (Enriched_Node_Type(i_N,i_C) .eq.1)then   !Tip node	
                  c_POS(i_N,i_C) = (Num_Node + n_h_Node + 
     &                           n_t_Node*4 + n_j_Node) + 1
                  n_t_Node = n_t_Node + 1
		      end if
          end do
      end do     
      
      !Total degrees of freedom. 
      Total_Freedom = 2*(Num_Node + n_h_Node + n_t_Node*4)
      
 
      !save posi file
      print *,'    Saving ennd file......'
      write(temp,'(I4)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'posi'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown',form='unformatted')         
      do i=1,Num_Node
              write(101) (c_POS(i,j),j=1,num_Crack)
      end do
      close(101)       
                
      
      RETURN
      END SUBROUTINE Number_Enriched_Nodes
