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
 
      SUBROUTINE Location_Element_Stiff_Matrix(i_E,i_C,POS_c_Crack,
     &                                       Location_ESM_C_Crack,
     &                                       num_Loc_ESM_C_Crack)
c     get location of element stiffness matrix

      use Global_Crack
      use Global_Model
      use Global_Common
   
      implicit none
      
      integer,intent(in)::i_E,i_C,POS_c_Crack(Num_Node)
      integer,intent(out)::Location_ESM_C_Crack(40)    
      integer,intent(out)::num_Loc_ESM_C_Crack
      
      integer::Location_FEM(8)     
      double precision c_X_NODES(4),c_Y_NODES(4)
      integer c_NN(4),i_f,i_N      
      integer Enriched_Node_c_Crack(Num_Node)
      integer num_t,num_h,num_j,cnt
      integer:: Location_XFEM(32)     
      integer tem
      
      Enriched_Node_c_Crack = Enriched_Node_Type(:,i_C)
      
      c_NN    = G_NN(i_E,:)
      c_X_NODES = G_X_NODES(i_E,:)
      c_Y_NODES = G_X_NODES(i_E,:)
      
      !initial
      Location_FEM(1:8)=0
      Location_XFEM(1:32)=0
      Location_ESM_C_Crack(1:40)=0
      
      num_Loc_ESM_C_Crack = 0
      
      if (i_C .eq. 1)then
          do i_N =1,4
              Location_FEM(2*i_N-1) = 2*c_NN(i_N) - 1
              Location_FEM(2*i_N)   = 2*c_NN(i_N)
	      end do
      end if
      
      !If the element has no enriched nodes, then:
      if(maxval(Enriched_Node_c_Crack(c_NN)).eq.0 .and.
     &   minval(Enriched_Node_c_Crack(c_NN)).eq.0) then
          if (i_C .eq. 1)then
              Location_ESM_C_Crack(1:8) = Location_FEM
              num_Loc_ESM_C_Crack =8
          end if
      !If the element has enriched nodes, then:      
      else
          !Get the number of the tip enriched nodes of the element.
          num_t = count(Enriched_Node_c_Crack(c_NN).eq.1)
	      !Get the number of the Heaviside enriched nodes of the element.
          num_h = count(Enriched_Node_c_Crack(c_NN).eq.2)
	      tem = 2*(num_h*1 + num_t*4)            
          cnt = 0
          
          do i_N = 1,4
              if (Enriched_Node_c_Crack(c_NN(i_N)) .eq. 2  )then   ! Heaviside enriched node
			      cnt = cnt + 1
			      Location_XFEM(2*cnt-1) = 2*POS_c_Crack(c_NN(i_N)) - 1
			      Location_XFEM(2*cnt  ) = 2*POS_c_Crack(c_NN(i_N))
		       elseif( Enriched_Node_c_Crack(c_NN(i_N))  .eq.  1)then ! Tip enriched node
                   do i_f=1,4
                       cnt = cnt + 1
                       Location_XFEM(2*cnt-1) = 
     &                          2*(POS_c_Crack(c_NN(i_N))+i_f-1) - 1
                       Location_XFEM(2*cnt  ) = 
     &                          2*(POS_c_Crack(c_NN(i_N))+i_f-1)
			      end do
		       end if
          end do   
          if (i_C .eq. 1)then
              Location_ESM_C_Crack(1:8)     = Location_FEM
              Location_ESM_C_Crack(9:8+tem) = Location_XFEM(1:tem)
              num_Loc_ESM_C_Crack           = 8+tem      
          else
              Location_ESM_C_Crack(1:tem)   = Location_XFEM(1:tem)
              num_Loc_ESM_C_Crack           = tem              
          end if
      end if
      
      RETURN
      END SUBROUTINE Location_Element_Stiff_Matrix
