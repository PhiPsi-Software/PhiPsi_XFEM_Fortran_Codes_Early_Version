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
 
      subroutine Cal_B_Matrix(kesi,yita,
     &                        i_C,i_E,i_G,
     &                        c_NN,c_X_NODES,c_Y_NODES,
     &                        tem_B,num_tem_B)    
c     calculate B matrix.

      use Global_Crack
      use Global_Model
      use Global_Filename
      use Global_Common
      use Global_Material
      
      implicit none
      
      integer,intent(in)::i_C,i_E,i_G
      double precision,intent(in)::c_X_NODES(4),c_Y_NODES(4)
      integer,intent(in)::c_NN(4)
      double precision,intent(in)::kesi,yita
      double precision,intent(out)::tem_B(3,40)
      integer,intent(out)::num_tem_B        
      double precision detJ, J(2,2), Inverse_J(2,2)
      double precision N(2,8),dNdkesi(4,2),dNdx(4,2),Coor_AB(2,2)
      integer mat_num,c_mat_type
      double precision B_FEM(3,8),B_XFEM(3,32),BI_enr(3,2)
      integer num_B_XFEM
      integer i_N,ref_elem,i_F
      logical Flag_Check_Kink_Tip,Yes_Gauss_in_BCD,Yes_in_KAB
      double precision distance_Vertex,Global_coor_Gauss(2),
     &           Tri_KAB(3,2),Closed_Tri_KAB(4,2),
     &           H _i,omega, Coor_Tip(2),Trans_Matrix(2,2),
     &           H_Vertex,distance_Gauss,H_gp,distance_Node
      integer  i_tem
      double precision  J_gp, J_i,c_Segment(2)
      integer find_element(8)    
      double precision Global_coor_Gauss_Doc(2),xp_Gauss(2),r_Gauss,
     &                 theta_Gauss,xp_Node(2),r_Node,theta_Node
      double precision dFdx(4),dFdy(4),F_Gauss(4),F_Node(4)  
      double precision c_N(4),aa,bb,B_enr(3,2),BI_enr_Tip(3,8)
      integer num_BI_enr_Tip
          
      call Cal_N_dNdkesi_J_detJ(kesi,yita,
     &                           c_X_NODES,c_Y_NODES,
     &                           detJ,J,N,dNdkesi)    
     
      call Matrix_Inverse(J,Inverse_J,2) 

      dNdx = MATMUL(dNdkesi,Inverse_J)
      
      Global_coor_Gauss(1) = DOT_PRODUCT(N(1,1:7:2),c_X_NODES(1:4))
      Global_coor_Gauss(2) = DOT_PRODUCT(N(1,1:7:2),c_Y_NODES(1:4))      
      
      mat_num    = Elem_Mat(i_E)
      
      ! Material type of the current element.
      c_mat_type = Material_Type(mat_num)     
      
      !B matrix for conventional element. 
      B_FEM(1:3,1:8) = 0.0D0
      
      if (i_C.eq.1) then
          B_FEM(1,1:8:2)   =  dNdx(:,1)
	      B_FEM(2,2:8:2)   =  dNdx(:,2)
	      B_FEM(3,1:8:2)   =  dNdx(:,2)
	      B_FEM(3,2:8:2)   =  dNdx(:,1)
      end if
      
      !If the element has no enriched nodes, then:
      if(maxval(Enriched_Node_Type(c_NN,i_C)).eq.0 .and.
     &   minval(Enriched_Node_Type(c_NN,i_C)).eq.0) then
          if (i_C.eq.1) then
              tem_B(1:3,1:8) = B_FEM
              num_tem_B = 8
          else
              num_tem_B = 0
          end if
      !If the element has enriched nodes, then:      
      else 
          B_XFEM(1:3,1:32) = 0.0D0
          num_B_XFEM = 0
          do i_N = 1,4
	          !Initialize flag of kinked crack of the tip segment.
              Flag_Check_Kink_Tip =  .False.
              Yes_Gauss_in_BCD    =  .False.
	          !-----------------------------------------------------------
		      !----------- Case 1: Heaviside enriched node ---------------
		      !-----------------------------------------------------------   
              if (Enriched_Node_Type(c_NN(i_N),i_C).eq.2)then     ! Heaviside enriched node   
                  ref_elem = 0
                  !if 2:fully cracked element without kink point or 3:fully cracked element with kink point, then:
		          if ((Elem_Type(i_E,i_C) .eq. 2 ).or. 
     &                                (Elem_Type(i_E,i_C) .eq. 3)) then
                      ref_elem = i_E
                      !Coordinates of the 2 intersections(A,B) of the crack segment and these 4 element(fully cracked) lines.
                      Coor_AB(1,:) = [Coors_Element_Crack(ref_elem,1),
     &                               Coors_Element_Crack(ref_elem,2)]
                      Coor_AB(2,:) = [Coors_Element_Crack(ref_elem,3),
     &                               Coors_Element_Crack(ref_elem,4)]  
                      if  (Elem_Type(i_E,i_C) .eq. 3) then
					      !Calculates the distance of the kink point(K) to the line(AB) composed of the 2 intersections(A,B).
                          call Cal_Signed_Distance(Coor_AB,
     &                                Coors_Vertex(ref_elem,:),
     &                                distance_Vertex)
                          call Cal_Sign(distance_Vertex,H_Vertex)
                          call Cal_Signed_Distance(Coor_AB,
     &                                Global_coor_Gauss,
     &                                distance_Gauss)
                          call Cal_Sign(distance_Gauss,H_gp)
                          
					      if (H_Vertex*H_gp <= 0) then
					          H_gp = H_gp
                          else
					          !The triangle KAB.
					          Tri_KAB(1,:) =Coor_AB(1,:)
                              Tri_KAB(2,:) =Coor_AB(2,:)
                              Tri_KAB(3,:) =Coors_Vertex(ref_elem,:)
                              Closed_Tri_KAB(1:3,:) =  Tri_KAB
                              Closed_Tri_KAB(4,:) =  Tri_KAB(1,:)
                              call Tool_Yes_In_Poly
     &                       (Global_coor_Gauss(1),Global_coor_Gauss(2),
     &                         Closed_Tri_KAB(:,1), Closed_Tri_KAB(:,2),
     &                         4,Yes_in_KAB)
     
						      if (Yes_in_KAB .eqv. .True.) then
						          H_gp = -H_gp
                              end if
                          end if       
                          call Cal_Signed_Distance(Coor_AB,
     &                                [c_X_NODES(i_N),c_Y_NODES(i_N)],
     &                                distance_Node)
                          call Cal_Sign(distance_Node,H_i)                   
                      else
                          call Cal_Signed_Distance(Coor_AB,
     &                                Global_coor_Gauss,
     &                                distance_Gauss)
                          call Cal_Sign(distance_Gauss,H_gp)
                          call Cal_Signed_Distance(Coor_AB,
     &                                [c_X_NODES(i_N),c_Y_NODES(i_N)],
     &                                distance_Node)
                          call Cal_Sign(distance_Node,H_i)
                      end if
                      !Calculate enriched B matrix.
				      BI_enr(1,:) = [dNdx(i_N,1)*(H_gp-H_i), 0.0D0]
                      BI_enr(2,:) = [0.0D0,dNdx(i_N,2)*(H_gp-H_i)]
                      BI_enr(3,:) = [dNdx(i_N,2)*(H_gp-H_i), 
     &                               dNdx(i_N,1)*(H_gp-H_i)]
                  else
                      BI_enr(1:3,1:2)=0.0D0
                  end if
			      !Update B_XFEM.
                  B_XFEM(1:3,num_B_XFEM+1:num_B_XFEM+2) = BI_enr
                  num_B_XFEM = num_B_XFEM + 2                 
                  
       	      !-----------------------------------------------------------
		      !----------- Case 2: Junction enriched node ----------------
		      !-----------------------------------------------------------	   
              elseif (Enriched_Node_Type(c_NN(i_N),i_C).eq.3)then ! Junction enriched node

              !-----------------------------------------------------------
		      !----------- Case 3: Tip enriched node ---------------------
		      !-----------------------------------------------------------
              elseif (Enriched_Node_Type(c_NN(i_N),i_C).eq.1)then ! Tip enriched node
                  ref_elem = 0
                  !Find the fully tip enriched element which the current tip enriched node attached to. 
			      !Fully tip enriched element. "ref_elem" means the element which contains the crack tip.              
                  if (Elem_Type(i_E,i_C) .eq. 1) then
                      ref_elem = i_E       ! The tip element.   
                  else
                      !Partly tip enriched element.
                      find_element = Node_Elements(c_NN(i_N),:)
                      
                      do i_tem=1,num_Node_Elements(c_NN(i_N))
                          if (Elem_Type(find_element(i_tem),i_C) .eq.1) 
     &                                                         then
                              !The current node belongs to the tip element.
                              ref_elem = find_element(i_tem)
                              Flag_Check_Kink_Tip = .True.
                          end if
                      end do
                  end if
			      !Coordinates of the 2 intersections(A,B) of the main crack segment and these 4 element(fully cracked) lines
			      !and the coordinates of the tip.
                  Coor_AB(1,:) = [Coors_Element_Crack(ref_elem,1),
     &                            Coors_Element_Crack(ref_elem,2)]
                  Coor_AB(2,:) = [Coors_Element_Crack(ref_elem,3),
     &                            Coors_Element_Crack(ref_elem,4)]  
                  c_Segment = Coor_AB(2,:) - Coor_AB(1,:)     
			      omega   = atan2(c_Segment(2),c_Segment(1))            
			      Coor_Tip = [Coor_AB(2,1),Coor_AB(2,2)]
			      !The transformation matrix from the local tip coordinates to the global coordinates. 
                  Trans_Matrix(1,:) = [  cos(omega),sin(omega)]
			      Trans_Matrix(2,:) = [ -sin(omega),cos(omega)]
                  
                  !Calculate the value(r,theta) of the gauss point in the local tip coordinates.
			      if (Yes_Gauss_in_BCD.eqv..False.) then
                      xp_Gauss= MATMUL(Trans_Matrix,
     &                           Global_coor_Gauss-Coor_Tip)
                      r_Gauss = sqrt(xp_Gauss(1)**2+xp_Gauss(2)**2)
                      theta_Gauss= atan2(xp_Gauss(2),xp_Gauss(1))
                  else
                      r_Gauss     = sqrt(Global_coor_Gauss_Doc(1)**2 + 
     &                                   Global_coor_Gauss_Doc(2)**2)
                      theta_Gauss = atan2(Global_coor_Gauss_Doc(2),
     &                              Global_coor_Gauss_Doc(1))
                  end if      
                  !Check theta.
                  if ((theta_Gauss >pi  .or. theta_Gauss < -pi)) then
                      print *,'    ********************************'
     &                             //'******************************' 
                      print *, '    Error :: Angle is wrong when' 
     &                          // ' calculates r and theta!'
                      print *,'    ********************************'
     &                             //'******************************'     
                  end if  
                  !Calculate the tip enrichment functions F and their derivative, dFdx and dFdy.          
                  call Cal_F_dFdx_dFdy(r_Gauss,theta_Gauss,omega,
     &                                 c_mat_type,
     &                                 F_Gauss,dFdx,dFdy)
		          !Calculate the value(r,theta) of node in the local tip coordinates.
                  xp_Node= MATMUL(Trans_Matrix,
     &                        [c_X_NODES(i_N),c_Y_NODES(i_N)]-Coor_Tip)
                  r_Node = sqrt(xp_Node(1)**2+xp_Node(2)**2)
                  theta_Node= atan2(xp_Node(2),xp_Node(1))            
                  call Cal_F(r_Node,theta_Node,omega,
     &                                 c_mat_type,F_Node)     
                  BI_enr(1,:)  = [dNdx(i_N,1)*(J_gp-J_i),     0.0D0]
                  BI_enr(2,:)  = [0.0D0,      dNdx(i_N,2)*(J_gp-J_i)]
                  BI_enr(3,:)  = [dNdx(i_N,2)*(J_gp-J_i),
     &                            dNdx(i_N,1)*(J_gp-J_i)]
                  num_BI_enr_Tip = 0
                  do i_F =1,4
                      c_N(1) = N(1,1)
                      c_N(2) = N(1,3)
                      c_N(3) = N(1,5)
                      c_N(4) = N(1,7)     
			          aa = dNdx(i_N,1)*(F_Gauss(i_F)-F_Node(i_F)) 
     &                                   + c_N(i_N)*dFdx(i_F)
			          bb = dNdx(i_N,2)*(F_Gauss(i_F)-F_Node(i_F)) 
     &                                   + c_N(i_N)*dFdy(i_F)
                      B_enr(1,:) = [aa,  0.0D0]
                      B_enr(2,:) = [0.0D0,  bb]
                      B_enr(3,:) = [bb,     aa]
                       BI_enr_Tip(1:3,num_BI_enr_Tip+1:num_BI_enr_Tip+2) 
     &                                                           = B_enr  
                       num_BI_enr_Tip = num_BI_enr_Tip +2
                  end do
                  !Update B_XFEM.
                  B_XFEM(1:3,num_B_XFEM+1:num_B_XFEM+8) = BI_enr_Tip
                  num_B_XFEM = num_B_XFEM + 8    
              end if
          end do
          !B matrix of enriched element.
          if (i_C.eq.1) then
              tem_B(1:3,1:8)            = B_FEM
              tem_B(1:3,9:8+num_B_XFEM) = B_XFEM
              num_tem_B = 8 + num_B_XFEM
          else
              tem_B(1:3,1:+num_B_XFEM) = B_XFEM
              num_tem_B = num_B_XFEM
          end if
      end if   

      RETURN
      END SUBROUTINE Cal_B_Matrix