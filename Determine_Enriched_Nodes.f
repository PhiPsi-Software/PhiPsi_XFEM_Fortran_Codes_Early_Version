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
 
      SUBROUTINE Determine_Enriched_Nodes(isub)
c     Determine enriched nodes.

      use Global_Common   
      use Global_Filename
      use Global_Model
      use Global_Elem_Area
      use Global_Crack
      
      implicit none
      integer isub,i_C,i_S,i_O,i_Side,i_E,i_N,c_Edge_Elem
      logical Yes_in_Outline1,Yes_in_Outline2
      double precision crack_p1(2),crack_p2(2),
     &                 Outline_p1(2),Outline_p2(2)
      integer EndByEnd_Outline(size(Outline,1)+1) 
      double precision c_X,c_Y,c_Line_AB(2,2),shorted_Line_AB(2,2),
     &                 new_Point(2)
      logical c_Yes_Cross,Flag_1,Flag_2
      integer N1,N2,N3,N4,NN(4),NN_L(5),n_node_1,n_node_2,Num_Intersect
      double precision point_side_1(2),point_side_2(2),
     &                 X_el_inter,Y_el_inter,crack_inter(5,2),
     &                 point_seg_1(2),point_seg_2(2),
     &                 X_NODES_L(5),Y_NODES_L(5)
      integer point_seg_inElemnt_Flag(Num_Elem)
      double precision point_seg_inElemnt(Num_Elem,2)
      integer i_inter,i_Segment,Crack_Seg_Num(5),
     &        Uniqued_Cr_Seg_Num(5),Uniqued_Num_Cr_Seg
      double precision x_inter,y_inter
      logical Yes_on_Segment,Yes_On_Line
      integer j_C,j_S,tem_i_1,int_Flag_1,int_Flag_2
      logical Yes_Ap_has_GP,Yes_Am_has_GP,Yes_Cutthrough
      integer DOMAIN_Outline(20,2),m_DOMAIN_Outline,
     &        Domain_El(10),n_Domain_El
      double precision tem_01,AB(2,2),CD(2,2),tem_02(2),tem_03(2),
     &                 tem_point_seg_1(2),tem_point_seg_2(2)
      
c     ------------------------------           
c     Variables specification:
c     ------------------------------     
C     Elem_Type(Num_Elem,num_crack);     ! Element type: 1--Tip element;
C                                        !               2--Fully cracked element without kink point;
C                                        !               3--Fully cracked element with kink point;
C									     !               4--Junction element.
C     Enriched_Node_Type(Num_Node,num_crack); ! Type of enriched nodes:1--Tip enrichment;
C                                             !                        2--Heaviside enrichment;
C										      !                        3--Junction enrichment.
C
C     Coors_Element_Crack(Num_Elem,4);        ! Coors of intersections of crack segment and sides of 
C                                             ! element, and coors of crack tip, 2 points at most.
C     Coors_Tip(Num_Elem,2);                  ! Coors of crack tip(1 at most) in each element
C     Coors_Vertex(Num_Elem,2);               ! Coors of vertex(kink point, 1 at most) of crack segment in each element
C     Coors_Junction(Num_Elem,4);             ! Coors of intersections of crack segment and sides of 
C                                             ! element, and coors of junction point, 2 points at most
C												
C     Crack_Tip_Type(num_crack,2);            ! Store the type of each tip of each crack, for example, junction 
C                                             ! crack tip, 1; tip near the edge of the model,-1;Edge crack, -2.
C     x_cr_tip_nodes(num_crack,Num_Node);     ! x coordinate of crack tip of tip enriched nodes.		
C     y_cr_tip_nodes(num_crack,Num_Node);     ! y coordinate of crack tip of tip enriched nodes.	
C     Flag_Node_Friction(num_crack,Num_Node); ! Flag: =1 if the enriched node is Frictional crack node.		

c     ------------------------------      
c     initialization.
c     ------------------------------  
      Elem_Type(1:Num_Elem,1:Max_Num_Crack) = 0
      Enriched_Node_Type(1:Num_Node,1:Max_Num_Crack) = 0
      Coors_Element_Crack(1:Num_Elem,1:4) = 0.0D0
      Coors_Tip(1:Num_Elem,1:2) = 0.0D0
      Coors_Vertex(1:Num_Elem,1:2) = 0.0D0
      Coors_Junction(1:Num_Elem,1:4) = 0.0D0
      x_cr_tip_nodes(1:Max_Num_Crack,1:Num_Node) = 0.0D0
      y_cr_tip_nodes(1:Max_Num_Crack,1:Num_Node) = 0.0D0
      Crack_Tip_Type(1:Max_Num_Crack,1:2) = 0
      Edge_Disposed_Crack(1:Max_Num_Crack,1:Max_Num_Cr_P,1:2) = 0.0D0
      
 1001 FORMAT('     Caution :: tip 1 of crack ',I3,' is edge crack.') 
 1002 FORMAT('     Caution :: tip 2 of crack ',I3,' is edge crack.')     
 
      print *,'    Determine enriched node......'
      
c     ----------------------------------------------------------------------------------
c     ---------------------   Step 1: determine the element type   ---------------------
c     ----------------------------------------------------------------------------------
      Edge_Disposed_Crack = Crack_Coor
      do i_C=1,num_Crack
          c_Edge_Elem = 0  
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !loop through each crack segment
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i_S = 1,Each_Cr_Poi_Num(i_C)-1
              !******************************************
              !coordinates of the current crack segment
              !******************************************
              crack_p1 = [Crack_Coor(i_C,i_S,1),Crack_Coor(i_C,i_S,2)]
              crack_p2 = [Crack_Coor(i_C,i_S+1,1),
     &                    Crack_Coor(i_C,i_S+1,2)]
              EndByEnd_Outline(1:size(Outline,1)) = 
     &                                 Outline(:,1)
              EndByEnd_Outline(size(Outline,1)+1) = 
     &                                 Outline(1,1)        
              !***************************************
              !check if point in EndByEnd_Outline
              !***************************************
              Call Tool_Yes_In_Poly
     &             (crack_p1(1),crack_p1(2),
     &              Coor(EndByEnd_Outline,1),Coor(EndByEnd_Outline,2),
     &              size(EndByEnd_Outline),Yes_in_Outline1)
              Call Tool_Yes_In_Poly
     &             (crack_p2(1),crack_p2(2),
     &              Coor(EndByEnd_Outline,1),Coor(EndByEnd_Outline,2),
     &              size(EndByEnd_Outline),Yes_in_Outline2)
              !*************************************************************************
              !find out the ntersections of the current crack segment has with outline
              !*************************************************************************
              do i_O = 1,size(Outline,1)
                 Outline_p1 = Coor(Outline(i_O,1),:)
                 Outline_p2 = Coor(Outline(i_O,2),:)   
                 call Tool_Cal_Intersection(crack_p1,crack_p2,
     &                                      Outline_p1,Outline_p2,
     &                                c_X,c_Y,c_Yes_Cross)  
                 !*************************************************************
                 !if the current crack segment has intersections with outline
                 !*************************************************************
                  if(c_Yes_Cross.eqv..true.) then
                      call Cal_Ele_Num_by_Coors(c_X,c_Y,c_Edge_Elem)
                      !#########################################
                      !the first crack segment
                      !#########################################
                      if ((i_S .eq. 1).and.
     &                    (Yes_in_Outline1.eqv..False.))then
                          !if crack tip is edge crack
                          WRITE(*,1001) i_C
                          Crack_Tip_Type(i_C,1) = -2
                          Flag_Crack_Tip_Out_Mol(i_C,1) = .True.
                          c_Line_AB(1,:) = [c_X,c_Y]
                          c_Line_AB(2,:) = [crack_p2]
                          tem_01 = -1.0D-3*Ave_Elem_L
                          call Cal_Shorten_or_Extend_Line
     &                              (c_Line_AB,tem_01,'A',
     &                              shorted_Line_AB,new_Point)
                          Edge_Disposed_Crack(i_C,1,:)=new_Point
                          !Get the intersection point of the shorted first crack segment 
                          !'shorted_Line_AB' and the 4 element sides.
      					  N1  = Elem_Node(c_Edge_Elem,1)                                         
					      N2  = Elem_Node(c_Edge_Elem,2)                                             
					      N3  = Elem_Node(c_Edge_Elem,3)                                             
					      N4  = Elem_Node(c_Edge_Elem,4)                                            
					      NN_L= [N1,N2,N3,N4,N1]
					      do i_Side = 1,4
					          n_node_1 = NN_L(i_Side)
                              n_node_2 = NN_L(i_Side+1)
						      !Get the coordinate of the two points of the current element side.
						      point_side_1 = [Coor(n_node_1,1),Coor(n_node_1,2)]
						      point_side_2 = [Coor(n_node_2,1),Coor(n_node_2,2)]
                              !Check if they intersects and calculate intersection point.
                              call Tool_Cal_Intersection(
     &                        shorted_Line_AB(1,:),shorted_Line_AB(2,:),
     &                        point_side_1,point_side_2,
     &                        c_X,c_Y,c_Yes_Cross) 
                              if (c_Yes_Cross.eqv..True.) then
                                  X_el_inter = c_X
                                  Y_el_inter = c_Y
                                  exit
                              end if
                          end do  
                          Coors_Element_Crack(c_Edge_Elem,:) 
     &                               = [c_X,c_Y,X_el_inter,Y_el_inter]
                      !##########################
                      !last crack segment
                      !##########################
                      elseif ((i_S .eq. Each_Cr_Poi_Num(i_C)-1).and.
     &                    (Yes_in_Outline2.eqv..False.))then
                          !if crack tip is edge crack
                          WRITE(*,1002) i_C
                          Crack_Tip_Type(i_C,2) = -2
                          Flag_Crack_Tip_Out_Mol(i_C,2) = .True.
                          c_Line_AB(1,:) = [crack_p1]
                          c_Line_AB(2,:) = [c_X,c_Y]
                          tem_01 = -1.0D-3*Ave_Elem_L
                          call Cal_Shorten_or_Extend_Line
     &                              (c_Line_AB,tem_01,'B',
     &                              shorted_Line_AB,new_Point)
                          Edge_Disposed_Crack(i_C,
     &                                 Each_Cr_Poi_Num(i_C),:)=new_Point
                          !Get the intersection point of the shorted first crack segment 
                          !'shorted_Line_AB' and the 4 element sides.
      					  N1  = Elem_Node(c_Edge_Elem,1)                                         
					      N2  = Elem_Node(c_Edge_Elem,2)                                             
					      N3  = Elem_Node(c_Edge_Elem,3)                                             
					      N4  = Elem_Node(c_Edge_Elem,4)                                            
					      NN_L= [N1,N2,N3,N4,N1]
					      do i_Side = 1,4
					          n_node_1 = NN_L(i_Side)
                              n_node_2 = NN_L(i_Side+1)
						      !Get the coordinate of the two points of the current element side.
						      point_side_1 = [Coor(n_node_1,1),Coor(n_node_1,2)]
						      point_side_2 = [Coor(n_node_2,1),Coor(n_node_2,2)]
                              !Check if they intersects and calculate intersection point.
                              call Tool_Cal_Intersection(
     &                        shorted_Line_AB(1,:),shorted_Line_AB(2,:),
     &                        point_side_1,point_side_2,
     &                        c_X,c_Y,c_Yes_Cross) 
                              if (c_Yes_Cross.eqv..True.) then
                                  X_el_inter = c_X
                                  Y_el_inter = c_Y
                                  exit
                              end if
                          end do  
                          Coors_Element_Crack(c_Edge_Elem,:) 
     &                               = [c_X,c_Y,X_el_inter,Y_el_inter] 
                      end if
                  end if
                  
              end do
          end do
          
          !!!!!!!!!!!!!!!!!!!!!!!!
          !loop through elements
          !!!!!!!!!!!!!!!!!!!!!!!!         
          do i_E = 1,Num_Elem
              N1  = Elem_Node(i_E,1)                                         
              N2  = Elem_Node(i_E,2)                                             
              N3  = Elem_Node(i_E,3)                                             
              N4  = Elem_Node(i_E,4)  
              NN  = [N1,N2,N3,N4]
              NN_L= [N1,N2,N3,N4,N1]
              X_NODES_L = Coor(NN_L,1)
              Y_NODES_L = Coor(NN_L,2)
              
              Num_Intersect = 0
      		  Flag_1 = .False.         !Flag if point_seg_1 is inside the current element. 
		      Flag_2 = .False.         !Flag if point_seg_2 is inside the current element.
              
              !***************************
              !loop through crack sgements
              !***************************
              do i_S = 1,Each_Cr_Poi_Num(i_C)-1
                  !Get the coordinate of the two points of the current segment.
                  point_seg_1=[Crack_Coor(i_C,i_S,1),
     &                         Crack_Coor(i_C,i_S,2)]
                  point_seg_2=[Crack_Coor(i_C,i_S+1,1),
     &                         Crack_Coor(i_C,i_S+1,2)] 
                  do i_Side = 1,4
				      n_node_1 = NN_L(i_Side)
                      n_node_2 = NN_L(i_Side+1)
				      !Get the coordinate of the two points of the current side.
			          point_side_1 = [Coor(n_node_1,1),Coor(n_node_1,2)]
			          point_side_2 = [Coor(n_node_2,1),Coor(n_node_2,2)]            
                      !Check if they intersects and calculate intersection point.
                      call Tool_Cal_Intersection(
     &                point_seg_1,point_seg_2,
     &                point_side_1,point_side_2,
     &                c_X,c_Y,c_Yes_Cross) 
                      if(c_Yes_Cross.eqv..True.)then
                          Num_Intersect  = Num_Intersect + 1                        
                          crack_inter(Num_Intersect,:) = [c_X,c_Y]
       					  !Check if point_seg_1 is inside the current element.
                          call Tool_Yes_In_Poly
     &                              (point_seg_1(1),point_seg_1(2),
     &                               X_NODES_L,Y_NODES_L,5,Flag_1)
					      !Check if point_seg_2 is inside the current element.
                          call Tool_Yes_In_Poly
     &                              (point_seg_2(1),point_seg_2(2),
     &                               X_NODES_L,Y_NODES_L,5,Flag_2)

                          !check if a crack segment inside a element
					      if ((Flag_1.eqv..True.) .and. (Flag_2.eqv..True.)) then
                              print *,"   "
                              print *,"   ?????????????????????????????"
                              print *, "   Error :: A "
     &                            // "complete crack segment"
     &                            // " was found in only one element."
                              print *,"   ?????????????????????????????"
                              print *,"   Notice :: Message produced in"
     &                        // "SUBROUTINE Determine_Enriched_Nodes."
                              print *,"   ?????????????????????????????"
                              print *,"   "
					      end if
 					      !Get the coordinates of the point of the crack segment 
                          !in the element, actually, it's point_seg_1 or point_seg_2. 
                          if (Flag_1.eqv..True.) then
                              int_Flag_1 = 1
                          else
                              int_Flag_1 = 0
                          end if
                          if (Flag_2.eqv..True.) then
                              int_Flag_2 = 1
                          else
                              int_Flag_2 = 0
                          end if
					      point_seg_inElemnt(i_E,:) = 
     &                          point_seg_1*int_Flag_1 + 
     &                          point_seg_2*int_Flag_2
					      if (Flag_1.eqv..True.) then
					          point_seg_inElemnt_Flag(i_E) = 1
					      elseif (Flag_2.eqv..True.) then
					          point_seg_inElemnt_Flag(i_E) = 2
					      end if
                      end if
                  end do
              end do

              !##################################################
              !(Type 2):Fully cracked element without kink point.
              !##################################################
              if((Num_Intersect.eq.2).
     &            and.(Flag_1.eqv..False.).and.(Flag_2.eqv..False.))then
                  Elem_Type(i_E,i_C) = 2
                  !Check if the order of the two points of "crack_inter" need to be changed.
                  !Make sure that the order of "crack_inter" are the same as the crack.
			      !Loop to find how many crack segments have "crack_inter" points, actually 1.
                  tem_i_1 = 0
                  Crack_Seg_Num(1:5) =0
                  do i_inter=1,2
                      x_inter = crack_inter(i_inter,1)
                      y_inter = crack_inter(i_inter,2)
                      do i_Segment = 1,size(Crack_Coor,2)-1
                          point_seg_1 = Crack_Coor(i_C,i_Segment,:)
                          point_seg_2 = Crack_Coor(i_C,i_Segment+1,:)
                          call Tool_Yes_On_Line(x_inter,y_inter,
     &                          point_seg_1,point_seg_2,Yes_on_Segment)
					      if (Yes_on_Segment) then
                                tem_i_1 = tem_i_1 + 1
                                Crack_Seg_Num(tem_i_1)  = i_Segment
					      end if
                      end do                     
                  end do
                  call Vector_Unique_Int(5,tem_i_1,Crack_Seg_Num,
     &                            Uniqued_Cr_Seg_Num,Uniqued_Num_Cr_Seg) 
                  point_seg_1 = Crack_Coor(i_C,Uniqued_Cr_Seg_Num(1),:)
			      point_seg_2 = Crack_Coor(i_C,Uniqued_Cr_Seg_Num(1)+1,:)
			      AB(1,:) = point_seg_1
                  AB(2,:) = point_seg_2
			      CD(1,:) = [crack_inter(1,1),crack_inter(1,2)]
                  CD(2,:) = [crack_inter(2,1),crack_inter(2,2)]
                  tem_02  = AB(2,:)-AB(1,:)
                  tem_03  = CD(2,:)-CD(1,:)
                  if(dot_product(tem_02,tem_03) .GT. 0.0D0) then
                       Coors_Element_Crack(i_E,:) = 
     &                         [crack_inter(1,1),crack_inter(1,2),
     &                          crack_inter(2,1),crack_inter(2,2)]
                  else
                  !Change the order of "crack_inter" to make sure that the order are the same as the crack.
                        Coors_Element_Crack(i_E,:) = 
     &                         [crack_inter(2,1),crack_inter(2,2),
     &                          crack_inter(1,1),crack_inter(1,2)]                     
                  end if
              end if
              !
              !##################################################
              !(Type 3):Fully cracked element with kink point.    
              !##################################################
              if((Num_Intersect.eq.2).and.((Flag_1).or.(Flag_2)))then       
                  Elem_Type(i_E,i_C) = 3  
                  ! Check if the order of the two points of crack_inter need to be changed.
			      ! Make sure that the order of "crack_inter" are the same as the crack.
			      ! Loop to find how many crack segments have "crack_inter" points, actually 2.
                  tem_i_1 = 0
                  Crack_Seg_Num(1:5) =0
                  do i_inter=1,2
                      x_inter = crack_inter(i_inter,1)
                      y_inter = crack_inter(i_inter,2)
                      do i_Segment = 1,Each_Cr_Poi_Num(i_C)-1
                          point_seg_1 = Crack_Coor(i_C,i_Segment,:)
                          point_seg_2 = Crack_Coor(i_C,i_Segment+1,:)
                          call Tool_Yes_On_Line(x_inter,y_inter,
     &                          point_seg_1,point_seg_2,Yes_on_Segment)
					      if (Yes_on_Segment) then
                                tem_i_1 = tem_i_1 + 1
                                Crack_Seg_Num(tem_i_1)  = i_Segment
					      end if
                      end do                     
                  end do
                  call Vector_Unique_Int(5,tem_i_1,Crack_Seg_Num,
     &                            Uniqued_Cr_Seg_Num,Uniqued_Num_Cr_Seg) 
                  point_seg_1 = Crack_Coor(i_C,Uniqued_Cr_Seg_Num(1),:)
			      point_seg_2 = Crack_Coor(i_C,Uniqued_Cr_Seg_Num(2)+1,:)
			      AB(1,:) = point_seg_1
                  AB(2,:) = point_seg_2
			      CD(1,:) = [crack_inter(1,1),crack_inter(1,2)]
                  CD(2,:) = [crack_inter(2,1),crack_inter(2,2)]
                  tem_02  = AB(2,:)-AB(1,:)
                  tem_03  = CD(2,:)-CD(1,:)
                  if(dot_product(tem_02,tem_03) .GT. 0.0D0) then
                       Coors_Element_Crack(i_E,:) = 
     &                         [crack_inter(1,1),crack_inter(1,2),
     &                          crack_inter(2,1),crack_inter(2,2)]
                  else
                  !Change the order of "crack_inter" to make sure that the order are the same as the crack.
                        Coors_Element_Crack(i_E,:) = 
     &                         [crack_inter(2,1),crack_inter(2,2),
     &                          crack_inter(1,1),crack_inter(1,2)]                     
                  end if   
                  Coors_Vertex(i_E,:)  = point_seg_inElemnt(i_E,:)
                  !print *,point_seg_inElemnt(i_E,:)
              end if
              !##################################################
              !(Type 1):Tip element or (Type 4):Junction element.    
              !##################################################
              if(Num_Intersect.eq.1)then 
                  Elem_Type(i_E,i_C) = 1
                  !The coordinates of the crack tip.
			      Coors_Tip(i_E,:) = point_seg_inElemnt(i_E,:)
                  !Crack tip of tip enriched nodes.	
                  x_cr_tip_nodes(i_C,N1) = Coors_Tip(i_E,1)
                  x_cr_tip_nodes(i_C,N2) = Coors_Tip(i_E,1)
                  x_cr_tip_nodes(i_C,N3) = Coors_Tip(i_E,1)
                  x_cr_tip_nodes(i_C,N4) = Coors_Tip(i_E,1)
                  y_cr_tip_nodes(i_C,N1) = Coors_Tip(i_E,2)
                  y_cr_tip_nodes(i_C,N2) = Coors_Tip(i_E,2)
                  y_cr_tip_nodes(i_C,N3) = Coors_Tip(i_E,2)
                  y_cr_tip_nodes(i_C,N4) = Coors_Tip(i_E,2)
              
                  do j_C=1,num_Crack
                      if(j_C.ne.i_C) then
                          !Loop through each segment of j_C. 
                          do j_S = 1,Each_Cr_Poi_Num(j_C)-1
                              !Get the coordinate of the two points of the current segment.
                              tem_point_seg_1 = Crack_Coor(j_C,j_S,:)
                              tem_point_seg_2 = Crack_Coor(j_C,j_S+1,:)

                              !Check that if the crack tip is on the segment or not.
                              call Tool_Yes_On_Line
     &                            (Coors_Tip(i_E,1),Coors_Tip(i_E,2),
     &                             tem_point_seg_1,tem_point_seg_2,
     &                             Yes_On_Line)
                              !If yes, then:Junction enriched elements
						      if(Yes_On_Line) then
						          Elem_Type(i_E,i_C) = 4
						      end if
                              
                          end do
                      end if
                  end do
			      !tip enriched elements.
			      if (Elem_Type(i_E,i_C) .ne. 4) then
                      Elem_Type(i_E,i_C) = 1 
                      Coors_Element_Crack(i_E,1:2) = crack_inter(1,:) 
                      Coors_Element_Crack(i_E,3:4) = 
     &                              [Coors_Tip(i_E,1),Coors_Tip(i_E,2)]
				      !If there is edge crak,then the crack tip should be changed to an edge.
                      if (i_E .eq. c_Edge_Elem) then
                          Elem_Type(c_Edge_Elem,i_C) = 2
                      end if
                      
                  end if
              end if
          end do
      end do
      

c     ----------------------------------------------------------------------------------
c     ---------------------   Step 2: determine the enriched nodes  --------------------
c     ----------------------------------------------------------------------------------
      do i_C=1,num_Crack       
          do i_E = 1,Num_Elem
              N1  = Elem_Node(i_E,1)                                         
              N2  = Elem_Node(i_E,2)                                             
              N3  = Elem_Node(i_E,3)                                             
              N4  = Elem_Node(i_E,4)  
              NN  = [N1,N2,N3,N4]
              NN_L= [N1,N2,N3,N4,N1]
              X_NODES_L = Coor(NN_L,1)
              Y_NODES_L = Coor(NN_L,2)
              !If type-1 element, i.e. tip element, then:
		      if (Elem_Type(i_E,i_C) .eq. 1) then
		          Enriched_Node_Type(NN,i_C) =1                !Tip nodes  
		      !If type-2 element, i.e. fully cracked element without kink point, then:
		      elseif (Elem_Type(i_E,i_C) .eq. 2) then
                  !Loop through each node of the current element.
                  do i_N = 1,4
			          !If the node has not been enriched yet, then:
                      if (Enriched_Node_Type(NN(i_N),i_C).eq. 0) then
				          !Get the support domain of the current node.
					      call Cal_Support_Domain_of_Node(NN(i_N),
     &                                  DOMAIN_Outline,
     &                                  m_DOMAIN_Outline,
     &                                  Domain_El,
     &                                  n_Domain_El)
					      !Get the A+ and A- Domain.
                          call Cal_Ap_and_Am(Crack_Coor(i_C,:,:),
     &                         Each_Cr_Poi_Num(i_C),
     &                         DOMAIN_Outline(1:m_DOMAIN_Outline,:),
     &                         m_DOMAIN_Outline,
     &                         Domain_El(1:n_Domain_El),n_Domain_El,
     &                         i_E,NN(i_N),
     &                         Yes_Cutthrough,
     &                         Yes_Ap_has_GP,Yes_Am_has_GP)
					      !If both A_Plus and A_Minus has gauss points, then:
					      if (Yes_Ap_has_GP .and. Yes_Am_has_GP) then
					          Enriched_Node_Type(NN(i_N),i_C) = 2   ! Heaviside nodes
					      end if
                      end if
			      end do  
              !If type-3 element, i.e. Fully cracked element with kink point, then:
		      elseif (Elem_Type(i_E,i_C) .eq. 3) then
                  !Loop through each node of the current element.
                  do i_N = 1,4
			          !If the node has not been enriched yet, then:
                      if (Enriched_Node_Type(NN(i_N),i_C).eq. 0) then
                          !print *,i_E,i_C
				          !Get the support domain of the current node.
					      call Cal_Support_Domain_of_Node(NN(i_N),
     &                                  DOMAIN_Outline,
     &                                  m_DOMAIN_Outline,
     &                                  Domain_El,
     &                                  n_Domain_El)
					      !Get the A+ and A- Domain.
                          call Cal_Ap_and_Am(Crack_Coor(i_C,:,:),
     &                         Each_Cr_Poi_Num(i_C),
     &                         DOMAIN_Outline(1:m_DOMAIN_Outline,:),
     &                         m_DOMAIN_Outline,
     &                         Domain_El(1:n_Domain_El),n_Domain_El,
     &                         i_E,NN(i_N),
     &                         Yes_Cutthrough,
     &                         Yes_Ap_has_GP,Yes_Am_has_GP)
                          !If both A_Plus and A_Minus has gauss points, then:
					      if (Yes_Ap_has_GP .and. Yes_Am_has_GP) then
					          Enriched_Node_Type(NN(i_N),i_C) = 2   ! Heaviside nodes
					      end if
                      end if
			      end do  
              end if
          end do
      end do

      RETURN
      END SUBROUTINE Determine_Enriched_Nodes
