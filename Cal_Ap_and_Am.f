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
 
      subroutine Cal_Ap_and_Am(c_Crack,c_Cr_Point,
     &                    DOMAIN_Outline,m_DOMAIN_Outline,
     &                    Domain_El,n_Domain_El,c_E,c_Node,
     &                    Yes_Cutthrough,
     &                    Yes_Ap_has_GP,Yes_Am_has_GP)
c     This function calculates the A+ and A- domain of the support domain of node, used for
c     the selection of the Heaviside node. If either A+ or A- does not contain any gauss point,
c     then the Heaviside node should be removed. Yes_Cutthrough equals to 1 if the support domain
c     is cut through by the crack. 


      use Global_Model
      use Global_Elem_Area
      use Global_Common
      use Global_Crack
      
      implicit none     
      
      integer,intent(in)::c_E,c_Node,c_Cr_Point
      integer,intent(in)::m_DOMAIN_Outline,n_Domain_El
      integer,intent(in)::DOMAIN_Outline(m_DOMAIN_Outline,2)
      integer,intent(in)::Domain_El(n_Domain_El)
      double precision,intent(in)::c_Crack(Max_Num_Cr_P,2)
      
      logical,intent(out)::Yes_Cutthrough,Yes_Ap_has_GP,Yes_Am_has_GP
      !double precision,intent(out):: Area_Ap,Area_Am
      
      integer i_Seg,i_Side,i,j,k,i_E,i_G
      integer tem_Num_Node,tem_count
      double precision Inters_and_Vex_point(10,2),Inters_point(10,2)
      integer new_DOMAIN_Outline(20,2),m_new_DOMAIN_Outline
      integer Uniq_new_DOMAIN_Outline(20,2)
      integer Num_In_and_Vex_Point,Num_Inters_point
      integer Cut_Segments(10,2),m_Cut_Segments
      integer EndByEnd_Cut_Segments(10,2),cou
      integer changed_DOMAIN_Side_Flag(10),num_Cha_DOM_S_Flag
      integer Uniq_changed_DOMAIN_S_F(10)
      double precision  point_seg_1(2),point_seg_2(2)
      double precision  point_side_1(2),point_side_2(2),Ver_point(2)
      double precision  X,Y,Ap_or_Am
      integer Uniq_num_Cha_DOM_S_Flag
      logical Yes_Cross,Yes_point_seg_1_in,Yes_point_seg_2_in
      integer EndByEnd_DOMAIN_Outline(m_DOMAIN_Outline+1) 
      double precision tem_ALL_Coor(Num_Node+10,2)
      logical Yes
      integer Uni_Mat_Count(20),Uniq_m_new_DOMAIN_Outline
      integer num_node_1,num_node_2,num_node_3,Location
      integer PP(3),Uniqued_PP(3),n_Uniqued_PP
      integer Ap_Outline(20,2),Am_Outline(20,2)
      integer Sorted_Ap_Outline(20,2),Sorted_Am_Outline(20,2)
      integer m_Ap_Outline,m_Am_Outline,c_point
      double precision Signed_Distance(2,10),tol,c_point_Coor(2)
      double precision Line_AB(2,2)
      integer cou_1
      double precision,ALLOCATABLE:: Coor_Ap_Closed(:,:),
     &                               Coor_Am_Closed(:,:) 
      double precision kesi(Num_Gauss_Points),
     &                 yita(Num_Gauss_Points),
     &                 weight(Num_Gauss_Points)
      logical In_Ap,In_Am
      integer NODES_iE(4)
      double precision X_NODES(4),Y_NODES(4)
      
      m_Ap_Outline         = 0
      m_Am_Outline         = 0
      Num_In_and_Vex_Point = 0
      Num_Inters_point     = 0
      m_Cut_Segments       = 0
      num_Cha_DOM_S_Flag   = 0
      tem_count            = 0
      tem_Num_Node         = Num_Node
      m_new_DOMAIN_Outline = 0
      num_node_1           = 0
      num_node_2           = 0
      num_node_3           = 0
      cou_1                = 0
      
      In_Ap = .False.
      In_Am = .False.
      
      changed_DOMAIN_Side_Flag(1:10) = 0
      new_DOMAIN_Outline(1:20,1:2)   =0
      Uniq_new_DOMAIN_Outline(1:20,1:2)=0
      Cut_Segments(1:10,1:2)  =0
      EndByEnd_Cut_Segments(1:10,1:2)=0
      Signed_Distance(1:2,1:10)  =0.0D0
      Ap_Outline(1:20,1:2)   =0
      Am_Outline(1:20,1:2)   =0
      Sorted_Ap_Outline(1:20,1:2)   =0
      Sorted_Am_Outline(1:20,1:2)   =0
      Inters_and_Vex_point(1:10,1:2)   =0.0D0
      
      kesi(1:Num_Gauss_Points)=0.0D0
      yita(1:Num_Gauss_Points)=0.0D0
      weight(1:Num_Gauss_Points)=0.0D0
      
      
      Yes_Cutthrough =.False.
      Yes_Ap_has_GP = .False.
      Yes_Am_has_GP = .False.
          
      tem_ALL_Coor(1:Num_Node,:) = Coor

c     Firstly, get the intersection point and kink point of the crack and the edges of the Domain_El.
c     Loop through each segment of the crack. 
      do i_Seg = 1,c_Cr_Point-1
          !Get the coordinate of the two points of the current segment.
          point_seg_1 = c_Crack(i_Seg,:)
          point_seg_2 = c_Crack(i_Seg+1,:)
	      !Loop through each edge of the support domain
          do i_Side = 1,m_DOMAIN_Outline
              point_side_1 = [Coor(DOMAIN_Outline(i_Side,1),1),
     &                        Coor(DOMAIN_Outline(i_Side,1),2)]
              point_side_2 = [Coor(DOMAIN_Outline(i_Side,2),1),
     &                        Coor(DOMAIN_Outline(i_Side,2),2)]    
     	      !Check if they intersects and calculate intersection point.
              call Tool_Cal_Intersection(
     &            point_seg_1,point_seg_2,
     &            point_side_1,point_side_2,
     &            X,Y,Yes_Cross) 
              if(Yes_Cross.eqv..True.)then
                  !Intersection points.
                  tem_count = tem_count + 1
                  Num_In_and_Vex_Point = Num_In_and_Vex_Point +1
                  !print *,Num_In_and_Vex_Point
                  Inters_and_Vex_point(Num_In_and_Vex_Point,:) = [X,Y]
                  Num_Inters_point = Num_Inters_point +1
                  Inters_point(Num_Inters_point,:) = [X,Y]
                  tem_Num_Node = tem_Num_Node+1
                  tem_ALL_Coor(tem_Num_Node,:) = [X,Y]
                  !print *,tem_Num_Node
                  m_new_DOMAIN_Outline = m_new_DOMAIN_Outline +1
			      new_DOMAIN_Outline(m_new_DOMAIN_Outline,:) =
     &                    [DOMAIN_Outline(i_Side,1),tem_Num_Node]
                  m_new_DOMAIN_Outline = m_new_DOMAIN_Outline +1
			      new_DOMAIN_Outline(m_new_DOMAIN_Outline,:) =
     &                    [DOMAIN_Outline(i_Side,2),tem_Num_Node]      
                  num_Cha_DOM_S_Flag = num_Cha_DOM_S_Flag +1 
                  changed_DOMAIN_Side_Flag(num_Cha_DOM_S_Flag)=i_Side                 
              end if
              
          end do
          !If there are intersection points, then:
          if (tem_count.ne.0) then
              Yes_point_seg_1_in = .False.
              Yes_point_seg_2_in = .False.
              if (i_Seg.ne. 1) then
                  !change Outline to EndByEnd_DOMAIN_Outline which is end by end
                  EndByEnd_DOMAIN_Outline(1:m_DOMAIN_Outline) = 
     &                                 DOMAIN_Outline(:,1)
                  EndByEnd_DOMAIN_Outline(m_DOMAIN_Outline+1) = 
     &                                 DOMAIN_Outline(1,1)
                  !check if point inside EndByEnd_DOMAIN_Outline
                  Call Tool_Yes_In_Poly
     &                 (point_seg_1(1),point_seg_1(2),
     &                  Coor(EndByEnd_DOMAIN_Outline,1),
     &                  Coor(EndByEnd_DOMAIN_Outline,2),
     &                  size(EndByEnd_DOMAIN_Outline),
     &                  Yes_point_seg_1_in)
              end if
              if (i_Seg .ne. (c_Cr_Point-1)) then
                  EndByEnd_DOMAIN_Outline(1:m_DOMAIN_Outline) = 
     &                                 DOMAIN_Outline(:,1)
                  EndByEnd_DOMAIN_Outline(m_DOMAIN_Outline+1) = 
     &                                 DOMAIN_Outline(1,1)
                  Call Tool_Yes_In_Poly
     &                 (point_seg_2(1),point_seg_2(2),
     &                  Coor(EndByEnd_DOMAIN_Outline,1),
     &                  Coor(EndByEnd_DOMAIN_Outline,2),
     &                  size(EndByEnd_DOMAIN_Outline),
     &                  Yes_point_seg_2_in)                   
              end if
              if(Yes_point_seg_1_in)then  
                  call Vector_belongs_Matrix_Is_Dou(10,2,
     &                       Inters_and_Vex_point,
     &                       [point_seg_1(1),point_seg_1(2)],
     &                       Location,Yes)
                  if    (Yes.eqv..False.) then
                      Num_In_and_Vex_Point = Num_In_and_Vex_Point +1
                      Inters_and_Vex_point(Num_In_and_Vex_Point,:) = 
     &                                   [point_seg_1(1),point_seg_1(2)]
                      Ver_point = [point_seg_1(1),point_seg_1(2)]
                  end if
              end if
              if(Yes_point_seg_2_in)then 
                  call Vector_belongs_Matrix_Is_Dou(10,2,
     &                       Inters_and_Vex_point,
     &                       [point_seg_2(1),point_seg_2(2)],
     &                       Location,Yes)
                  if(Yes .eqv..False.) then
                      Num_In_and_Vex_Point = Num_In_and_Vex_Point +1
                      Inters_and_Vex_point(Num_In_and_Vex_Point,:) = 
     &                                   [point_seg_2(1),point_seg_2(2)]
                      Ver_point = [point_seg_2(1),point_seg_2(2)]
                  end if                  
              end if
          end if
      end do
      
c     If has no or has only one intersection point, then end the function.
c     If has 3 intersection point, the Domain_El may be a Concave Polygon, it is more complex, 
c     so here we give up and end the function.
      if ((Num_In_and_Vex_Point.eq.0)  .or.
     &     (Num_Inters_point  .eq. 1)  .or.
     &     (Num_Inters_point  .ge. 3)) then
	      Yes_Cutthrough =.False.
	      Yes_Ap_has_GP = .True.
	      Yes_Am_has_GP = .True.
	      return
      end if
     
      !delete reduplicated content in changed_DOMAIN_Side_Flag
      call Vector_Unique_Int(10,
     &                     num_Cha_DOM_S_Flag,changed_DOMAIN_Side_Flag,
     &                     Uniq_changed_DOMAIN_S_F,
     &                     Uniq_num_Cha_DOM_S_Flag)
      
      do i_Side =1,m_DOMAIN_Outline
          if (any(i_Side .eq. changed_DOMAIN_Side_Flag
     &                  (1:Uniq_num_Cha_DOM_S_Flag)).eqv..False.) then
              m_new_DOMAIN_Outline = m_new_DOMAIN_Outline +1
              new_DOMAIN_Outline(m_new_DOMAIN_Outline,:) = 
     &                              DOMAIN_Outline(i_Side,:)
	     end if
      end do
      
      !delete reduplicated content
      call Matrix_Unique_Row_Int(20,2,m_new_DOMAIN_Outline,
     &                             new_DOMAIN_Outline,
     &            Uniq_new_DOMAIN_Outline,Uniq_m_new_DOMAIN_Outline,
     &                             Uni_Mat_Count)    
 
      !Get Cut_Segments for different case.
      select case(Num_In_and_Vex_Point)
      !Case 1:
      case(2)
          call Vector_belongs_Matrix_Is_Dou(Num_Node+10,2,tem_ALL_Coor,
     &                              Inters_and_Vex_point(1,:),
     &                                    num_node_1,Yes) 
          call Vector_belongs_Matrix_Is_Dou(Num_Node+10,2,tem_ALL_Coor,
     &                              Inters_and_Vex_point(2,:),
     &                                    num_node_2,Yes)   
          m_Cut_Segments = m_Cut_Segments +1
  	      Cut_Segments(m_Cut_Segments,:) = [num_node_1,num_node_2]
      !Case 2:	
      case(3)
 	      tem_Num_Node = tem_Num_Node + 1
          tem_ALL_Coor(tem_Num_Node,:) = Ver_point
          call Vector_belongs_Matrix_Is_Dou(Num_Node+10,2,tem_ALL_Coor,
     &                              Inters_and_Vex_point(1,:),
     &                                    num_node_1,Yes) 
          call Vector_belongs_Matrix_Is_Dou(Num_Node+10,2,tem_ALL_Coor,
     &                              Inters_and_Vex_point(2,:),
     &                                    num_node_2,Yes) 
          call Vector_belongs_Matrix_Is_Dou(Num_Node+10,2,tem_ALL_Coor,
     &                              Inters_and_Vex_point(3,:),
     &                                    num_node_3,Yes) 
 	      PP = [num_node_1,num_node_2,num_node_3]
          call Vector_Unique_Int(3,3,PP,
     &                             Uniqued_PP,n_Uniqued_PP) 
          call Vector_Sort_Int(n_Uniqued_PP,Uniqued_PP(1:n_Uniqued_PP))  
          m_Cut_Segments = m_Cut_Segments +1
          Cut_Segments(m_Cut_Segments,:) = [Uniqued_PP(1),tem_Num_Node]
          m_Cut_Segments = m_Cut_Segments +1
          Cut_Segments(m_Cut_Segments,:) = [Uniqued_PP(2),tem_Num_Node]
      !Case 3:	
      !In this case, obviously, both Ap and Am has gauss point, because 
      !a full crack increment delta_L is inside the Domain_El.	
      case(4)
	      Yes_Cutthrough= .True.
	      Yes_Ap_has_GP =.True.
	      Yes_Am_has_GP =.True.
      !Case 4:	
      case default
	     Yes_Cutthrough=.False.
	     Yes_Ap_has_GP=.False.
	     Yes_Am_has_GP=.False.
      end select
      
      if(((Num_In_and_Vex_Point.eq.2) .or. (Num_In_and_Vex_Point.eq.3))
     &   .and. (m_Cut_Segments.ne.0))     then
          Yes_Cutthrough = .True.
          call Tool_Sort_by_End_to_End(10,m_Cut_Segments,Cut_Segments,
     &                                      EndByEnd_Cut_Segments,cou)
          tol = sqrt(Ave_Elem_Area)*1.0D-8
          do i=1,m_new_DOMAIN_Outline
		      !Loop through two points of the new_DOMAIN_Outline.
		      !Signed_Distance = zeros(2,size(Inters_point,1)-1)
              do j=1,2
                  c_point = Uniq_new_DOMAIN_Outline(i,j)
	              c_point_Coor = tem_ALL_Coor(c_point,:)
			      do k=1,Num_Inters_point-1
                      point_seg_1 = Inters_point(k,:)
	                  point_seg_2 = Inters_point(k+1,:)
                      Line_AB(1,:)= point_seg_1
                      Line_AB(2,:)= point_seg_2 
	                  call Cal_Signed_Distance
     &                          (Line_AB,c_point_Coor,
     &                           Signed_Distance(j,k))
 	                  if (abs(Signed_Distance(j,k)) .le. tol) then
 					      Signed_Distance(j,k) = 0.0D0
 	                  end if
	              end do
              end do
              Ap_or_Am = maxval(Signed_Distance(:,1:Num_Inters_point-1))
		      if (Ap_or_Am .gt. tol) then
                  m_Ap_Outline = m_Ap_Outline + 1
		          Ap_Outline(m_Ap_Outline,:)=Uniq_new_DOMAIN_Outline(i,:)
              else
                  m_Am_Outline = m_Am_Outline + 1
		          Am_Outline(m_Am_Outline,:)=Uniq_new_DOMAIN_Outline(i,:)
		      end if
	      end do
          
          do i=1,m_Cut_Segments
              m_Ap_Outline = m_Ap_Outline + 1
              Ap_Outline(m_Ap_Outline,:) = EndByEnd_Cut_Segments(i,:)
              m_Am_Outline = m_Am_Outline + 1
              Am_Outline(m_Am_Outline,:) = EndByEnd_Cut_Segments(i,:)
	      end do  
          
          call Tool_Sort_by_End_to_End(20,m_Ap_Outline,Ap_Outline,
     &                 Sorted_Ap_Outline,cou_1)
          call Tool_Sort_by_End_to_End(20,m_Am_Outline,Am_Outline,
     &                 Sorted_Am_Outline,cou_1)
           allocate(Coor_Ap_Closed(m_Ap_Outline+1,2))
           allocate(Coor_Am_Closed(m_Am_Outline+1,2))
          Coor_Ap_Closed(1:m_Ap_Outline,:) = 
     &               tem_ALL_Coor(Sorted_Ap_Outline(1:m_Ap_Outline,1),:)
          Coor_Ap_Closed(m_Ap_Outline+1,:) = 
     &               tem_ALL_Coor(Sorted_Ap_Outline(1,1),:)
     
          Coor_Am_Closed(1:m_Am_Outline,:) = 
     &               tem_ALL_Coor(Sorted_Am_Outline(1:m_Am_Outline,1),:)   
          Coor_Am_Closed(m_Am_Outline+1,:) = 
     &               tem_ALL_Coor(Sorted_Am_Outline(1,1),:) 

	      !Check if Ap has gauss points.
          call Cal_Gauss_Points_QUAD(Num_Gauss_Points,
     &                                   kesi,yita,weight)
	      do i_E = 1,n_Domain_El
              NODES_iE = [Elem_Node(Domain_El(i_E),1),
     &                    Elem_Node(Domain_El(i_E),2),
     &                    Elem_Node(Domain_El(i_E),3),
     &                    Elem_Node(Domain_El(i_E),4)]             
	          X_NODES = Coor(NODES_iE,1)                           
	          Y_NODES = Coor(NODES_iE,2)
	          do i_G = 1,Num_Gauss_Points
                  call Cal_Coor_by_KesiYita(kesi(i_G),yita(i_G),
     &                                      X_NODES,Y_NODES,
     &                                      x,y)
                  call Tool_Yes_In_Poly(x,y,
     &                        Coor_Ap_Closed(:,1),
     &                        Coor_Ap_Closed(:,2),
     &                        m_Ap_Outline+1,In_Ap)
			      if (In_Ap) then
                     Yes_Ap_has_GP = .True.
                     goto 201
	              end if
	          end do
          end do        
  201     continue
	      !Check if Ap has gauss points.
	      do i_E = 1,n_Domain_El
              NODES_iE = [Elem_Node(Domain_El(i_E),1),
     &                    Elem_Node(Domain_El(i_E),2),
     &                    Elem_Node(Domain_El(i_E),3),
     &                    Elem_Node(Domain_El(i_E),4)]             
	          X_NODES = Coor(NODES_iE,1)                           
	          Y_NODES = Coor(NODES_iE,2)
	          do i_G = 1,Num_Gauss_Points
                  call Cal_Coor_by_KesiYita(kesi(i_G),yita(i_G),
     &                                      X_NODES,Y_NODES,
     &                                      x,y)
                  call Tool_Yes_In_Poly(x,y,
     &                        Coor_Am_Closed(1:,1),
     &                        Coor_Am_Closed(:,2),
     &                        m_Am_Outline+1,In_Am)
			      if (In_Am) then
                     Yes_Am_has_GP = .True.
                     goto 202
	              end if
	          end do
          end do        
  202     continue          
      end if
      
      return
      end SUBROUTINE Cal_Ap_and_Am                         
