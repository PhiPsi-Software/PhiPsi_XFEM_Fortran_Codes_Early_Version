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
 
      SUBROUTINE Read_Geo
      !read model info    
          
      use Global_Filename
      use Global_Model
      use Global_Elem_Area
      use Global_Crack

      implicit none
      LOGICAL ALIVE
      character*200 temp_name
      integer Tool_Count_Lines
      double precision,ALLOCATABLE::Temp_DATA(:,:)
      logical Flag_Blank   
      integer,ALLOCATABLE::tem1(:,:)
      integer,ALLOCATABLE::tem2(:)
      integer i,j,all_num_Outline,Out_num_Outline
      integer,ALLOCATABLE::All_Outline(:,:)   
      integer,ALLOCATABLE::Temp_Outline(:,:)  
      integer N1,N2,N3,N4,NN(4)
      double precision area
      
 1001 FORMAT('     !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!')
 
c     -----------------------      
c     read *.node file
c     -----------------------
      print *, "    Trying to read nodal files...." 
      temp_name = trim(trim(Full_Pathname)//'.node')
      inquire(file=temp_name, exist=alive)  !如果存在，alive赋值为True
      if(alive.EQV..FALSE.)then
	      WRITE(*,1001)  
          print *, "    Error :: Can not find nodal files!!!" 
          call Warning_Message('S',' ',0) 
      else
          !count number of row of node file
          Num_Node = Tool_Count_Lines(temp_name) !get node number
          ALLOCATE(Coor(Num_Node,2))
          ALLOCATE( Temp_DATA(Num_Node,2))
          Call Tool_Read_File(temp_name,"node",Num_Node,2,Temp_DATA,
     &                        Flag_Blank)
          Coor = Temp_DATA
          DEALLOCATE(Temp_DATA)
          !range of model
          Max_X_Coor = maxval(Coor(:,1))
          Min_X_Coor = minval(Coor(:,1))
          Max_Y_Coor = maxval(Coor(:,2))
          Min_Y_Coor = minval(Coor(:,2))
      endif
      
c     -----------------------     
c     read *.elem file
c     -----------------------
      print *, "    Trying to read element files...." 
      temp_name = trim(trim(Full_Pathname)//'.elem')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
          WRITE(*,1001)  
          print *, "    Error :: Can not find element files!!!" 
          call Warning_Message('S',' ',0) 
      else
          Num_Elem = Tool_Count_Lines(temp_name) !element number
          ALLOCATE( Elem_Node(Num_Elem,4))
          ALLOCATE( Elem_Mat(Num_Elem)) 
          ALLOCATE( Temp_DATA(Num_Elem,5))
          Call Tool_Read_File(temp_name,"elem",Num_Elem,5,Temp_DATA,
     &                        Flag_Blank)
          Elem_Node = int(Temp_DATA(:,1:4))
          Elem_Mat  = int(Temp_DATA(:,5))
          num_of_Material = MaxVal(Elem_Mat)     !material munber of the element
          DEALLOCATE(Temp_DATA)
      endif  
      
c     ------------------------------      
c     read *.boux file
c     ------------------------------
      print *, "    Trying to read boux files...." 
      temp_name = trim(trim(Full_Pathname)//'.boux')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
          print *, "    Warning :: Can not find boux files!!!" 
      else
          Num_Bou_x = Tool_Count_Lines(temp_name) 
          ALLOCATE( Bou_x(Num_Bou_x))
          ALLOCATE( Temp_DATA(Num_Bou_x,1))
          Call Tool_Read_File(temp_name,"boux",Num_Bou_x,1,Temp_DATA,
     &                        Flag_Blank)
          Bou_x  = int(Temp_DATA(:,1))
          DEALLOCATE(Temp_DATA)
      endif  
      
c     ------------------------------      
c     read *.bouy file
c     ------------------------------
      print *, "    Trying to read bouy files...." 
      temp_name = trim(trim(Full_Pathname)//'.bouy')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
          print *, "    Warning :: Can not find bouy files!!!" 
      else
          Num_Bou_y = Tool_Count_Lines(temp_name) 
          ALLOCATE( Bou_y(Num_Bou_y))
          ALLOCATE( Temp_DATA(Num_Bou_y,1))
          Call Tool_Read_File(temp_name,"bouy",Num_Bou_y,1,Temp_DATA,
     &                        Flag_Blank)
          Bou_y  = int(Temp_DATA(:,1))
          DEALLOCATE(Temp_DATA)
      endif    
      
c     ------------------------------      
c     read *.focx file
c     ------------------------------
      print *, "    Trying to read focx files...." 
      temp_name = trim(trim(Full_Pathname)//'.focx')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
          print *, "    Warning :: Can not find focx files!!!" 
      else
          Num_Foc_x = Tool_Count_Lines(temp_name) 
          ALLOCATE( Foc_x(Num_Foc_x,2))
          ALLOCATE( Temp_DATA(Num_Foc_x,2))
          Call Tool_Read_File(temp_name,"focx",Num_Foc_x,2,Temp_DATA,
     &                        Flag_Blank)
          Foc_x  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
      endif  
      
c     ------------------------------      
c     read *.focy file
c     ------------------------------
      print *, "    Trying to read focy files...." 
      temp_name = trim(trim(Full_Pathname)//'.focy')
      inquire(file=temp_name, exist=alive)  
      if(alive.EQV..FALSE.)then
          print *, "    Warning :: Can not find focy files!!!" 
      else
          Num_Foc_y = Tool_Count_Lines(temp_name) 
          ALLOCATE( Foc_y(Num_Foc_y,2))
          ALLOCATE( Temp_DATA(Num_Foc_y,2))
          Call Tool_Read_File(temp_name,"focy",Num_Foc_y,2,Temp_DATA,
     &                        Flag_Blank)
          Foc_y  = Temp_DATA(:,:)
          DEALLOCATE(Temp_DATA)
      endif 
      
c     ---------------------------------------------------------      
c     get the boundary of the model, and store in Outline
c     ---------------------------------------------------------
      ALLOCATE(tem1(4*Num_Elem,2))
      ALLOCATE(tem2(4*Num_Elem))
      tem1(1:Num_Elem,1)              = Elem_Node(1:Num_Elem,1)
      tem1(1:Num_Elem,2)              = Elem_Node(1:Num_Elem,2)
      tem1(Num_Elem+1:2*Num_Elem,1)   = Elem_Node(1:Num_Elem,2)
      tem1(Num_Elem+1:2*Num_Elem,2)   = Elem_Node(1:Num_Elem,3)
      tem1(2*Num_Elem+1:3*Num_Elem,1) = Elem_Node(1:Num_Elem,3)
      tem1(2*Num_Elem+1:3*Num_Elem,2) = Elem_Node(1:Num_Elem,4)
      tem1(3*Num_Elem+1:4*Num_Elem,1) = Elem_Node(1:Num_Elem,4)
      tem1(3*Num_Elem+1:4*Num_Elem,2) = Elem_Node(1:Num_Elem,1)
      call Matrix_Sort_Int(4*Num_Elem,2,tem1)
      call Matrix_Count_Row_Int(4*Num_Elem,2,tem1,tem2,all_num_Outline)
      ALLOCATE(all_Outline(all_num_Outline,2))
      ALLOCATE(Temp_Outline(all_num_Outline,2))
      j=0
      do i=1,4*Num_Elem
          if (tem2(i).eq.1)then
              j=j+1
              all_Outline(j,:) = tem1(i,:)
          end if
      end do
      call Tool_Sort_by_End_to_End(all_num_Outline,all_num_Outline,
     &                              All_Outline,Temp_Outline,
     &                              Out_num_Outline)
      ALLOCATE(Outline(Out_num_Outline,2))
      do i=1,Out_num_Outline
          Outline(i,:) = Temp_Outline(i,:) 
      end do
      
c     -----------------------------------------------
c     get coordinates of elements
c     -----------------------------------------------
      ALLOCATE(G_NN(Num_Elem,4))     
      ALLOCATE(G_X_NODES(Num_Elem,4))
      ALLOCATE(G_Y_NODES(Num_Elem,4))
      ALLOCATE(Elem_Area(Num_Elem))
      do i=1,Num_Elem
          N1  = Elem_Node(i,1)                                                
          N2  = Elem_Node(i,2)                                              
          N3  = Elem_Node(i,3)                                             
          N4  = Elem_Node(i,4)                                                
	      NN  = [N1,N2,N3,N4]                                                
          G_NN(i,:)      = NN
          G_X_NODES(i,:) = Coor(NN,1)
          G_Y_NODES(i,:) = Coor(NN,2) 
          call Tool_Area_Polygon(Coor(NN,1),Coor(NN,2),4,
     &                           area)
          Elem_Area(i) = area
      end do
      Max_Elem_Area = MaxVal(Elem_Area)
      Min_Elem_Area = MinVal(Elem_Area)
      Ave_Elem_Area = sum(Elem_Area)/size(Elem_Area)
      Ave_Elem_L    = sqrt(Ave_Elem_Area)
      
c     -----------------------------------------------
c     get elements of each node
c     -----------------------------------------------
      ALLOCATE(Node_Elements(Num_Node,8))     
      ALLOCATE(num_Node_Elements(Num_Node))   
      do i=1,Num_Node
          num_Node_Elements(i) = 0
          do j=1,Num_Elem
              if (any(Elem_Node(j,:) == i))then
                  num_Node_Elements(i) = num_Node_Elements(i) +1
                  Node_Elements(i,num_Node_Elements(i)) = j
              end if
	      end do
      end do      
      
c     ----------------------------
c     allocate temporary data
c     ----------------------------
      ALLOCATE(Elem_Type(Num_Elem,Max_Num_Crack))
      ALLOCATE(Enriched_Node_Type(Num_Node,Max_Num_Crack))
      ALLOCATE(Coors_Element_Crack(Num_Elem,4))
      ALLOCATE(Coors_Tip(Num_Elem,2))
      ALLOCATE(Coors_Vertex(Num_Elem,2))
      ALLOCATE(Coors_Junction(Num_Elem,4))
      ALLOCATE(x_cr_tip_nodes(Max_Num_Crack,Num_Node))
      ALLOCATE(y_cr_tip_nodes(Max_Num_Crack,Num_Node))
      ALLOCATE(c_POS(Num_Node,Max_Num_Crack))
      
      RETURN
      END SUBROUTINE Read_Geo