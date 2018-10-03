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
 
      module Global_Common
          implicit none
          double precision pi,Con_M,Con_G
          double precision Delta_Factor,Fac_Pro_Crack
          integer Max_Itera
          parameter(pi    =    3.1415926D0)  !PI
          parameter(Con_M =    1.0D6)        !M,10^6
          parameter(Con_G =    1.0D9)        !G,10^9
          parameter (Delta_Factor    = 0.001D0 )              
          parameter (Max_Itera       = 20000 )            
          integer tip_Order, split_Order,vertex_Order,junction_Order
          integer Num_Gauss_Points,Num_Substeps
          integer Key_Initiation,Key_Propagation,Key_Problem
          double precision Coff_a_AMPT
          integer,ALLOCATABLE::Material_Type(:)             !material type
          double precision,ALLOCATABLE:: Material_Para(:,:) !material parameters
          !integer Key_SIFs_Method,Key_Force_Control,Key_Dynamic
          integer num_Crack_Log(Max_Itera)                  !log of crack
          integer Key_Gravity                               !gravity
          integer Key_SLOE                                  !linear system solver

      end module Global_Common
      
C     -----------------------------      
C     2.global parameters
C     -----------------------------    
      module Global_Model
          implicit none    
          integer Num_Node, Num_Elem, num_of_Material
          integer Num_Bou_x,Num_Bou_y,Num_Foc_x,Num_Foc_y
          double precision Min_X_Coor,Max_X_Coor,Min_Y_Coor,Max_Y_Coor 
          double precision,ALLOCATABLE::Coor(:,:)
          integer,ALLOCATABLE::Bou_x(:)
          integer,ALLOCATABLE::Bou_y(:)
          double precision,ALLOCATABLE::Foc_x(:,:)
          double precision,ALLOCATABLE::Foc_y(:,:)      
          integer,ALLOCATABLE::Elem_Node(:,:)         !nodes of each element
          integer,ALLOCATABLE::Node_Elements(:,:)     !elements of each nodes
          integer,ALLOCATABLE::num_Node_Elements(:)   !number of elements of each nodes
          integer,ALLOCATABLE::Elem_Mat(:)
          integer,ALLOCATABLE::Outline(:,:)
          integer,ALLOCATABLE:: G_NN(:,:)
          double precision,ALLOCATABLE:: G_X_NODES(:,:),G_Y_NODES(:,:) 
          double precision,ALLOCATABLE:: Elem_Area(:)
          integer Total_Freedom
      end module Global_Model

C     -----------------------------      
C     3.filename
C     -----------------------------
      module Global_Filename
          implicit none    
          character(100) Filename,Work_Dirctory
          character(200) Full_Pathname
          character(200) PhiPsi_Current_Directory
      end module Global_Filename     

C     -----------------------------      
C     5.crack
C     -----------------------------
      module Global_Crack
          implicit none 
          integer Max_Num_Crack,Max_Num_Cr_P,
     &    Max_Num_Cr_CalP,Max_Num_Seg_CalP
          parameter (Max_Num_Crack    = 10   )      !max number of cracks
          parameter (Max_Num_Cr_P     = 10   )      !max number of points of each crack
          parameter (Max_Num_Cr_CalP  = 1000 )      !max number of calculation points of each crack
          parameter (Max_Num_Seg_CalP = 200  )      !max number of calculation points of each crack segment
          double precision Crack_Coor(Max_Num_Crack,Max_Num_Cr_P,2) !crack coordinates
          integer num_Crack                                         !number of crack
          integer Each_Cr_Poi_Num(Max_Num_Crack)                    !number of points of each crack
          integer,ALLOCATABLE:: Elem_Type(:,:)
          integer ,ALLOCATABLE:: c_POS(:,:)
          integer,ALLOCATABLE:: Enriched_Node_Type(:,:)
          double precision,ALLOCATABLE::  Coors_Element_Crack(:,:)
          double precision,ALLOCATABLE::  Coors_Tip(:,:)
          double precision,ALLOCATABLE::  Coors_Vertex(:,:)
          double precision,ALLOCATABLE::  Coors_Junction(:,:)
          double precision,ALLOCATABLE::  x_cr_tip_nodes(:,:)
          double precision,ALLOCATABLE::  y_cr_tip_nodes(:,:)
          integer Crack_Tip_Type(Max_Num_Crack,2)
          double precision Edge_Disposed_Crack(Max_Num_Crack,
     &                                         Max_Num_Cr_P,2) 
          logical Flag_Crack_Tip_Out_Mol(Max_Num_Crack,2)    
          integer n_h_Node,n_t_Node,n_j_Node
      end module Global_Crack 
      
C     -----------------------------      
C     6.element area
C     -----------------------------      
      module Global_Elem_Area
          double precision :: Max_Elem_Area,Min_Elem_Area,
     &                        Ave_Elem_Area,Ave_Elem_L
      end module Global_Elem_Area
      
C     -----------------------------      
C     7.material
C     -----------------------------    
      module Global_Material
          double precision,ALLOCATABLE:: D(:,:,:)   
          double precision,ALLOCATABLE:: S(:,:,:)    !inverse matrix of D
          double precision,ALLOCATABLE:: St(:,:)
          double precision,ALLOCATABLE:: Sc(:,:)
          double precision,ALLOCATABLE:: KIc(:,:)
          double precision,ALLOCATABLE:: E(:,:)
          double precision,ALLOCATABLE:: v(:,:)
          double precision,ALLOCATABLE:: density(:)
          double precision,ALLOCATABLE:: thick(:)
      end module Global_Material
C     -----------------------------      
C     8.dispalcement
C     -----------------------------  
      module Global_DISP
          double precision,ALLOCATABLE:: DISP(:)
      end module Global_DISP