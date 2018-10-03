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
 
      SUBROUTINE PhiPsi2D_Input
            
      use Global_Common
      use Global_Filename
      use Global_Crack

      implicit none
      
      print *, " >> Pre-Processing...." 
      print *, "    Reading input control file...." 
      
c     Filename and Work_Dirctory

      Filename      =  'tension_test_data'                                                                     
      Work_Dirctory =  trim(PhiPsi_Current_Directory)  !PhiPsi_Current_Directory is current work directory
     &                  //'\'//'tension_test_data'                 

      Full_Pathname =  trim(Work_Dirctory)//'\'//trim(Filename) 
      
      
C     define initial crack
      num_Crack = 1                        
      num_Crack_Log(1)    = num_Crack
      Each_Cr_Poi_Num(1)  = 3               !number of points of crack

c     initial crack  
      Crack_Coor(1,1,1:2) = [-2.84D0,3.5D0]  
      Crack_Coor(1,2,1:2) = [ 1.41D0,3.5D0]   
      Crack_Coor(1,3,1:2) = [ 2.41D0,3.5D0]

c     contral parameters
      tip_Order        = 7
      split_Order      = 4   
      vertex_Order     = 4   
      junction_Order   = 4 
      Key_SLOE         = 4       !Linear solver: itpack
      Num_Gauss_Points = 64      !number of gauss points
      
      
      Num_Substeps     = 1       !number of substeps
      Key_Problem      = 2       !problem type: =1,plane stress; =2,plane strain
      
      
c     material parameters
      ALLOCATE(Material_Type(1))
      ALLOCATE(Material_Para(1,15))
      Material_Type    = [1]   
      Material_Para(1,1:15) = 
     &       [30.0D9,0.2D0,2000.0D0,20.0D-3,-100.0D6,
     &        -100D6,  5D6,   1.0D0,  0.0D0,   0.0D0,
     &         0.0D0,0.0D0,   0.0D0,  0.0D0,   0.0D0]
      
      RETURN
      END SUBROUTINE PhiPsi2D_Input