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
 
      SUBROUTINE PhiPsi2D_Static
          
      use Global_Common   
      use Global_Filename
      use Global_Model
      use Global_Elem_Area
      use Global_Crack
      use Global_DISP
      
      implicit none
      integer isub
      logical Yes_Last_Growth  
      double precision Lambda
      integer,ALLOCATABLE::freeDOF(:)
      double precision,ALLOCATABLE:: globalK(:,:),globalF(:)
      double precision,ALLOCATABLE::tem_DISP(:)
      integer num_freeDOF
      double precision Max_Disp_x,Max_Disp_y,Min_Disp_x,Min_Disp_y
      
      Yes_Last_Growth = .False.
      
 1001 FORMAT(' >>  Iteration ',I5,' of ',I5,' started:')   
 1002 FORMAT('     Force factor is ',F5.3)   
 1021 FORMAT(5X,'Range of displacement x:   ',F12.8,
     &                 ' m to ',F12.8,' m')  
 1022 FORMAT(5X,'Range of displacement y:   ',F12.8,
     &                 ' m to ',F12.8,' m')  
c     ------------------------------        
c     get material info   
c     ------------------------------ 
      call Get_Material_Matrix
      
c     ------------------------------        
c     Loop over load steps   
c     ------------------------------  
      print *, "  " 
      do isub = 1,Num_Substeps
          !get force factor
          call Force_Factor(Lambda,isub,Yes_Last_Growth)
          WRITE(*,1001) isub,Num_Substeps
          WRITE(*,1002) Lambda
C         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C         %%%%%%%%    if crack exists (XFEM)  %%%%%%%%%%
C         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if(num_Crack_Log(isub).ne.0)then
              !determine enriched nodes   
              call Determine_Enriched_Nodes(isub)
              
              !number enriched nodes 
              call Number_Enriched_Nodes(isub)     
              print *,'    Total_Freedom:',Total_Freedom        
              
              !save crack file
              call Save_Files_Crack(isub)
              
              !assemble stiffness matrix 
              ALLOCATE(globalK(Total_Freedom,Total_Freedom))
              call Assemble_Stiffness_Matrix_XFEM(isub,globalK)
              print *,'    Sum of globalK: ',  sum(globalK)        

              !assemble force vector
              ALLOCATE(globalF(Total_Freedom))
              call Force_Vector(Total_Freedom,isub,Lambda,globalF)

              !deal with boundary condition
              ALLOCATE(freeDOF(Total_Freedom))
              call Boundary_Cond(Total_Freedom,isub,
     &                         freeDOF,num_freeDOF)
              
              !solve KU=F to get displacements of nodes
              ALLOCATE(DISP(Total_Freedom))
              DISP(1:Total_Freedom) = 0.0D0
              ALLOCATE(tem_DISP(num_freeDOF))
              print *,'    Solving displacements......'   
              call Matrix_Solve_LSOE(Key_SLOE,
     &                     globalK(freeDOF(1:num_freeDOF),
     &                     freeDOF(1:num_freeDOF)),
     &                     globalF(freeDOF(1:num_freeDOF)),
     &                     tem_DISP,
     &                     num_freeDOF)
              DISP(freeDOF(1:num_freeDOF)) = tem_DISP
	          Max_Disp_x = maxval(DISP(1:2*Num_Node:2))
              Min_Disp_x = minval(DISP(1:2*Num_Node:2))
	          Max_Disp_y = maxval(DISP(2:2*Num_Node:2))
              Min_Disp_y = minval(DISP(2:2*Num_Node:2))
              WRITE(*,1021) Min_Disp_x,Max_Disp_x
              WRITE(*,1022) Min_Disp_y,Max_Disp_y
              
              !save displacement file
              call Save_Disp(isub)      

              !update crack number
              num_Crack              = 1
              num_Crack_Log(isub+1)  = num_Crack
              
              !clear temporary data
              deallocate(globalK)
              deallocate(globalF)
              deallocate(freeDOF)
              deallocate(DISP)
              deallocate(tem_DISP)
C         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C         %%%%%%%%       if no crack (FEM)       %%%%%%%
C         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          elseif (num_Crack_Log(isub).eq.0)then 
              !total number of defs
              Total_Freedom = 2*Num_Node
              print *,'    Total_Freedom:',Total_Freedom  
              
              !assemble stiffness matrix 
              ALLOCATE(globalK(Total_Freedom,Total_Freedom))
              call Assemble_Stiffness_Matrix_FEM(isub,globalK)    
              print *,'    Sum of globalK: ',  sum(globalK)
              
              !assemble force vector
              ALLOCATE(globalF(Total_Freedom))
              call Force_Vector(Total_Freedom,isub,Lambda,globalF)
              
              !deal with boundary condition
              ALLOCATE(freeDOF(Total_Freedom))
              call Boundary_Cond(Total_Freedom,isub,
     &                         freeDOF,num_freeDOF)    
     
              !solve KU=F to get displacements of nodes
              ALLOCATE(DISP(Total_Freedom))
              DISP(1:Total_Freedom) = 0.0D0
              ALLOCATE(tem_DISP(num_freeDOF))
              print *,'    Solving displacements......'   
              call Matrix_Solve_LSOE(Key_SLOE,
     &                     globalK(freeDOF(1:num_freeDOF),
     &                     freeDOF(1:num_freeDOF)),
     &                     globalF(freeDOF(1:num_freeDOF)),
     &                     tem_DISP,
     &                     num_freeDOF)
              DISP(freeDOF(1:num_freeDOF)) = tem_DISP
	          Max_Disp_x = maxval(DISP(1:Total_Freedom:2))
              Min_Disp_x = minval(DISP(1:Total_Freedom:2))
	          Max_Disp_y = maxval(DISP(2:Total_Freedom:2))
              Min_Disp_y = minval(DISP(2:Total_Freedom:2))
              WRITE(*,1021) Min_Disp_x,Max_Disp_x
              WRITE(*,1022) Min_Disp_y,Max_Disp_y
              
              !save displacement file
              call Save_Disp(isub)     

              
              !clear temporary data           
              deallocate(globalK)
              deallocate(globalF)
              deallocate(freeDOF)
              deallocate(DISP)
              deallocate(tem_DISP)
          end if
      end do
      

 
      RETURN
      END SUBROUTINE PhiPsi2D_Static
