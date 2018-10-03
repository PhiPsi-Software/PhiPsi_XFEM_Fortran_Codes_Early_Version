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
 
      SUBROUTINE Input_Check_Display
c     Check input infomation.  
      
      use Global_Common   
      use Global_Filename
      use Global_Model
      use Global_Elem_Area
      use Global_Crack
      
      implicit none
      integer i
      character*2 tem_char
      
 5001 FORMAT('     !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!')  
 

      print *, "    -----------------------------------" 
      print *, "     < Checking input information....>" 


      print *, "    -----------------------------------" 
      print *, '    Filename:          ' // Filename 
      print *, '    Working directory: '  // Work_Dirctory 
     
      
      !***************
      !问题类型
      !***************
      if (Key_Problem .eq. 1)then
          print *, '    Problem type:      plane stress'
      elseif (Key_Problem .eq. 2)then
          print *, '    Problem type:      plane strain'
      else
          call Warning_Message('E','Key_Problem',11) 
      end if

      !***************
      !   material
      !***************
 1001 FORMAT(23X,'Elastic modulus    -- ',F8.3,' GPa')
 1002 FORMAT(23X,'Poisson ratio      -- ',F8.3)
 1003 FORMAT(23X,'Density            -- ',F8.3,' Kg/m^3')
 1004 FORMAT(23X,'Thickness          -- ',F8.3,' m')
 1005 FORMAT(23X,'Tensile strength   -- ',F8.3,' MPa')
 1006 FORMAT(23X,'Fracture toughness -- ',F8.3,' MPa.m^1/2')
      do i=1,size(Material_Type)
          if (Material_Type(i) .eq. 1) then
              WRITE (tem_char,"(i2)") i 
              print *, '    Material type:          '
              print *, '                  *****************************'
     &                                 // '**************'
              print *, '                     Material ' //tem_char//  
     &                                    ' is isotropic material' 
              print *, '                  *****************************'
     &                                 // '**************'
              WRITE(*,1001) Material_Para(i,1)/Con_G  
              WRITE(*,1002) Material_Para(i,2)
              WRITE(*,1003) Material_Para(i,3)
              WRITE(*,1004) Material_Para(i,4)
              WRITE(*,1005) Material_Para(i,5)/Con_M
              WRITE(*,1006) Material_Para(i,6)/Con_M
          end if
      end do
      !*********************
      !Key_Initiation
      !*********************
      if (Key_Initiation .eq. 0)then
          print *, '    Crack initiation:  not allowed'
      elseif (Key_Initiation .eq. 1)then
          print *, '    Crack initiation:  allowed'
      else 
          call Warning_Message('E','Key_Initiation',14)    
      end if  
 

      !**********************************
      !Key_SLOE
      !**********************************
      if     (Key_SLOE .eq. 1)then
          print *, '    SLOE:  solved directly, D = K/F'
      elseif (Key_SLOE .eq. 2)then
          print *, '    SLOE:  solved by Gauss elimination method'
      elseif (Key_SLOE .eq. 3)then
          print *, '    SLOE:  solved by LU factorization method'
      elseif (Key_SLOE .eq. 4)then
          print *, '    SLOE:  solved by iterative solver ITPACK'
      else
          call Warning_Message('E','Key_SLOE',8) 
      end if
      
      !****************************
      !elment and node number
      !****************************
 1011 FORMAT(5X,'Number of elements:   ',I8.3)
 1012 FORMAT(5X,'Number of nodes:      ',I8.3)   
 1013 FORMAT(5X,'Average area of elements:',F12.7,' m^2')   
      WRITE(*,1011) Num_Elem
      WRITE(*,1012) Num_Node
      WRITE(*,1013) Ave_Elem_Area
      
      !*********************
      !range of model
      !*********************
 1021 FORMAT(5X,'Range of x coordinates of the model:   ',F12.7,
     &                 ' m to ',F12.7,' m')  
 1022 FORMAT(5X,'Range of y coordinates of the model:   ',F12.7,
     &                 ' m to ',F12.7,' m')  
      WRITE(*,1021) Min_X_Coor,Max_X_Coor
      WRITE(*,1022) Min_Y_Coor,Max_Y_Coor  
      
      !*********************
      !initial cracks
      !*********************
 1031 FORMAT(5X,'Number of initial cracks:   ',I8.3) 
      WRITE(*,1031) num_Crack
         
      
      print *, "    Pre-Process done!" 
      
      RETURN
      END SUBROUTINE Input_Check_Display