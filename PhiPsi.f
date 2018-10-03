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
 
      PROGRAM PhiPsi2D
          
      use Global_Common
      use omp_lib
      use Global_Filename
      
      integer date_time(8)
      integer S_time,F_time
      character*10  current_data  
      
 1000 FORMAT('  ') 
 1001 FORMAT('     Start time£º',A8,', ',I2,':',I2,':',I2)   
 1002 FORMAT('     End time£º',A8,', ',I2,':',I2,':',I2)   
 1003 FORMAT('     Total elapsed time£º',I9,' ms, about ',F7.3,' mins')   
 
c     Get current work directory of PhiPsi
      CALL getcwd(PhiPsi_Current_Directory) 
      
c     show welcome infomation
      CALL Welcome
      
c     time initialization
      call Tool_Get_Current_Time(current_data,date_time,S_time)
      WRITE(*,1000)
      WRITE(*,1001) current_data,date_time(5),date_time(6),date_time(7)
      WRITE(*,1000)
      
c     read input info    
      CALL PhiPsi2D_Input

c     read model data   
      CALL Read_Geo
      
c     check input data
      CALL Input_Check_Display      
      
c     perform the calculation
      CALL PhiPsi2D_Static

c     get current time
      call Tool_Get_Current_Time(current_data,date_time,F_time)
      print *,' '
      WRITE(*,1002) current_data,date_time(5),date_time(6),date_time(7)
      WRITE(*,1003) F_time-S_time,(F_time-S_time)/1000.0D0/60.0D0

C     close
      print *,' '
      write( *, * ) '    Press any key to exist PhiPsi.' 
      read( *, * )    
    
      !pause
      RETURN
      END PROGRAM PhiPsi2D