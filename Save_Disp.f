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
 
      SUBROUTINE Save_Disp(isub)
c     save displacement.
      
      use Global_Common
      use Global_Model
      use Global_Filename
      use Global_DISP
      
      implicit none
      
      integer isub
      integer i
      character(200) c_File_name_1
      character(5) temp
      
      print *,'    Saving disp......'
      write(temp,'(I4)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &                 //'disp'//'_'//ADJUSTL(temp)   
    
      open(101,file=c_File_name_1,status='unknown')     
      do i=1,Total_Freedom
          write(101, '(1E20.12)') DISP(i)
      end do
      close(101)          
      
      write(temp,'(I4)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'disn'//'_'//ADJUSTL(temp)    
      open(101,file=c_File_name_1,status='unknown')     
      do i=1,Total_Freedom
          write(101, '(I8,2E20.12)') i,DISP(2*i-1),DISP(2*i)
      end do
      close(101)  

 
      RETURN
      END SUBROUTINE Save_Disp
