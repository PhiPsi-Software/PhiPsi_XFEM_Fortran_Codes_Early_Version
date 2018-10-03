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
 
      SUBROUTINE Welcome
c     This function display information about the program.


      write(*,"(A60)") "=============================================="
      write(*,"(A60)") "****                                      ****"
      write(*,"(A60)") "****         PhiPsi Version 0.1.0         ****"
      write(*,"(A60)") "****                                      ****"
      write(*,"(A60)") "=============================================="
      write(*,"(A60)") "> Features:                                   "
      write(*,"(A60)") "  PhiPsi is a numerical simulation program   ."
      write(*,"(A60)") "  Author: SHI Fang, University of Science and "
      write(*,"(A60)") "                          Technology of China "
      write(*,"(A60)") "> Release date: July 26, 2015                 "
      write(*,"(A60)") "> Website: http://phipsi.top                  "
      write(*,"(A60)") "> Email: phipsi@sina.cn                       "
      write(*,"(A60)") "=============================================="
      
      RETURN
      END SUBROUTINE Welcome
