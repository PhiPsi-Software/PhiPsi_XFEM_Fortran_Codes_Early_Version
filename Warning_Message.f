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
 
      SUBROUTINE Warning_Message(Mess_Type,Keywords,size_Keyword)
c     Show error and warning message.
      
      integer size_Keyword
      character*1  Mess_Type
      character*20 Keywords
      
      character*100 Mess_Content
      
 1001 FORMAT('     ?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?-?')      
      select case(Mess_Type)
          
      case('E')   ! Error
          Mess_Content =  '       '//'ERROR :: Keyword *'
     &                    //trim(Keywords(1:size_Keyword))//' Wrong!'
          WRITE(*,1001)  
          print *,trim(Mess_Content)
          WRITE(*,1001)   
          write( *, * ) '    Close the Program and Check Input Files' 
     &                      //' or Press Enter to Continue.' 
          read( *, * )        
          
      case('W')   ! Warning
          Mess_Content =  '       '//'Warning :: Keyword *'
     &                    //trim(Keywords(1:size_Keyword))//' Wrong!'      
          WRITE(*,1001)  
          print *,trim(Mess_Content)
          WRITE(*,1001)  
          
      case('S')   ! Error and just stop
          Mess_Content =  '       '//'An error occurred.'  
          WRITE(*,1001)  
          print *,trim(Mess_Content)
          WRITE(*,1001)   
          write( *, * ) '    Close the Program and Check Input Files' 
     &                      //' or Press Enter to Continue.' 
          read( *, * )               
      end select
      

      
      RETURN
      END SUBROUTINE Warning_Message
