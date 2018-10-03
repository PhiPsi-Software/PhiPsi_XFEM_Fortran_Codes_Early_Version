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
 
      SUBROUTINE Save_Files_Crack(isub)
c     save crack files for post-process in matlab.

      use Global_Common
      use Global_Crack
      use Global_Model
      use Global_Filename
      
      implicit none
      
      integer isub
      integer i,j
      character(200) c_File_name_1,c_File_name_2
      character(5) temp  
 

      print *,'    Saving coordinates of cracks......'
      write(temp,'(I4)') isub
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'crax'//'_'//ADJUSTL(temp)  
      c_File_name_2   =  trim(Full_Pathname)//'.'
     &             //'cray'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')    
      open(102,file=c_File_name_2,status='unknown')  
      do i=1,num_Crack
          write(101, '(50E20.12)') (Edge_Disposed_Crack(i,
     &                               j,1),j=1,Each_Cr_Poi_Num(i))
          write(102, '(50E20.12)') (Edge_Disposed_Crack(i,
     &                               j,2),j=1,Each_Cr_Poi_Num(i))
      end do
      close(101)
      close(102)            
      !ennd file
      print *,'    Saving ennd file......'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'ennd'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')         
      do i=1,Num_Node
          write(101, '(20I10)') (Enriched_Node_Type(i,
     &                               j),j=1,num_Crack)
      end do
      close(101)           
      !elty file
      print *,'    Saving elty file......'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'elty'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')         
      do i=1,Num_Elem
          write(101, '(20I10)') (Elem_Type(i,j),j=1,num_Crack)
      end do
      close(101)    
      !celc file
      print *,'    Saving celc file......'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'celc'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')         
      do i=1,Num_Elem
          write(101, '(20E20.12)') (Coors_Element_Crack(i,j),j=1,4)
      end do
      close(101)    
      !celv file
      print *,'    Saving celv file......'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'celv'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')         
      do i=1,Num_Elem
          write(101, '(20E20.12)') (Coors_Vertex(i,j),j=1,2)
      end do
      close(101)   
      !celj file
      print *,'    Saving celj file......'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'celj'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')         
      do i=1,Num_Elem
          write(101, '(20E20.12)') (Coors_Junction(i,j),j=1,4)
      end do
      close(101)         
      !celt file
      print *,'    Saving celt file......'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'celt'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')         
      do i=1,Num_Elem
          write(101, '(20E20.12)') (Coors_Tip(i,j),j=1,2)
      end do
      close(101)   
      !ctty file
      print *,'    Saving ctty file......'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'ctty'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')         
      do i=1,num_Crack
          write(101, '(20I10)') (Crack_Tip_Type(i,j),j=1,2)
      end do
      close(101)         
      !posi file
      print *,'    Saving posi file......'
      c_File_name_1   =  trim(Full_Pathname)//'.'
     &             //'posi'//'_'//ADJUSTL(temp)      
      open(101,file=c_File_name_1,status='unknown')         
      do i=1,Num_Node
          write(101, '(200I10)') (c_POS(i,j),j=1,num_Crack)
      end do
      close(101) 

      RETURN
      END SUBROUTINE Save_Files_Crack
