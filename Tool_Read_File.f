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
 
      subroutine Tool_Read_File(Full_Name,Case_Type,m,n,Temp_DATA,
     &                          Flag_Blank)
C     read file
      implicit none

      character*200 Full_Name
      character*10  Case_Type
      integer i,j,m,n,istat
      logical Flag_Blank
      double precision Temp_DATA(m,n)
      
      select case(Case_Type(1:4))
      case('node')    
          print *, "    Reading nodal files...."
      case('elem')    
          print *, "    Reading element files...." 
      case('boux')    
          print *, "    Reading boux files...." 
      case('bouy')    
          print *, "    Reading bouy files...." 
      case('focx')    
          print *, "    Reading focx files...." 
      case('focy')    
          print *, "    Reading focy files...." 
      end select
      
      open(11,file=Full_Name,status='old')
      Read(11,*,IOSTAT=istat)
      close(11)
      if ( istat /= 0 ) then 
          Flag_Blank = .True.
      else                     
          Flag_Blank = .False.
          open(12,file=Full_Name,status='old')
          read(12,*)((Temp_DATA(i,j),j=1,n),i=1,m)
          close(12)
      end if
      
      RETURN
      END SUBROUTINE Tool_Read_File