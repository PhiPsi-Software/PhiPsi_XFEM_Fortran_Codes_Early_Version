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
 
      subroutine Cal_F(r,theta,omega,c_mat_type,F)
     
      use Global_Elem_Area
      implicit none
      
      double precision,intent(in)::r,theta,omega
      integer,intent(in)::c_mat_type
      double precision,intent(out)::F(4)
      
      double precision fac,st2,ct2,s3t2,c3t2,st,ct,
     &                 r2,c_theta
     
      c_theta = theta
      if (r.ne.0.0D0) then
          r2 = sqrt(r)
      else
          r2 = sqrt(Ave_Elem_Area)*0.1D-4
	      c_theta = 0.0d0
      end if
      
      ! ------------------------------
      !         ISO material.
      ! ------------------------------
      if (c_mat_type ==1 ) then
	      fac = 0.5D0/r2
	      st2 = sin(c_theta/2.0D0)
	      ct2 = cos(c_theta/2.0D0)
	      s3t2= sin(1.5D0*c_theta)
	      c3t2= cos(1.5D0*c_theta)
	      st  = sin(c_theta)
	      ct  = cos(c_theta)

          ! Tip enrichment functions F.
	
		  F(1) = r2*st2
		  F(2) = r2*ct2
		  F(3) = r2*st2*st
		  F(4) = r2*ct2*st
      ! ------------------------------
      !       Orthotropic material.
      ! ------------------------------
      elseif (c_mat_type ==2 .or. c_mat_type ==3) then
      
      end if
      
      return 
      end SUBROUTINE Cal_F                
