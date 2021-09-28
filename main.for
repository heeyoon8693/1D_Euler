! =====================================================================
! 					1D Euler Sod shock tube problem
! 				Compile by silverfrost FTN95 personal edition 
! 				computer spec : windows7, 32bit 
!   				 TVD-MUSCL(3rd) with AUSM	
!					made by Hee Yoon(Republic of Korea) 
! 					recent update : 2021.09.13
! =====================================================================
! 						declare program 
!----------------------------------------------------------------------
      program Main
      use parameter_mod 
      call initial 
      iter = 0
10	  continue

	  iter = iter + 1

      call AUSM 
      call update

      if(iter.ne.iend)  goto 10
         
      call output
      stop
      end

! =====================================================================
! 							output data to Excel
! =====================================================================
	  subroutine output 
      use parameter_mod
      character*1 z
      z = char(9)

      open(2,file='TEST_AUSM_20210913.xls')

      write(2,*)'x',z,'r',z, 'u',z, 'p',z, 'e',z, 'ru'
      do 100 i=1,im
100   write(2,200)x(i),z,rho(i),z,vel(i),z,pre(i),z,e(i),z,q(i,2)
200   format(5(1pe11.4,a1),1pe11.4)
      close(2)
      return
      end


! =====================================================================
! 							declare initial condition
! =====================================================================
      subroutine initial 
      use parameter_mod

      dx = leng_x/real(im1)
      do i =1, im
        x(i) = dx*dble(i-1)
      continue
      if(x(i).le.x_diap) then
        rho(i) = rhol_ini
        vel(i) = vell_ini
        pre(i) = prel_ini
        e(i)   = prel_ini/gm1 +.5*rhol_ini*(vell_ini**2)
      else
       	rho(i) = rhor_ini
        vel(i) = velr_ini
        pre(i) = prer_ini
        e(i)   = prer_ini/gm1 +.5*rhol_ini*(velr_ini**2)
      end if
      continue
	    q(i,1) = rho(i)
        q(i,2) = rho(i)*vel(i)
        q(i,3) = e(i)
      end do 
      return
      end subroutine initial 

! =====================================================================
! 								update parameters
! =====================================================================
      subroutine update
      use parameter_mod

      do i = 2, im1
        rho(i) = q(i,1)
        vel(i) = q(i,2)/q(i,1)
        e(i)   = q(i,3)
        pre(i) = gm1*(e(i)-.5*rho(i)*vel(i)*vel(i))
      end do 
      return
      end subroutine update

! =====================================================================
! 									main solver
! =====================================================================
      subroutine AUSM
      use parameter_mod
      dimension fp(3), fm(3)

      call flux_AUSM(1,fm)

      do 100 i = 2, im1
        call flux_AUSM(i,fp)
        do 10 k = 1, 3
       	 q(i,k) = q(i,k) - dt/dx*(fp(k)-fm(k))
         fm(k)  = fp(k)
10		continue
100	  continue
      
      return
      end subroutine AUSM


! =====================================================================
! 				    	Flux calculation with MUSCL
! =====================================================================
      subroutine flux_AUSM(i,f)
      use parameter_mod
      real :: rl, rr, ul, ul, pl, pr
      dimension f(3)

      imm = i-1
      ipp = i+1
      ip2 = i+2

      if(i.eq.1 .or. i.eq.im1) then
       rl = rho(i)
       ul = vel(i)
       pl = pre(i)
       rr = rho(ipp)
       ur = vel(ipp)
       pr = pre(ipp)
      else
       call muscl(rho(imm),rho(i),rho(ipp),rho(ip2),  rl,  rr)
       call muscl(vel(imm),vel(i),vel(ipp),vel(ip2),  ul,  ur)
       call muscl(pre(imm),pre(i),pre(ipp),pre(ip2),  pl,  pr)
      endif


      el  = pl/gm1+.5*rl*(ul*ul)
      er  = pr/gm1+.5*rr*(ur*ur)
      hl  = el+pl
      hr  = er+pr
      al  = sqrt(gm*pl/rl)
      ar  = sqrt(gm*pr/rr)
      aml = ul/al
      amr = ur/ar

      if(abs(aml).le.1.) then
       amlp = .25*(aml+1.)**2
       plp  = .25*pl*(1.+aml)**2*(2-aml)
      else
       amlp = .5*(aml+abs(aml))
       plp  = .5*pl*(aml+abs(aml))/aml
      endif

      if(abs(amr).le.1.) then
       amrm = -.25*(amr-1.)**2
       prm  = .25*pr*(1.-amr)**2*(2+amr)
      else
       amrm = .5*(amr-abs(amr))
       prm  = .5*pr*(amr-abs(amr))/amr
      endif

      amc = amlp+amrm
      rul = rl*ul
      rur = rr*ur
      f(1)=.5*amc*( rr*ar+ rl*al)-.5*abs(amc)*( rr*ar- rl*al)
      f(2)=.5*amc*(rur*ar+rul*al)-.5*abs(amc)*(rur*ar-rul*al)+(plp+prm)
      f(3)=.5*amc*( hr*ar+ hl*al)-.5*abs(amc)*( hr*ar- hl*al)

      return
      end

! =======================================================================
!	3rd order TVD-MUSCL scheme with minmod limiter 
!  Reference "Yamamoto, S. and Daiguji, H. (1993). Higher-order-accurate 
!             upwind schemes for solving the compressible Euler and 
!             Navier-Stokes equations. Computers & Fluids, 22(2/3), 259-270"
! =======================================================================
 	  subroutine muscl(qm,qi,qp,qp2,ql,qr)
      use parameter_mod
        
      if(iorder.eq.1) then
       phi = 0.
       kpa = 0.
      endif
      if(iorder.eq.2) then
       phi = 1.
       kpa =-1.
      endif
      if(iorder.eq.3) then
       phi = 1.
       kpa = 1./3.
      endif

      beta = (3.d0-kpa)/(1.d0-kpa)
      phi4 = phi*.25
      kp   = 1.+kpa
      km   = 1.-kpa

      dq1  = qi -qm
      dq2  = qp -qi
      dq3  = qp2-qp

      gbq  = aminmod(dq1,beta*dq2)
      dbq  = aminmod(dq2,beta*dq1)
      ql   = qi+phi4*(km*gbq+kp*dbq)

      gbq  = aminmod(dq2,beta*dq3)
      dbq  = aminmod(dq3,beta*dq2)
      qr   = qp-phi4*(kp*gbq+km*dbq)

      return
      end

!=======================================================================
!	Calculate minmod    
!=======================================================================
 	  function   aminmod(xx,yy)
      real xx,yy,aminmod,aa,bb,cc,dd

      aa = sign(1.,xx)
      bb = abs(xx)
      cc = aa*yy
      dd = min(bb,cc)
      aminmod = aa*max(0.d0,dd)

      return
      end

