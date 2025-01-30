c--------------------------------------------------------------------c
c                                                                    c
c CAI STORAGE                                                        c
c v18                                                                c
c                                                                    c
c Steve Desch                                                        c
c August 28, 2017                                                    c
c                                                                    c
c Calculates radial transport of solids before and after formation   c
c  of Jupiter.                                                       c
c                                                                    c
c 1 = small, CI composition             a = 1 um -- processed        c
c 2 = small, CI composition             a = 1 um -- not processed    c
c 3 = small, CI-CAI composition         a = 1 um                     c 
c 4 = large, CI-CAI composition (AOA?)  a = 300 um                   c 
c 5 = small CAIs                        a = 400 um                   c 
c 6 = medium CAIs                       a = 800 um                   c 
c 7 = large CAIs                        a =1600 um                   c 
c                                                                    c
c--------------------------------------------------------------------c

       implicit none 
       integer ir,nr,nrmax,ir1,ir3,im2
       parameter(nrmax=500)
       double precision Rin,Rout,r(0:nrmax),rmid(nrmax),Area(nrmax)
       double precision Sigma(nrmax),rhog(nrmax),P(nrmax)  
       double precision eta(0:nrmax) 
       double precision Mass,Mass0
       double precision Temp0,q0,Tpassive,mbar,Temp(nrmax),C(nrmax)
       double precision Mstar,Mdisk,Omega(nrmax),H(nrmax)
       double precision alphalow,alphahi,alphagap,alphaMRI,salpha
       double precision alpha(nrmax),nu(nrmax) 
       double precision Mdot(0:nrmax),Vgr(0:nrmax),Q(0:nrmax) 
       double precision gamma,R1,tvisc,tau
       double precision MJup,RJ,dR,x,rCAI,T1,tJup,xp,tauJup
       double precision Sfloor,cfloor,ravg,Sigmaavg,Stavg,Diffavg
       double precision Stmin,Stmax,Sigmacavg,ccavg,Sc0,Sigma0
       integer is,ns,nsmax
       parameter(nsmax=10)
       double precision Massc(nsmax),Massc0(nsmax),Mouter,Minner
       double precision agr(nsmax),rhop(nsmax),cc0(nsmax)
       double precision St(nsmax,nrmax)
       double precision Sigmac(nsmax,nrmax),cc(nsmax,nrmax)
       double precision Diffc(nsmax,nrmax)
       double precision dumax,Mdotc(nsmax,0:nrmax),du(nsmax,0:nrmax)   
       double precision ur(0:nrmax),Vgr0,dlnTdlnr,dlnasdlnr 
       double precision rd,sigrd,xi,Mlost,Moutside,Mcell,Mdotpe,A0
       double precision rzin,rzout 
       integer it,jt,nt,np 
       double precision time,dtime,tprint,tsim,dtmax 
       double precision AU,Msol,MEarth,yr,Myr,um
       double precision ssb,Gconst,k,amu,pi,kappa 
       parameter(AU=1.4959787d+13,um=1.0d-04) 
       parameter(Msol=1.9891d+33,Mearth=5.98e+27) 
       parameter(yr=3.155d+07,Myr=3.155d+13)
       parameter(Gconst=6.67408d-08,k=1.38065d-16,amu=1.66054d-24)
       parameter(ssb=5.671d-05) 
       parameter(pi=3.14159265d0)
       character junk*25, filename*25
       integer yesno

c--------------------------------------------------------------------c
c Initialize                                                         c 
c                                                                    c 
c Code reads input parameters from file. Option at end of file is to c 
c  either initialize surface density from scratch ("2") or read in   c 
c  values from a file ("1"), name of input file specified.           c 
c                                                                    c 
c--------------------------------------------------------------------c

c READ IN INPUT PARAMETERS 
       open(21,file='code18input1.dat',form='formatted',status='old')
       read(21,503) junk,Mstar
        Mstar = Mstar * Msol 
       read(21,503) junk,Mdisk 
         Mdisk = Mdisk * Msol 
       read(21,501) junk,nr
       read(21,502) junk,Rin
        Rin  = Rin *AU
       read(21,502) junk,Rout
       Rout = Rout*AU
       read(21,501) junk,ns 
       do is=1,ns
        read(21,503) junk,agr(is) 
        agr(is) = agr(is)*um 
        read(21,503) junk,rhop(is) 
        read(21,503) junk,cc0(is) 
       end do 
       read(21,503) junk,time
        time = time*Myr 
       read(21,503) junk,tsim
        tsim = tsim*Myr 
       read(21,503) junk,dtime
        dtime = dtime*yr 
       read(21,503) junk,tprint
        tprint = tprint*yr 
       read(21,501) filename,yesno
501    format(A25,I3) 
502    format(A25,F6.3) 
503    format(A25,E9.3) 
511    format(E10.4,2X,I3,10(2X,E10.4)) 
       close(21) 

c CREATE RADIAL GRID 
c ir1 = index at which r closest to 1 AU
c ir3 = index at which r closest to 3 AU
       r(0) = Rin 
       do ir=1,nr
        r(ir) = Rin*dexp(+(dble(ir)/dble(nr))*dlog(Rout/Rin))
        rmid(ir) = dsqrt(r(ir-1)*r(ir)) 
        Area(ir) = pi*(r(ir)*r(ir)-r(ir-1)*r(ir-1)) 
       end do
       do ir=1,nr-1 
        if ( (rmid(ir  ).le.(1.*AU)) .and. 
     &       (rmid(ir+1).ge.(1.*AU)) ) then 
         ir1 = ir 
        end if 
        if ( (rmid(ir  ).le.(3.*AU)) .and. 
     &       (rmid(ir+1).ge.(3.*AU)) ) then 
         ir3 = ir 
        end if 
       end do 

c IF OPTION 1, READ IN Temp, Sigma, Sigmac
       if (yesno.eq.1) then
        open(23,file=filename,form='formatted',status='old')
         write(*,*) ' Opening file ',filename
         write(*,*) ' nr = ',nr 
         do ir=1,nr
          read(23,511) T1,ir1,R1,Sigma(ir),Temp(ir),
     & cc(1,ir),cc(2,ir),cc(3,ir),cc(4,ir),cc(5,ir),cc(6,ir),cc(7,ir)  
          do is=1,ns
c          cc(is,ir) = Sigmac(is,ir) / Sigma(ir) 
           Sigmac(is,ir) = cc(is,ir)*Sigma(ir) 
          end do 
c         write(*,*) ir, 
c    &                 T1,ir1,
c    &                 R1,
c    &                ' Sigma = ',Sigma(ir),
c    &                ' T = ',Temp(ir)
c    &                ' S1 = ',Sigmac(1,ir),
c    &                ' S2 = ',Sigmac(2,ir),
c    &                ' S3 = ',Sigmac(3,ir)  
c    &                ' S4 = ',Sigmac(4,ir)  
         end do 
        close(23) 
c       stop 
       end if

c IF OPTION 2, GENERATE Temp, Sigma, Sigmac
c Initialize Sigma with self-similar profile of Hartmann et al. (1998).
c Assume uniform distribution of CAIs at cc0, out to radius rCAI.
c Initialize Temp with a simple passive disk form like Chiang & Goldreich (1997).

       if (yesno.eq.2) then 
        Temp0 = 150.
        q0    = 0.428571d0
        gamma = 1.5d0 -q0
        do ir=1,nr
         Tpassive = Temp0*(rmid(ir)/AU)**(-q0) 
         Temp(ir) = Tpassive 
        end do 
c       R1    = r(ir1)
c DEBUG
        R1 = 1.*AU 
        tau = 1.0d0 
c would need to define nu(ir1) first! 
c       tvisc = R1*R1/nu(ir1) / ( 3.*(2.-gamma)*(2.-gamma) )
c       tau   = 1.0d0 + time/tvisc 
c       write(*,*) Mdisk,Mdisk*(2.-gamma)/(2.*pi*R1*R1) 
        do ir=1,nr 
         Sigma(ir) = Mdisk * (2.-gamma) / (2.*pi*R1*R1) 
     &                * (rmid(ir)/R1)**(-gamma)
     &                * tau**(-(2.5-gamma)/(2.0-gamma)) 
     &                * dexp( -(rmid(ir)/R1)**(2.-gamma) / tau ) 
         Sfloor = 1.0d-19
         if (Sigma(ir).lt.Sfloor) then
          Sigma(ir) = Sfloor
         end if 
         write(*,*) ir,rmid(ir)/AU,Sigma(ir) 
         do is=1,ns
          Sigmac(is,ir) = cc0(is)*Sigma(ir) 
          cc(is,ir) = cc0(is)
         end do 
        end do 
       end if 
c INITIALIZE OTHER VARIABLES 

       q0    = 0.428571d0
       MJup  = 30.*MEarth
       mbar  = 2.33*amu 
c      kappa = 10.  
c      kappa = 8.0
       kappa = 5.0 

c ALPHA 
c      alphahi =  3.0d-04
c      alphalow = 3.0d-05 
c      alphaMRI = 2.0d-03 
c      alphahi  = 1.25d-04 
       alphahi = 5.00d-04
c      alphahi = 1.00d-03
       alphalow = 1.00d-05 
       alphaMRI = 5.00d-04

       alphagap = 1.0d-02 
       tJup = 0.20*Myr 
       RJ       = (3.0d0)*AU
       rd = 55.*AU
c       tauJup = (6.0d+04)*yr
c       tauJup = (1.25d+05)*yr 
        tauJup = (1.00d+05)*yr 
        tauJup = (1.10d+05)*yr 

       do ir=1,nr
        C(ir)    = dsqrt(k*Temp(ir)/mbar) 
        Omega(ir) = dsqrt(Gconst*Mstar)*rmid(ir)**(-1.5d0) 
        H(ir)     = C(ir) / Omega(ir) 
        P(ir)    = Sigma(ir)*Omega(ir)*C(ir)/dsqrt(2.*pi) 
        rhog(ir) = Sigma(ir)*Omega(ir)/C(ir)/dsqrt(2.*pi)
       end do 
 
       Mass0 = 0.
       do is=1,ns
        Massc0(is) = 0. 
       end do 
       do ir=1,nr
        Mass0 = Mass0 + Sigma(ir)*Area(ir)
        do is=1,ns
         Massc0(is) = Massc0(is) + Sigmac(is,ir)*Area(ir) 
        end do
        write(*,514) rmid(ir)/AU,Sigma(ir),Temp(ir),
     & cc(1,ir),cc(2,ir),cc(3,ir),cc(4,ir),cc(5,ir),cc(6,ir),cc(7,ir)
514     format(E10.4,9(2X,E10.4)) 
       end do 

c--------------------------------------------------------------------c
c Output                                                             c 
c--------------------------------------------------------------------c

       open(22,file='code18output1.dat',form='formatted',status='old')
       open(24,file='code18jupiter1.dat',form='formatted',status='old')
       do ir=1,nr
        write(22,511) time/Myr,ir,rmid(ir)/AU,Sigma(ir),Temp(ir),
     & cc(1,ir),cc(2,ir),cc(3,ir),cc(4,ir),cc(5,ir),cc(6,ir),cc(7,ir) 
       end do 
       write(24,524) time/Myr,MJup/MEarth  
524    format(E12.6,2X,E12.6) 

c--------------------------------------------------------------------c
c Evolve over time                                                   c 
c--------------------------------------------------------------------c

       np     = int(tprint/dtime+0.001)
       nt     = int(tsim/tprint+0.001)

       Minner = 0.0d0
       Mouter = 0.0d0

       do it=1,nt
        do jt=1,np

        time = time + dtime 

c Calculate masses 
        Mass = 0.
        do ir=1,nr
         Mass = Mass + Sigma(ir)*Area(ir)
        end do 

        do is=1,ns
         Massc(is) = 0.0d0
         do ir=1,nr
          Massc(is) = Massc(is) +Sigmac(is,ir)*Area(ir)
         end do 
        end do 

c Calculate temperatures and viscosities 
c Passive disk temperature using Chiang & Goldreich (1997) and a 
c  luminosity that decreases in time, from Baraffe et al. (2002)
c 
c DEBUG - changed order of alphahi, alphalo, made limits 1.0 -> 2.0 AU
       do ir=1,nr
        rzin = (1.0)*AU
        rzout= (10.0)*AU
        Tpassive = (171.4)*(time/(1.*Myr))**(-0.142857d0)
     &                  *(rmid(ir)/AU)**(-q0) 
        if (rmid(ir).le.(rzin)) then
         alpha(ir) = alphahi
        end if 
        if (rmid(ir).gt.(rzout)) then 
         alpha(ir) = alphalow
        end if 
        if ((rmid(ir).gt.(rzin)).and.(rmid(ir).le.(rzout))) then 
         salpha = dlog(alphalow/alphahi) / dlog(rzout/rzin)
         alpha(ir) = alphahi*(rmid(ir)/(rzin))**salpha 
        end if 

        Temp(ir) = ( (27./128.)*Sigma(ir)*Sigma(ir)*kappa*alpha(ir)*  
     &               k*Omega(ir)/(mbar*ssb) )**(1./3.) 
        Temp(ir) = ( Temp(ir)**4. + Tpassive**4. )**(0.25d0) 

c-> CAI FACTORY
         if (Temp(ir).gt.1000.) then
          do is=1,10
           x = (Temp(ir)-1000.)/400.
           x = min(x,1.0d0) 
c DEBUG -> changed alphalow to alphahi
           alpha(ir) = alphahi*(alphaMRI/alphahi)**x
           Temp(ir) =((27./128.)*Sigma(ir)*Sigma(ir)*kappa*alpha(ir)*    
     &                 k*Omega(ir)/(mbar*ssb) )**(1./3.) 
           Temp(ir) = ( Temp(ir)**4. + Tpassive**4. )**(0.25d0) 
          end do 
         end if 

c         do is=1,10
cc         x = (Temp(ir)-1000.) / 400.
c          if (x.gt.1.0) then
c           x = 1.0
c          end if 
cc OK because in CAI Factory, alpha = alphalow
c          alpha(ir) = alphalow*(alphaMRI/alphalow)**x
c          Temp(ir) = ( (27./128.)*Sigma(ir)*Sigma(ir)*kappa*alpha(ir)*  
c     &                k*Omega(ir)/(mbar*ssb) )**(1./3.) 
c          Temp(ir) = ( Temp(ir)**4. + Tpassive**4. )**(0.25d0) 
c         end do   
c        end if 
c
cc       if (Temp(ir).lt.Tpassive) then
cc        Temp(ir) = Tpassive
c        end if 
c

        if (Temp(ir).ge.1400.) then
         Temp(ir) = 1400. 
        end if 

       end do 

c Gas pressure and density 
       do ir=1,nr
        C(ir)    = dsqrt(k*Temp(ir)/mbar) 
        H(ir)     = C(ir) / Omega(ir) 
        P(ir)    = Sigma(ir)*Omega(ir)*C(ir)/dsqrt(2.*pi) 
        rhog(ir) = Sigma(ir)*Omega(ir)/C(ir)/dsqrt(2.*pi)
       end do 
       do ir=1,nr-1
        eta(ir) = - ( P(ir+1)-P(ir) )/( rmid(ir+1)-rmid(ir) ) 
     &              / dsqrt(rhog(ir+1)*rhog(ir)) 
     &            *r(ir)*r(ir) / (Gconst*Mstar)
       end do 
       eta(0)  = eta(1)
       eta(nr) = eta(nr-1) 

c CALCULATE VISCOSITIES, INCLUDING GAP 
c If t > 1 Myr we assume alpha ~ 0.01 within a few Hill radii around
c  where Jupiter purportedly forms at 3 AU.  

        dR       = 1.* RJ*(MJup/3./Mstar)**(1./3.) 
        if (time.ge.(tJup)) then
         do ir=1,nr
          x = (rmid(ir)-RJ) / dR 
          alpha(ir) = alpha(ir) + (alphagap-alpha(ir))*dexp(-x*x) 
         end do 
        end if  

        do ir=1,nr
         nu(ir) = alpha(ir)*C(ir)*H(ir) 
        end do 

c UPDATE WITH GAS ADVECTION 
c Vgr calculated using formula with Q
c Donor cell approximation used 
c

       do ir=1,nr-1 
        Q(ir) = log((Sigma(ir+1)*nu(ir+1))/(Sigma(ir)*nu(ir))) 
     &         / log( rmid(ir+1)/rmid(ir) ) 
        Vgr(ir) = -1.5*((nu(ir+1)+nu(ir))/2.)/r(ir)*(1.+2.*Q(ir))  
        dtmax = (r(ir)-r(ir-1))/dabs(Vgr(ir)) 
        ravg  = r(ir) -Vgr(ir)*dtime
        x = (ravg  - rmid(ir)) / (rmid(ir+1)-rmid(ir))
         if ((x.lt.0.).or.(x.gt.1.)) then
          write(*,*) ' Courant condition violated w/ Vgr'
          stop
         end if 
c       Sigmaavg  = Sigma(ir)  + (Sigma(ir+1) -Sigma(ir) )*x 
        if (Vgr(ir).ge.0.) then 
         Sigmaavg = Sigma(ir)
         Mdot(ir)  = -2.*pi*r(ir)*Vgr(ir)*Sigmaavg 
         do is=1,ns
          ccavg = cc(is,ir)
          Mdotc(is,ir) = Mdot(ir)*ccavg 
         end do 
        else
         Sigmaavg = Sigma(ir+1)
         do is=1,ns
          ccavg = cc(is,ir+1)
          Mdotc(is,ir) = Mdot(ir)*ccavg 
         end do 
        end if 
        Mdot(ir)  = -2.*pi*r(ir)*Vgr(ir)*Sigmaavg 
c       do is=1,ns
c Should I be doing this? 
c        ccavg = cc(is,ir) +(cc(is,ir+1)-cc(is,ir))*x    
c        Mdotc(is,ir) = Mdot(ir)*ccavg 
c       end do 
       end do 

c At inner boundary, particles advected inward with gas 
       Vgr(0) = -1.5*(nu(1))/r(0)*(1.+2.*Q(1))
       if (Vgr(0).lt.0.) then
        Mdot(0)  = -2.*pi*r(0)*Sigma(1)*Vgr(0) 
        do is=1,ns
         Mdotc(is,0) = Mdot(0)*cc(is,1) 
        end do 
       else
        Mdot(0) = 0.0d0 
        do is=1,ns
         Mdotc(is,0) = 0.0d0 
        end do 
       end if 

c At outer boundary, particles advected outward with gas
c      Vgr(nr) = -1.5*(2.*nu(nr)-nu(nr-1))/r(nr)*(1.+2.*Q(nr-1))  
c      Vgr(nr) = -1.5*(nu(nr))/r(nr)*(1.+2.*Q(nr-1))
c DEBUG 
c Because of photoevaporation, we don't need to let them cross the outer boundary.
       Vgr(nr) = 0. 
       Mdot(nr) = 0.
       do is=1,ns
        Mdotc(is,ir) = 0. 
       end do 
c
c      if (Vgr(nr).lt.0.) then
c       Mdot(nr) = 0.
c       do is=1,ns
c        Mdotc(is,nr) = 0.
c       end do 
c      else 
c       Mdot(nr) = -2.*pi*r(nr)*Sigma(nr)*Vgr(nr)
c       do is=1,ns
c        Mdotc(is,nr) = Mdot(nr)*cc(is,nr)
c       end do 
c      end if 

c Update surface densities accordingly 
       do ir=1,nr
        Sigma(ir)   = Sigma(ir) +dtime*(Mdot(ir)-Mdot(ir-1))/Area(ir)  
        Sfloor = 1.0d-19
        if (Sigma(ir).lt.Sfloor) then
         Sigma(ir) = Sfloor 
        end if 
        do is=1,ns
         Sigmac(is,ir) = Sigmac(is,ir) 
     &                  +dtime*(Mdotc(is,ir)-Mdotc(is,ir-1))/Area(ir)
         if (Sigmac(is,ir).lt.0.) then
          Sigmac(is,ir) = 0.
         end if
        end do 
       end do 

       Minner = Minner + Mdot(0) *dtime 

c UPDATE FOR PARTICLE DRIFT 

       Sc0  = 0.7d0
       dumax = 1.0d+04 

       do is=1,ns
        do ir=1,nr
         St(is,ir) = Omega(ir)*rhop(is)*agr(is)/rhog(ir)/C(ir)  
         Stmin = 1.0d-05 
         Stmax = 1.0d+03 
         if (St(is,ir).ge.Stmax) then
          St(is,ir) = Stmax
         end if
         Diffc(is,ir) = nu(ir)/Sc0 
     &             / ( 1.0d0 + St(is,ir)*St(is,ir) ) 
        end do 
       end do 

       do ir=1,nr-1
        dlnasdlnr = log( (alpha(ir+1)*Sigma(ir+1)) / 
     &                   (alpha(ir  )*Sigma(ir  )) ) / 
     &              log( rmid(ir+1) / rmid(ir) ) 
        dlnTdlnr = log( Temp(ir+1) / Temp(ir) ) / 
     &              log( rmid(ir+1) / rmid(ir) ) 
        ur(ir) = 0.5*(nu(ir)+nu(ir+1))/r(ir)*
     &           (-1.5d0 -3.*dlnasdlnr -0.5*dlnTdlnr)
       end do 
        ur(1)  = 0.0
        ur(nr) = 0.0

c Vgr0 accounts for meridional flow 
       do is=1,ns
        do ir=1,nr-1
         Stavg = (St(is,ir)+St(is,ir+1))/2. 
         Vgr0 = 0.0d0
         if (Stavg.gt.0.18*alpha(ir)) then 
          Vgr0 = (ur(ir) - Vgr(ir)) / (1.+Stavg*Stavg) 
         end if 
         du(is,ir) = -Stavg/(1.+Stavg*Stavg) *
     &        ( Vgr(ir)*Stavg +eta(ir)*dsqrt(Gconst*Mstar/r(ir)) )
     &      + Vgr0 
         if (du(is,ir).gt.0.) then
          Sigmacavg = Sigmac(is,ir)
         else 
          Sigmacavg = Sigmac(is,ir+1)
         end if 
         Mdotc(is,ir) = -2.*pi*r(ir)*du(is,ir)*Sigmacavg 
        end do

c Inner boundary 
        du(is,0) = du(is,1) 
        if (du(is,0).lt.0.) then
         Mdotc(is,0) = -2.*pi*r(ir)*du(is,0)*Sigmac(is,1)
        else
         Mdotc(is,0) = 0. 
        end if 
c Outer boundary 
c       du(is,nr) = du(is,nr-1)
        du(is,nr) = 0. 
        Mdotc(is,nr) = 0. 
c       if (du(is,0).gt.0.) then
c        Mdotc(is,nr) = -2.*pi*r(nr)*du(is,nr)*Sigmac(is,nr)
c       else
c        Mdotc(is,nr) = 0.
c       end if 

c Update surface densities 
         do ir=1,nr
          Sigmac(is,ir) = Sigmac(is,ir) 
     &                   +dtime*(Mdotc(is,ir)-Mdotc(is,ir-1))/Area(ir)
          if (Sigmac(is,ir).lt.0.) then
c          write(*,*) ' negative density because of du at ',is,ir 
           Sigmac(is,ir) = 0.
          end if 
         end do 

       end do

c UPDATE FOR PARTICLE DIFFUSION 
c For simplicity we don't allow particles to flow across the inner or outer boundaries 
c  due to diffusion (just advection and drift). 

       do is=1,ns

        do ir=1,nr-1
         Diffavg  = (Diffc (is,ir)+Diffc (is,ir+1))/2. 
         Sigmaavg = (Sigma(ir)+Sigma(ir+1))/2. 
         Mdotc(is,ir) = 2.*pi*r(ir)*Sigmaavg*Diffavg*  
     &           (cc(is,ir+1)-cc(is,ir)) / (rmid(ir+1)-rmid(ir))   
        end do 
        Mdotc(is,0) = 0. 
        Mdotc(is,nr) = 0. 
        do ir=1,nr
         Sigmac(is,ir) = Sigmac(is,ir) 
     &                  +dtime*(Mdotc(is,ir)-Mdotc(is,ir-1))/Area(ir)
         if (Sigmac(is,ir).lt.0.) then
c         write(*,*) ' negative density because of diffusion at ',is,ir 
          Sigmac(is,ir) = 0.
         end if 
        end do 

       end do 

c ACCOUNT FOR JUPITER FORMING 
        if (time.ge.(tJup)) then
c DEBUG 
        dR       = 1.* RJ*(MJup/3./Mstar)**(1./3.) 
c
         do ir=1,nr
          x = (rmid(ir)-RJ) / dR 
          MJup = MJup +Sigma(ir)*dtime/tauJup*Area(ir)*dexp(-x*x) 
          Sigma(ir) = Sigma(ir) -Sigma(ir)*dtime/tauJup*dexp(-x*x) 
          do is=1,ns
           Sigmac(is,ir) = Sigmac(is,ir) 
     &                    -Sigmac(is,ir)*dtime/tauJup*dexp(-x*x) 
           MJup =MJup +Sigmac(is,ir)*dtime/tauJup*Area(ir)*dexp(-x*x)    
          end do 
         end do 
        end if  
c      do ir=1,nr
c      end do 

c CONVERT CI MATERIAL INTO AOAs + CAIs + depleted material 
c DEBUG 
         do ir=1,nr
          if (Temp(ir).ge.1399.) then
c          if (time.le.(0.1*Myr)) then 
           Sigma0 = Sigmac(1,ir) 
           Sigmac(1,ir) = Sigmac(1,ir) - 0.1*Sigmac(1,ir)
           Sigmac(2,ir) = Sigmac(2,ir) 
           Sigmac(3,ir) = Sigmac(3,ir) + 0.1*Sigma0*0.89     
           Sigmac(4,ir) = Sigmac(4,ir) + 0.1*Sigma0*0.03 
           Sigmac(5,ir) = Sigmac(5,ir) + 0.1*Sigma0*0.08/3. 
           Sigmac(6,ir) = Sigmac(6,ir) + 0.1*Sigma0*0.08/3. 
           Sigmac(7,ir) = Sigmac(7,ir) + 0.1*Sigma0*0.08/3.
c         end if 
          end if 
         end do 

c PHOTOEVAPORATION?

c        x = (333.*AU)/rd 
c        Mdotpe = ((1.5d-07)*Msol/yr)*dexp(-x)*x 
c        x = (375.*AU)/rd
c        Mdotpe = ((1.57d-07)*Msol/yr)*x*dexp(-x) 
cc       x = (450.*AU)/rd
cc       Mdotpe = ((1.72d-07)*Msol/yr)*x*dexp(-x)
         x = (750.*AU)/rd
         Mdotpe = ((2.23d-07)*Msol/yr)*x*dexp(-x) 
c        x = (1500.*AU)/rd
c        Mdotpe = ((3.16d-07)*Msol/yr)*x*dexp(-x) 
c        x = (2500.*AU)/rd
c        Mdotpe = ((4.08d-07)*Msol/yr)*x*dexp(-x) 
 
c
         Mlost = Mdotpe*dtime
         Moutside = 0.0d0
         ir = nr+1
620      continue
         ir = ir-1
         Mcell = Sigma(ir)*Area(ir) 
         Moutside = Moutside + Mcell
         if (Moutside.gt.Mlost) then
          sigrd = (Moutside-Mlost)/Area(ir)
          xi = sigrd / Sigma(ir)
          rd = dsqrt( xi*r(ir)*r(ir) + (1.0-xi)*r(ir-1)*r(ir-1) ) 
          im2 = ir
         else 
          if (ir.gt.1) then
           goto 620
          end if
         end if

         do is=1,ns
c         Sigmac(is,im2) = (Sigmac(is,im2)/Sigma(im2))*sigrd
c         Sigmac(is,im2) = cc(is,im2)*sigrd 
c DEBUG
          Sigmac(is,im2) = 0. 
         end do
          Sigma(im2) = sigrd
         do ir=im2+1,nr
          Sigma(ir) = 0.
          do is=1,ns
           Sigmac(is,ir) = 0.
          end do 
         end do 

         Mouter = Mouter + Mlost 
c        write(*,*) time/yr,rd/AU,Mlost/Msol,Mouter/Msol 
c        write(*,*) time/yr,rd/AU,Mdotpe/Msol*yr 
      

c UPDATE CONCENTRATIONS
        Sfloor = 1.0d-19
        cfloor = 1.0d-24 
        do ir=1,nr
         if (Sigma(ir).lt.Sfloor) then
          Sigma(ir) = Sfloor 
         end if 
         do is=1,ns
          cc(is,ir) = Sigmac(is,ir) / Sigma(ir)
          if (cc(is,ir).lt.cfloor) then
           cc(is,ir) = cfloor 
          end if 
         end do 
        end do 
c DEBUG 
c        write(*,*) time/yr, Minner/Msol, Mouter/Msol,
c    &      Vgr(nr-2),Vgr(nr-1),Vgr(nr)
c    &      Mdot(nr)/Msol*yr 
c    &      Mdot(232)/Msol*yr 
c        write(*,*) time/yr, rd/AU 

c--------------------------------------------------------------------c
c End inner time loop                                                c 
c--------------------------------------------------------------------c

c       do ir=1,nr
c        write(*,*) time/yr,ir,
c    &       cc(1,ir),cc(5,ir)
c    &       cc(1,ir),cc(2,ir),cc(3,ir),cc(4,ir),
c    &                 cc(5,ir),cc(6,ir),cc(7,ir)
c       end do 
c        read(*,*) is 

        end do

c--------------------------------------------------------------------c
c Output                                                             c 
c--------------------------------------------------------------------c

        write(*,*) ' ' 
        write(*,*) ' time = ',time/Myr,' Myr'        
        write(*,*) ' Mass = ',Mass/(Msol)
        write(*,*) ' Min  = ',Minner/(Msol)
        write(*,*) ' Mout = ',Mouter/(Msol)
        write(*,*) ' MJup / ME = ',MJup/MEarth 
        write(*,*) ' Mtot = ',(Minner+Mouter+Mass+MJup)/(Msol)
        write(*,*) ' rd = ',rd/AU,' AU'
c       do is=1,ns
c        write(*,*) ' Mcai / Mcai0 = ',Massc(is)/Mass0 
c       end do 

       do ir=1,nr
        write(22,511) time/Myr,ir,rmid(ir)/AU,Sigma(ir),Temp(ir),
     & cc(1,ir),cc(2,ir),cc(3,ir),cc(4,ir),cc(5,ir),cc(6,ir),cc(7,ir)
       end do 
        write(24,524) time/Myr,MJup/MEarth 

c--------------------------------------------------------------------c
c End time loop                                                      c 
c--------------------------------------------------------------------c
       end do 
        
      write(*,*) ' time = ',time/Myr,' Myr' 
      do ir=1,nr
        write(*,514) rmid(ir)/AU,Sigma(ir),Temp(ir),
     &  cc(1,ir),cc(2,ir),cc(3,ir),cc(4,ir),cc(5,ir),cc(6,ir),cc(7,ir)
      end do 

       close(22) 
       close(24) 




      end 

