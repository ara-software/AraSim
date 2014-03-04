c Given Q(z) longitudinal development of excess charge
c computes vector potential at observer position (Fresnel) 
c in time domain. 
c 
c Calculates x,y,z components of vector potential
c Accounts for polarization dependence on shower depth
c
c 2D-model i.e. the lateral distribution is accounted for
c by a parameterisation of the vector potential. 
c 
c S.I. units
c Based on "Practical and accurate calculations of Askaryan radiation" PRD in press (2012).
c J. Alvarez-Muniz, A. Romero-Wolf, E. Zas
c
c Uses ZHAireS input
      program model
      implicit double precision (a-h,o-z)
      parameter(pi=3.1415926)
      parameter(c=2.99792458e8)    ! speed of light m^s-1
      parameter(xmu=12.566370e-7)  ! magnetic permeability N A^-2 
      parameter(e=1.602177e-19)    ! electron charge C 

c Separated files for electron and positron longitudinal profiles
      character*80 fileine,fileinp

      common/filename/fileine,fileinp

      common/antennae/xa(100),ya(100),za(100),rhoa(100),Ra(100)
      common/antenna_index/iant

      common/depthmax/zend
      common/ice/xn,x0,rho
      common/cherenkov/cher
      common/const/factor
      common/time/tobs_sec
      common/charge/xntot
      common/shower/E_TeV
      common/units/s_to_ns,xns_to_s

c Vegas routine commons
      CHARACTER NAMEOU*12
      COMMON/BVEG1/NCALL,ITMX,NPRN,NDEV,XL(18),XU(18),ACC
      COMMON/NAME/NAMEOU
      external vegas,xintegrand_x,xintegrand_y,xintegrand_z

      ncall=10000
      itmx=30

c Vegasd
      nameou='vegas.out'

c Output file containing waveforms of vector potential vs observer's time
      open (unit=22,status='unknown',file='output.dat')    

      write(*,*) 'File containing long shower profile (electrons)'
      read(*,*) fileine
      write(*,*) fileine
      write(*,*) 'File containing long shower profile (electrons)'
      read(*,*) fileinp
      write(*,*) fileinp

      write(*,*) 'Shower energy [TeV]'
      read(*,*) E_TeV
      write(*,*) E_TeV

c ICE      
      xn=1.78   ! refractive index 
      x0=36.08  ! radiation length g cm^-2
      rho=0.924 ! density g cm^-3

      write(*,*) 'Number of antennas'
      read(*,*) nant
      write(*,*) nant


c Position of shower maximum in g cm^-2 if antenna
c needs to be referred to zmax
      zmax=0.
      zmax=zmax/rho/100.               ! m

      write(*,*) 'Position of antenna i'
      write(*,*) '(x,y,z) in m'
      write(*,*) 'Primary injected at (0,0,0)'
      do ia=1,nant
        read(*,*) xa(ia),ya(ia),za(ia)
        write(*,*) xa(ia),ya(ia),za(ia)
c Convert antenna positions from g cm^-2 to meters      
c        xa(ia)=xa(ia)/rho/100.               ! m
c        ya(ia)=ya(ia)/rho/100.               ! m
c        za(ia)=za(ia)/rho/100.               ! m
        rhoa(ia)=sqrt(xa(ia)*xa(ia)+ya(ia)*ya(ia))  ! m
        Ra(ia)=sqrt(rhoa(ia)*rhoa(ia)+
     #              (za(ia)-zmax)*(za(ia)-zmax))  ! m
c        write(55,*) ia-1,Ra(ia)
      end do

c First call to interpolation to initialize zend and xntot 
      dum=xnep(10.)
      write(*,*) 'Integral[N(z)dz] ',xntot

c Limits in integration
c Vegasd
c Integration limits in shower depth z' 
      xl(1)=0.
      xu(1)=zend

c Assumes shower travels at speed of light
      beta=1.
      cher=acos(1./xn)

c Factor in Eq.(22) PRD paper
      factor=-xmu/(4.*pi)  

      s_to_ns=1.e9
      xns_to_s=1.e-9

      write(*,*) 'Time resolution [ns]'
      read(*,*) dtbin
      write(*,*) dtbin
      write(*,*) 'Time min. and max. [ns]'
      read(*,*) tmin,tmax
      itmax=int((tmax-tmin)/dtbin)

c Loop in observer's time
c Compute vector potential at each antenna position
c Electric field can be computed by a simple derivative w.r.t. time.
      do iant=1,nant

        write(*,*) 
     #'# Time [ns] Vector potential [Vs] antenna position (x,y,z) m',
     #xa(iant),ya(iant),za(iant)
        write(22,*) 
     #'# Time [ns] Vector potential [Vs] antenna position (x,y,z) m',
     #xa(iant),ya(iant),za(iant)

        vp_x=0.
        vp_y=0.
        vp_z=0.

c Calculate x,y,z components of vector potential
        do it=1,itmax   
          tobs_plot=tmin+dtbin*(it-1)+dtbin/2.
          tobs=tobs_plot+(xn*Ra(iant)/c)*s_to_ns
          tobs_sec=tobs*xns_to_s

c x-component of Eq.(22) PRD paper
          if(xa(iant).eq.0.) then 
            vp_x=0.
          else
            call vegas(1,xintegrand_x,vp_x,sal_x,schi_x)
          end if

c y-component of Eq.(22) PRD paper
          if(ya(iant).eq.0.) then 
            vp_y=0.
          else
            call vegas(1,xintegrand_y,vp_y,sal_y,schi_y)
          end if

c z-component of Eq.(22) PRD paper
          call vegas(1,xintegrand_z,vp_z,sal_z,schi_z)

          vp_x=factor*vp_x
          vp_y=factor*vp_y
          vp_z=factor*vp_z

          write(*,*) tobs,vp_x,vp_y,vp_z
c          write(*,*) tobs_plot,vp_x,vp_y,vp_z
          write(22,*) tobs,vp_x,vp_y,vp_z
c          write(22,*) tobs_plot,vp_x,vp_y,vp_z
        end do
        write(*,*) 
        write(22,*) 

      end do

      end

c ----------------------------------------------------------------------------
c x-component of integral in Eq.(22) PRD paper
      function xintegrand_x(qww)
      implicit double precision (a-h,o-z)
      save
      common/x/x(1),xd(17)
      common/antennae/xa(100),ya(100),za(100),rhoa(100),Ra(100)
      common/antenna_index/iant
      common/cherenkov/cher
      common/time/tobs_sec
      common/dist/R

c Depth in m 
      z=x(1)  ! m

c Argument of F_p function. Eq.(22) PRD.
      arg=argument(z,tobs_sec)

c Distance from position in shower depth z' to each antenna. 
c Denominator in Eq. (22) PRD paper 
      R = sqrt(rhoa(iant)*rhoa(iant)+
     #        (za(iant)-z)*(za(iant)-z))  ! m

c x-component of v_perp
c Account for polarization dependence on shower position
c Eq. (22) PRD paper
      u_x=xa(iant)/R
      u_z=(za(iant)-z)/R
      beta_z=1.
      vperp_x=u_x*u_z*beta_z

      xintegrand_x=-vperp_x*xnep(z)*F_p(arg)/R

      return
      end

c ----------------------------------------------------------------------------
c y-component of integral in Eq.(22) PRD paper
      function xintegrand_y(qww)
      implicit double precision (a-h,o-z)
      save
      common/x/x(1),xd(17)
      common/antennae/xa(100),ya(100),za(100),rhoa(100),Ra(100)
      common/antenna_index/iant
      common/cherenkov/cher
      common/time/tobs_sec
      common/dist/R

c Depth in m 
      z=x(1)  ! m

c Argument of F_p function. Eq.(22) PRD.
      arg=argument(z,tobs_sec)

c Distance from position in shower depth z' to each antenna. 
c Denominator in Eq. (22) PRD paper
      R = sqrt(rhoa(iant)*rhoa(iant)+
     #        (za(iant)-z)*(za(iant)-z))  ! m

c y-component of v_perp
c Account for polarization dependence on shower position
c Eq.(22) PRD 
      u_y=ya(iant)/R
      u_z=(za(iant)-z)/R
      beta_z=1.
      vperp_y=u_y*u_z*beta_z

      xintegrand_y=-vperp_y*xnep(z)*F_p(arg)/R

      return
      end

c ----------------------------------------------------------------------------
c z-component of integral in Eq.(22) PRD paper
      function xintegrand_z(qww)
      implicit double precision (a-h,o-z)
      save
c Common for vegas routine
      common/x/x(1),xd(17)

      common/antennae/xa(100),ya(100),za(100),rhoa(100),Ra(100)
      common/antenna_index/iant
      common/cherenkov/cher
      common/time/tobs_sec
      common/dist/R

c Depth in m 
c Specify first and only variable to be integrated by vegas
      z=x(1)  ! m

c Argument of F_p function. Eq.(22) PRD.
      arg=argument(z,tobs_sec)

c Distance from position in shower depth z' to each antenna 
c Denominator in Eq. (22) PRD paper
      R = sqrt(rhoa(iant)*rhoa(iant)+
     #        (za(iant)-z)*(za(iant)-z))  ! m

c z-component of v_perp
c Eq.(22) PRD
      u_x=xa(iant)/R
      u_y=ya(iant)/R
      beta_z=1.
      vperp_z=-(u_x*u_x+u_y*u_y)*beta_z

      xintegrand_z=-vperp_z*xnep(z)*F_p(arg)/R

      return
      end


c ----------------------------------------------------------------------------- 
c Argument of F_p Eq.(22) PRD paper
c ----------------------------------------------------------------------------- 
      function argument(z,t)
      implicit double precision (a-h,o-z)
      save
      parameter(c=2.99792458e8)    ! speed of light m^s-1
      common/antennae/xa(100),ya(100),za(100),rhoa(100),Ra(100)
      common/antenna_index/iant
      common/cherenkov/cher
      common/ice/xn,x0,rho
      common/dist/R

c Assumes shower travels at speed of light
      beta=1.
      argument=z-(beta*c*t-xn*R)  ! m 

      return
      end

c ----------------------------------------------------------------------------- 
c Function F_p Eq.(15) PRD paper.
c ----------------------------------------------------------------------------- 
      function F_p(arg)
      implicit double precision (a-h,o-z)
      save
      parameter(pi=3.1415926)
      parameter(c=2.99792458e8)    ! speed of light m^s-1
      parameter(xmu=12.566370e-7)  ! magnetic permeability N A^-2 
      parameter(e=1.602177e-19)    ! electron charge C 

      common/cherenkov/cher
      common/charge/xntot
      common/shower/E_TeV
      common/units/s_to_ns,xns_to_s
      common/ice/xn,x0,rho

c Factor accompanying the F_p in Eq.(15) in PRD paper
      fc=4.*pi/(xmu*sin(cher))  

c Converts arg into time
      beta=1.
c Note that Acher peaks at tt=0 which corresponds to the observer time.
c The shift from tobs to tt=0 is done when defining argument
      tt=(-arg/(c*beta))*s_to_ns ! Parameterisation of A_Cherenkov with t in ns


c Parameterisation of vector potential at Cherenkov angle
c Eq.(16) PRD paper.
      Af=-4.5e-14  ! V s
      if (tt.ge.0.) then 
       Acher=Af*E_TeV*(exp(-abs(tt)/0.057)+(1.+2.87*abs(tt))**(-3))
      else
       Acher=Af*E_TeV*(exp(-abs(tt)/0.03)+(1.+3.05*abs(tt))**(-3.5))
      end if

c Cut fit above +/-5 ns
c      if (abs(tt).gt.5.) Acher=1.e-30

c Obtain "shape" of Lambda-function from vp at Cherenkov angle 
c xntot = LQ_tot in PRD paper
      F_p=Acher*fc/xntot

      return
      end

C-----------------------------------------------------------
c Simple interpolation in depth of longitudinal profile of excess charge.
C-----------------------------------------------------------
      function xnep(z)
      implicit double precision (a-h,o-z)
      save
      common/depthmax/zend
      dimension xx(1000),yexcess(1000),ylog(1000)
      dimension xxp(1000),ye(1000),yp(1000)
      character*80 fileine,fileinp
      logical firstgo
      data firstgo / .true. /
      common/filename/fileine,fileinp
      common/ice/xn,x0,rho
      common/charge/xntot
      nnint=1

      
      if (firstgo .eqv. .true.) then  
      xntot=0.
c Read longitudinal profile
      npoints=0
      open (unit=11,status='old',file=fileine)      
      open (unit=12,status='old',file=fileinp)      

      open (unit=14,status='unknown',file='long.dat')      

      do idum=1,34
        read(11,*) 
        read(12,*) 
      end do

      do j=1,1000
        read(11,*,err=11,end=11) dum1,xx(j),ye(j)
        npoints=npoints+1
      end do
11    close(unit=11)

      do j=1,1000
        read(12,*,err=12,end=12) dum1,xxp(j),yp(j)
      end do
12    close (unit=12)


      do j=1,npoints
        yexcess(j)=ye(j)-yp(j)
        write(14,*) xx(j),ye(j),yp(j),yexcess(j)
        xx0=1000. ! Showers in ZHAireS start at 1000 g cm^-2
        xx(j)=(xx(j)-xx0)/rho/100.   ! m
        xntot=xntot+yexcess(j)
      end do
      close (unit=14)

      xntot=xntot*(xx(2)-xx(1))   ! m
      npoints=j-1
      firstgo=.false.
      zend=xx(npoints)
      write(*,*) 'Depth max (m) ',zend

c      do j=1,npoints
c        if (yexcess(j).gt.0.) then 
c          ylog(j)=log(yexcess(j))
c        else 
c          ylog(j)=0.  
c        end if
c      end do

      end if

c      xnep=ddivdif(ylog,xx,npoints,z,nnint)
      xnep=ddivdif(yexcess,xx,npoints,z,nnint)

c      xnep=-exp(xnep)    ! Excess negative charge => - sign

      RETURN      
      END       

C ******************************************************************
C ******************************************************************
C ******************************************************************
C ******************************************************************
