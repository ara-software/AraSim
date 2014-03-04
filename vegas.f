C**************************************************************************
       block data
      implicit double precision(a-h,o-z)
      common /bveg1/ncall,itmx,nprn,ndev,xl(18),xu(18),acc
     1       /bveg2/it,ndo,si,swgt,schi,xi(50,18)
     2       /bveg3/alph,ndmx,mds
      data ncall/5000/,itmx/5/,nprn/5/,acc/-1./,
     1     xl/18*0./,
     3     alph/1.5/,ndmx/50/,mds/1/,ndev/55/,
     4     ndo/1/,xi/900*1./,it/0/,si,swgt,schi/3*0./
      end
C=========================================================================
c the following routines are essential to vegas()

c N.B. Integrals up to dimension 10 can be done with the code
c   as it stands. To work in higher dimensions change every 10
c   in the declarations to whatever maximum dimension you require,
c   and recompile.

      subroutine vegas(ndim,fxn,avgi,sd,chi2a)
c   ndim-dimensional Monte Carlo integration of fxn
c        - gp lepage  sept 1976/(rev)aug 1979  (j comp phys 27,192(1978))

      implicit double precision(a-h,o-z)
      common/bveg1/ncall,itmx,nprn,ndev,xl(18),xu(18),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,18)
      common/bveg3/alph,ndmx,mds
      dimension d(50,18),di(50,18),xin(50),r(50),dx(18),ia(18),
CCC  1          kg(10),dt(10),x(10)
     1          kg(18),dt(18)
      COMMON/XF/XF
      COMMON /X/x(18)
      CHARACTER NAMEOU*12
      COMMON /NAME/NAMEOU
      integer idum
      real ran2
      external fxn
      data idum/-1/
      sqrt(a)=dsqrt(a)
      alog(a)=dlog(a)
      abs(a)=dabs(a)

      ndo=1
      do 1 j=1,ndim
 1    xi(1,j)=1.

      entry vegas1(ndim,fxn,avgi,sd,chi2a)
c         - initializes cummulative variables, but not the grid
      it=0
      si=0.
      swgt=si
      schi=si

      entry vegas2(ndim,fxn,avgi,sd,chi2a)
c         - no initialization
      nd=ndmx
      ng=1
      if(mds.eq.0) go to 2
      ng=(ncall/2.)**(1./ndim)
      mds=1
      if((2*ng-ndmx).lt.0) go to 2
      mds=-1
      npg=ng/ndmx+1
      nd=ng/npg
      ng=npg*nd
 2    k=ng**ndim
      npg=ncall/k
      if(npg.lt.2) npg=2
      calls=npg*k
      dxg=1d0/ng
      dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-1d0)
      xnd=nd
      ndm=nd-1
      dxg=dxg*xnd
      xjac=1./calls
      do 3 j=1,ndim
      dx(j)=xu(j)-xl(j)
 3    xjac=xjac*dx(j)

c   rebin, preserving densities
      if(nd.eq.ndo) go to 8
      rc=ndo/xnd
      do 7 j=1,ndim
      k=0
      xn=0
      dr=xn
      i=k
 4    k=k+1
      dr=dr+1.
      xo=xn
      xn=xi(k,j)
 5    if(rc.gt.dr) go to 4
      i=i+1
      dr=dr-rc
      xin(i)=xn-(xn-xo)*dr
      if(i.lt.ndm) go to 5
      do 6 i=1,ndm
 6    xi(i,j)=xin(i)
 7    xi(nd,j)=1.
      ndo=nd
c     if (it.le.24) go to 1111
 8    if(nprn.ge.0) THEN
      open(ndev,file=nameou,access='append',status='UNKNOWN')
*      write(ndev,200) ndim,calls,it,itmx,acc,nprn,
      write(ndev,*) ndim,calls,it,itmx,acc,nprn,
     1                    alph,mds,nd,(j,xl(j),j,xu(j),j=1,ndim)
c      print      200, ndim,calls,it,itmx,acc,nprn,
c     1                    alph,mds,nd,(j,xl(j),j,xu(j),j=1,ndim)
      close(NDEV)
      ENDIF
1111  continue	

      entry vegas3(ndim,fxn,avgi,sd,chi2a)
c         - main integration loop
 9    it=it+1
      ti=0.
      tsi=ti
      do 10 j=1,ndim
      kg(j)=1
      do 10 i=1,nd
      d(i,j)=ti
 10   di(i,j)=ti

 11   fb=0.
      f2b=fb
      k=0
 12   k=k+1
      wgt=xjac
      do 15 j=1,ndim
      xn=(kg(j)-ran2(idum))*dxg+1.
      ia(j)=xn
      if(ia(j).gt.1) go to 13
      xo=xi(ia(j),j)
      rc=(xn-ia(j))*xo
      go to 14
 13   xo=xi(ia(j),j)-xi(ia(j)-1,j)
      rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
 14   x(j)=xl(j)+rc*dx(j)
 15   wgt=wgt*xo*xnd

      f=wgt
CCC   f=f*fxn(x,wgt)
      f=f*fxn(wgt)
      f2=f*f
      fb=fb+f
      f2b=f2b+f2
      do 16 j=1,ndim
      di(ia(j),j)=di(ia(j),j)+f
 16   if(mds.ge.0) d(ia(j),j)=d(ia(j),j)+f2
      if(k.lt.npg) go to 12

      f2b=sqrt(f2b*npg)
      f2b=(f2b-fb)*(f2b+fb)
      ti=ti+fb
      tsi=tsi+f2b
      if(mds.ge.0) go to 18
      do 17 j=1,ndim
 17   d(ia(j),j)=d(ia(j),j)+f2b
 18   k=ndim
 19   kg(k)=mod(kg(k),ng)+1
      if(kg(k).ne.1) go to 11
      k=k-1
      if(k.gt.0) go to 19

c   compute final results for this iteration
      tsi=tsi*dv2g
      ti2=ti*ti
      wgt=1./tsi
      si=si+ti*wgt
      swgt=swgt+wgt
      schi=schi+ti2*wgt
      avgi=si/swgt
      chi2a=(schi-si*avgi)/(it-.9999)
      sd=sqrt(1./swgt)
      tsi=sqrt(tsi)
c     if (it.le.24) go to 1112
      if(nprn.lt.0) go to 21
      open(ndev,file=nameou,access='append',status='UNKNOWN')
*      write(ndev,201) it,ti,tsi,avgi,sd,chi2a
      write(ndev,*) it,ti,tsi,avgi,sd,chi2a
c      print      201, it,ti,tsi,avgi,sd,chi2a
      close(NDEV)
1112  continue	
      if(nprn.eq.0) go to 21
      do 20 j=1,ndim
c 20   write(ndev,202) j,(xi(i,j),di(i,j),i=1+nprn/2,nd,nprn)
  20  i=1+nprn/2

c   refine grid
 21   do 23 j=1,ndim
      xo=d(1,j)
      xn=d(2,j)
      d(1,j)=(xo+xn)/2.
      dt(j)=d(1,j)
      do 22 i=2,ndm
      d(i,j)=xo+xn
      xo=xn
      xn=d(i+1,j)
      d(i,j)=(d(i,j)+xn)/3.
 22   dt(j)=dt(j)+d(i,j)
      d(nd,j)=(xo+xn)/2.
 23   dt(j)=dt(j)+d(nd,j)

      do 28 j=1,ndim
      rc=0.
      do 24 i=1,nd
      r(i)=0.
      if(d(i,j).le.0.) go to 24
      xo=dt(j)/d(i,j)
      r(i)=((xo-1.)/xo/alog(xo))**alph
 24   rc=rc+r(i)
      rc=rc/xnd
      k=0
      xn=0.
      dr=xn
      i=k
 25   k=k+1
      dr=dr+r(k)
      xo=xn
      xn=xi(k,j)
 26   if(rc.gt.dr) go to 25
      i=i+1
      dr=dr-rc
      xin(i)=xn-(xn-xo)*dr/r(k)
      if(i.lt.ndm) go to 26
      do 27 i=1,ndm
 27   xi(i,j)=xin(i)
 28   xi(nd,j)=1.

      if(it.lt.itmx.and.acc*abs(avgi).lt.sd) go to 9

c 200  format(/35h,'input params for vegas ndim=',i3,8h,'ncall=',f8.0,
c     1  / 28x,5h,' it=',i5,7h,'itmx=',i5,/ 28x,6h,' acc=',g10.3,
c     2  / 28x,7h, ' nprn=',i3,7h,' alph=',f5.2, / 28x,6h,' mds='i3,6h,
c     3  ' nd=',i4,
c     4  /(30x,3h,' xl(',i2,3h,')=', g11.4,4h ' xu(',i2,3h,')=' g11.4))
c 201  format(
c     1  /14h iteration no.i3,1h:9x,10hintegral =g14.7,4h+/- g9.2
c     2  /10h all iter:,
c     3  12h  integral =g14.7,4h+/- g9.2,
c     4    16h   chi**2/it'n =g9.2/)
c 202  format(/15h data for axis i2/25h    x       delta i       ,
c     1  24h   x       delta i       ,18h   x       delta i
c     2  /(1h f7.5,1x,g11.4,5x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4))
      return
      end

      function ran2(idum)
c   random numbers uniformly distributed between 0 and 1.
c   (code by A.T. Service, Harvard-Smithsonian Center for Astrophysics;
c   can be replaced by any other suitable generator for random numbers)
      common/ixtbl/ix1,ix2,ix3,rm1,rm2,r(99)
      data ia1,ic1,m1/1279,351762,1664557/
      data ia2,ic2,m2/2011,221592,1048583/
      data ia3,ic3,m3/15551,6150,29101/
      if(idum.ge.0) go to 2
      ix1=mod(-idum,m1)
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ix1,m2)
      ix1=mod(ia1*ix1+ic1,m1)
      ix3=mod(ix1,m3)
      rm1=1./float(m1)
      rm2=1./float(m2)
      do 1 j=1,99
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
 1    r(j)=(float(ix1)+float(ix2)*rm2)*rm1
 2    ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(99*ix3)/m3
      ran2=r(j)
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      idum=ix1
      return
      end
