      program formatBN
      implicit real*8 (a-h,o-z)
      parameter (N=450)
      character linek*120,t3*3,dum*1
      character*1 element(N)
      real*8 x(N),y(N),z(N),vx(N),vy(N),vz(N)
      call getarg(1,t3)
      read(t3(:3),11)it
11    format(i3)
      open(1,file='geo_stp.xyz')
      open(2,file='bnnt.gen')
      open(3,file='veloc.dat')
      cmx=0
      cmy=0
      cmz=0
      do i=0,1000000
      read(1,22,end=99)linek
22    format(a120)
      if(linek.eq."MD iter: 4000")then
      k=0
      do j=1,10000
      read(1,*,end=99)dum,xx,yy,zz,vxx,vyy,vzz
      r=dsqrt(xx**2+yy**2+zz**2)
      if(r.lt.30)then
      k=k+1
      element(k)=dum
      x(k)=xx
      y(k)=yy
      z(k)=zz
      vx(k)=vxx
      vy(k)=vyy
      vz(k)=vzz
      end if
      end do
      end if
      end do
99    continue
      n_atom=k
      amm=0
      do i=1,n_atom
      if(element(i).eq."B")am=10.811
      if(element(i).eq."N")am=14.007
      cmx=cmx+am*x(i)
      cmy=cmy+am*y(i)
      cmz=cmz+am*z(i)
      amm=amm+am
      end do
      cmx=cmx/amm
      cmy=cmy/amm
      cmz=cmz/amm
      do i=1,n_atom
      x(i)=x(i)-cmx 
      y(i)=y(i)-cmy 
      z(i)=z(i)-cmz 
      end do
      PI=4.D0*datan(1.D0)
      theta=(mod(it,11)+1)*PI/12.d0
      phi=(mod(it,5))*2*PI/5.d0
      theta2=(mod(it,7)+1)*PI/8.d0
      phi2=(mod(it,13))*2*PI/13.d0
      vel_BN=3.1702e-1*dsqrt(2000.d0)
      write(2,3)n_atom+2
3     format(i3," C")
      write(2,4)
4     format("B  N")
      do i=1,n_atom
      if(element(i).eq."B")kuk=1
      if(element(i).eq."N")kuk=2
      write(2,5)i,kuk,x(i),y(i),z(i)
5     format(i3,2x,i1,2x,3f12.6)
      write(3,6)vx(i),vy(i),vz(i)
6     format(3f12.6)
      end do
      write(2,7)n_atom+1,"  1  ",25.d0*dsin(theta)*dcos(phi),25.d0*
     1dsin(theta)*dsin(phi),25.d0*dcos(phi)
      write(2,7)n_atom+2,"  2  ",25.d0*dsin(theta)*dcos(phi)+1.287*
     1dsin(theta2)*dcos(phi2),25.d0*dsin(theta)*dsin(phi)+1.287*
     1dsin(theta2)*dsin(phi2),25.d0*dcos(phi)+1.287*dcos(theta2)
7     format(i3,a5,3f12.6)
      write(3,6)-vel_BN*dsin(theta)*dcos(phi),-vel_BN*dsin(theta)*
     1dsin(phi),-vel_BN*dcos(theta)
      write(3,6)-vel_BN*dsin(theta)*dcos(phi),-vel_BN*dsin(theta)*
     1dsin(phi),-vel_BN*dcos(theta)
      end
