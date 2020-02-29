!     filter following R.Vuilleumier
! "One often requires only the frequencies of the peaks of the VDOS and not the actual power spectrum.
!If this is the case, it is useful to calculate only the real (cosine) part of the Fourier transform since the peaks will be more well defined by these data."   t im Teatro.
! We follow Tim Teatro's idea, and only calculate the real part of the Fourier Transform.
!     1/03/2004

      implicit none
      integer nmax,i,ii,imax,j,ncor,num_tf,imax2
      parameter(nmax=300000)
      real*8 corr(nmax),corr2(nmax)
      real*8 tf(0:nmax),tf2(0:nmax)
      real*8 wmax,dw,dt,t,omega
      real*8 filter,cel,pi,filter2,filter3
      real*8 convert
      parameter(convert=0.52917725d-10)
      real*8 planck,kt,volume,factor,argumt,denom,class,sigma,num
      real*8 tyty, time


      open(10,file='correl.dat')
      do i=1,nmax
         read(10,*,end=100)ii,corr(i)
      enddo
 100  continue
      imax=i-1
      write(6,*)'imax=',imax
      close(10)
! normalisation pour VDOS
!      do i=2,imax
!         corr(i)=corr(i)/corr(1)
!      enddo
!      corr(1)=1.d0

      do i=1,nmax
         tf(i)=0.d0
         tf2(i)=0.d0
      enddo

      wmax=4000.d0              !cm-1
      dw=1.d0                   !cm-1
      num_tf=wmax/dw
      write(6,*)'num_tf=',num_tf
      if(wmax/dw >= nmax)then
         write(6,*)'!!! dimension table tf !!!'
         stop
      endif

!!      dt = 4.d0*0.024d0*10.d0*1.d-15           !s
! dt is the time-step in seconds (16.6*0.024 = 0.4 fs)
!      dt = 0.4*1.d-15 !s
!      dt = 16.6d0*0.024*1.d-15 !s
      dt = 0.5*1.d-15 !s

      cel = 3.d10 !cm/s
      pi=dacos(-1.d0)
      imax2=imax/2

! tf = real part of the Fourier Transform
! tf2 = imaginary part of the Fourier Transform
      do j=0,num_tf
         omega=2.d0*pi*cel*dw*dble(j)

         do i=1,imax
            t=dt*dble(i)
            tf(j)=tf(j)+dcos(omega*t)*corr(i)*filter(i,imax2)
            tf2(j)=tf2(j)+dsin(omega*t)*corr(i)*filter(i,imax2)
         enddo
         tf(j)=tf(j)*2.d0 !because of DFT for symetrical function (see Allen p.336)
      enddo




! ===============================================      
      open(10,file='spectra.dat')

      do j=0,num_tf

         sigma=dble(j)*dw       !cm-1
         write(10,*)sigma,tf(j)-tf(wmax/dw)

      enddo
      close(10)


      STOP
      END



      real*8 function filter(ind,ncor)
      implicit none
      integer ind,ncor
      
      filter =dexp(-0.5d0*20.d0*(dble(ind)/dble(ncor))*
     1     (dble(ind)/dble(ncor)))

      return
      end


      real*8 function filter2(ind,ncor)
      implicit none
      integer ind,ncor
      
      filter2 =1.0d0-dble(ind)/dble(ncor)

      return
      end


      real*8 function filter3(ind,ncor)
      implicit none
      integer ind,ncor
      real*8 pi

      pi=dacos(-1.d0)
      filter3=1.0+0.5*cos(pi*dble(ind)/dble(ncor))

      return
      end


