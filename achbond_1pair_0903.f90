!     2015/09/03
!=====================================
!对指定的 1 对原子，求其氢键布居算符之
!平均值和其自关联函数S_HB(t),S_d_HB(t),
!C_HB(t),C_d_HB(t)。
     ! input file: input
     ! the first line of input: nmovie
     ! the 2nd line: natoms
     ! the 3rd line: number of pairs of
     ! molecules, or npair 
     ! the 4th-(npair-3)th: the index of O,O,H
     ! atoms which will be H-bonded 
!=====================================
      
      program achbond_1pair
      implicit none
!==========
!parameters
!==========
      integer,parameter :: rk=4              
      integer,parameter :: nmax=775            ! max number of atoms
      integer,parameter :: nmovie_max=300000   ! max number of movie
      real,parameter :: rooc=29.16             ! cutoff distance of rOO (5.4**2 )
      real,parameter :: rohc=12.25             ! rOH (3.5**2)
      real,parameter :: hb_min=0.0000001       ! the condition that a pair is H-bonded 
      real,parameter :: cosphic=0.866          ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: delta_t=0.005 ! fs
      real(kind=rk)           :: r12,r13,r23,cosphi,pm
      integer        :: begin_time,end_time,rat,&
                        i,j,k,nmovie,natoms,iatom,& 
                        imovie,npair,m1,m2,m3,mt 
      real(kind=rk),allocatable,dimension (:,:) :: h,h_d,hh,hh_d
      real,allocatable,dimension (:,:)         :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:)   :: ndx_1, ndx_2, ndx_3
      real,allocatable,dimension (:)      :: hb,hb_d             !for calcualate the average HB population
      real(kind=rk),allocatable,dimension (:) :: corr_h, corr_h_d
      real(kind=rk) :: scalar 
      call system_clock(begin_time,rat)
!==================
!read data in input
!==================
      open(10,file='input')
      read(10,*)nmovie     !number of movie steps
      read(10,*)natoms     !number of atoms per molecules
      read(10,*)npair      !number of pair of water molecules to consider  
      if (nmovie >= nmovie_max)then
        write(6,*)'!! nmovie> nmovie_max !!'
        write(6,*)'!! stop !!'
        stop
      endif
      
      if (natoms >= nmax)then
        write(6,*)'!! natoms> nmax !!' 
        write(6,*)'!! stop !!'
        stop
      endif  
     ! write(6,*)'number of movie steps:',nmovie
     ! write(6,*)'number of atoms in molecule:',natoms
      allocate(ndx_1(npair))
      allocate(ndx_2(npair))
      allocate(ndx_3(npair))
      do k=1,npair
       read(10,*)ndx_1(k),ndx_2(k),ndx_3(k)
      enddo
      close(10)

      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      allocate(h(npair,nmovie))
      allocate(h_d(npair,nmovie))
      allocate(hh(npair,nmovie))
      allocate(hh_d(npair,nmovie))
      allocate(hb(npair))         !Average H-bonded population 
      allocate(hb_d(npair))

!=======================
!read in trajectory file 
!=======================
      open(10,file='traj_pos.xyz')     
      do imovie=1,nmovie
         read(10,*)                 !Neglect data of this line
         read(10,*)                 !Neglect data of this line
         do iatom= 1,natoms
         read (10,*)atom_type(iatom),x(iatom,imovie),y(iatom,imovie),&
                    z(iatom,imovie)
         enddo

      enddo
      close(10)
      write(6,*)'end of trajectory reading'
!===============================关键部分===============================
      write(6,*)'For test:','ndx_1(1)=',ndx_1(1),', ndx_2(1)=',ndx_2(1)
      do k=1,npair
        hb(k)=0
        hb_d(k)=0
        m1=ndx_1(k)
        m2=ndx_2(k)
        m3=ndx_3(k)

        do j =1, nmovie
          h(k,j)=0
          h_d(k,j)=0 
          r13= (x(m1,j)-x(m3,j))**2+  &
                    (y(m1,j)-y(m3,j))**2+  &
                    (z(m1,j)-z(m3,j))**2      !r:squra of distances
          r12= (x(m1,j)-x(m2,j))**2+  &
                    (y(m1,j)-y(m2,j))**2+  &
                    (z(m1,j)-z(m2,j))**2
          r23= (x(m2,j)-x(m3,j))**2+  &
                    (y(m2,j)-y(m3,j))**2+  &
                    (z(m2,j)-z(m3,j))**2
          pm= (x(m3,j)-x(m2,j))*   &
                   (x(m1,j)-x(m2,j))+   & 
                   (y(m3,j)-y(m2,j))*   & 
                   (y(m1,j)-y(m2,j))+   & 
                   (z(m3,j)-z(m2,j))*   &
                   (z(m1,j)-z(m2,j)) 
          cosphi= pm/(sqrt(r23*r12))            !pm: point multiplication.
          if (r13 .lt. rohc .and. r12.lt.rooc          & 
             .and. cosphi.gt. cosphic) then    
              h(k,j)=1.0 
              hb(k)=hb(k)+h(k,j)                          
          endif
          if (r13 .lt. rohc ) then                           
              h_d(k,j)=1.0
              hb_d(k)=hb_d(k)+h_d(k,j)                          
          endif
       enddo   
       hb(k)=hb(k)/nmovie
       hb_d(k)=hb_d(k)/nmovie
      enddo
      
      deallocate (x,y,z,atom_type)
!================
!Write the result
!================
      open(10,file='achbond_1pair_0903_nhb.dat')

      do k=1,npair
        if (hb(k)> hb_min)then
         write(10,*)'#pair(k) ','Ave HBonded population(h_k) ','hb_d_k'
         write(10,*) k, hb(k),hb_d(k)
         write (10,*)'#Step  ','time (ps)','Population(h_k(t))  ',&
                   'h_d_k(t)'
         do j =1,nmovie
            write(10,*)j,j*delta_t,h(k,j),h_d(k,j) 
         enddo
       endif 
      enddo

      write(6,*)'written in nhb_profile'
      close(10)
      deallocate (ndx_1,ndx_2,ndx_3)

!==================================
!Calculate autocorrelation function
!==================================
      allocate(corr_h(nmovie))
      allocate(corr_h_d(nmovie))

! calculate <h(0)h(t)>/<h> 
      ! Notice here <> is not average over
      ! different pairs of water molecules,
      ! but average over the time steps.
      do i=1, nmovie
        corr_h(i)=0
      enddo
      
      do mt=0,nmovie-1     ! time interval
          scalar=0.d0
          do j=1, nmovie-mt-1
              scalar=scalar+h(1,j)*h(1,j+mt)  ! 1: the first pair of water molecules
          enddo
          scalar=scalar/(nmovie-mt)
          corr_h(mt+1)=scalar/hb(1)
      enddo    
  
! calculate <h^d(0)h^d(t)>/<h^d> 
      ! Notice here <> is not average over
      ! different pairs of water molecules,
      ! but average over the time steps.
      do i=1, nmovie
        corr_h_d(i)=0
      enddo
      
      do mt=0,nmovie-1     ! time interval
          scalar=0.d0
          do j=1, nmovie-mt-1
              scalar=scalar+h_d(1,j)*h_d(1,j+mt)  ! 1: the first pair of water molecules
          enddo
          scalar=scalar/(nmovie-mt)
          corr_h_d(mt+1)=scalar/hb_d(1)
      enddo    
  
!======================
!Write the correlation
!C_HB(t) and C_HB_d(t)     
!======================
      open(10,file='achbond_1pair_0903_acf_h.dat')
        do i=1,int(nmovie*2/5)
            write(10,*)i-1,corr_h(i)
        enddo
        write(6,*)'written in achbond_1pair_acf_h.dat'
      close(10)

      open(10,file='achbond_1pair_0903_acf_h_d.dat')
        do i=1,int(nmovie*2/5)
            write(10,*)i-1,corr_h_d(i)
        enddo
        write(6,*)'written in achbond_1pair_0903_acf_h_d.dat'
      close(10)

!==========================
!Write the correlation
!ln(C_HB(t)) and lnC_HB_d(t)     
!==========================
      open(10,file='achbond_1pair_0903_ln_acf_h.dat')
        do i=1,int(nmovie*2/5)
            write(10,*)i-1,log(corr_h(i))
        enddo
        write(6,*)'written in achbond_1pair_0903_ln_acf.dat'
      close(10)

      open(10,file='achbond_1pair_0903_ln_acf_h_d.dat')
        do i=1,int(nmovie*2/5)
            write(10,*)i-1,log(corr_h_d(i))
        enddo
        write(6,*)'written in achbond_1pair_0903_ln_acf_h_d.dat'
      close(10)

      deallocate (h,h_d,hh,hh_d,hb,hb_d,corr_h,corr_h_d)
!==============================================================
      call system_clock(end_time,rat)
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END

