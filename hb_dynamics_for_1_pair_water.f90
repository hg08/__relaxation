!     2013/4/10
!     2015/08/27
!     对指定的一对原子，求其氢键布居算符之平均值和其自关联函数。
!=============================================================

      program hbond
      implicit none
!==========
!parameters
!==========      
      real :: hb
      integer,parameter :: nmax=775 !max number of atoms
      integer,parameter :: nmovie_max=30000 ! max number of movie
      real,parameter :: rooc=29.16 !cutoff distance of rOO (5.4**2 )
      real,parameter :: rohc=12.25 !...................rOH (3.5**2)
      real,parameter :: cosphic=0.866  ! 1.732/2; phiC=pi/6.
      integer        ::  j,k,nmovie,natoms,iatom,imovie 
!      integer,dimension(100) ::coun,accum  !此处，能用指针将coun 之维度指示到ndmax吗？
!      real, dimension (100) :: rcoun,raccum  
      real(kind=8),allocatable,dimension (:) :: h, r12, r13, r23
      real(kind=8),allocatable,dimension (:) :: cosphi, pm
!     便于将coun 数组中的整形赋与rcoun,以可计算概率，否则整型数据相除，结果为0.
!     accum(i)代表累积次数之值
      real,allocatable,dimension (:,:) :: x,y,z
      character(LEN=3),allocatable,dimension (:):: atom_type

!=============== 
!data in BOX_VEL
!=============== 
      open(10,file='BOX_VEL')
      read(10,*)nmovie   !number of movie steps
      read(10,*)natoms   !number of atoms per molecules
        
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
      write(6,*)
      close(10)

      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))

      allocate(h(nmovie))
      allocate(cosphi(nmovie))
      allocate(pm(nmovie))
      allocate(r12(nmovie))
      allocate(r13(nmovie))
      allocate(r23(nmovie))
!======================== 
! read in trajectory file 
!======================== 
      open(10,file='traj_pos.xyz')     
      do imovie=1,nmovie
         read(10,*) 
         read(10,*)
         do iatom= 1,natoms
         read (10,*)atom_type(iatom), x(iatom,imovie),y(iatom,imovie),&
                    z(iatom,imovie)
         enddo

      enddo
      close(10)
      write(6,*)'end of trajectory reading'
!===该部分为关键，尤其do ii循环=========
! 1: 2 
! 2: 162
! 3: 575
!===该部分为关键，尤其do ii循环=========
      hb=0
      do j =1, nmovie
         h(j)=0 
         r13(j)= (x(4,j)-x(485,j))*(x(4,j)-x(485,j))+    &  
                (y(4,j)-y(485,j))*(y(4,j)-y(485,j))+          & 
                (z(4,j)-z(485,j))*(z(4,j)-z(485,j))
         r12(j)= (x(4,j)-x(117,j))**2+            &
                 (y(4,j)-y(117,j))**2+                 &
                 (z(4,j)-z(117,j))**2
         r23(j)= (x(117,j)-x(485,j))**2+        &
                (y(117,j)-y(485,j))**2+              &
                (z(117,j)-z(485,j))**2
         pm(j)= (x(485,j)-x(117,j))*(x(4,j)-x(117,j))+     &   
                (y(485,j)-y(117,j))*(y(4,j)-y(117,j))+     & 
                (z(485,j)-z(117,j))*(z(4,j)-z(117,j)) 
         cosphi(j)= pm(j)/(sqrt(r23(j)*r12(j)))            !pm: point multiplication.
         if (r13(j) .lt. rohc .and. r12(j).lt.rooc         & 
            .and. cosphi(j).gt. cosphic) then    
          h(j)=1.0 
          hb=hb+1                          
         endif
      enddo   
      hb=hb/nmovie 
      deallocate (x,y,z,atom_type)
!============
!write result
!============
      open(10,file='nhb.dat')
      write (10,*)'# Step  ','Population  '

         write(10,*)'# Average HBonding probability:', hb
      do k =1,nmovie
         write(10,*)k, h(k)
      enddo
      write(6,*)'written in nhb_profile'
      close(10)
      deallocate (h,cosphi,pm,r12,r13,r23)
!================================================
      END

