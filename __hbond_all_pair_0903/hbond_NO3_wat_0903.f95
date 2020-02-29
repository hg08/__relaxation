!============
! 2015/10/14
! huang gang
!============

!=========================================      
! For specified k pairs of water molecules，
! to calculate the average and the auto-
! correlation function of the H-Bond 
! population operator:
! <h>,<h_d>,<hh>, <hh_d>;      
! S_HB(t),S_d_HB(t),
! C_HB(t),C_d_HB(t).

     ! input file: input
     ! the first line of input: nmovie
     ! the 2nd line: natoms
     ! the 3rd line: number of pairs of
     ! molecules, or npair 
     ! the 4th-(npair-3)th: the index 
     ! of O,O,H atoms which will be H-
     ! bonded. 
!=========================================
      program hbond_NO3_wat
      implicit none

!========================
!parameters and variables
!========================
      character(LEN=20) :: filename            ! For specific filename to analyzing data
      character(LEN=20) :: pos_filename            ! For specific filename to analyzing data
      integer,parameter :: rk=4              
      integer,parameter :: nmax=775            ! max number of atoms
      integer,parameter :: nmovie_max=2000000   ! max number of movie
      real,parameter :: rooc=12.25             ! cutoff distance of rOO (3.5**2 )
      real,parameter :: rohc=6.0025            ! rOH (2.45**2)
      real,parameter :: r_ohc=1.44           ! rOH (1.2**2)
      real,parameter :: r_ohc_min=0.25         ! rOH lower bound(0.5)
      real,parameter :: r_ON_min=9.0              ! r_ON lower bound (3.0 A)
      real,parameter :: r_ON_max=20.25           ! r_ON upper bound (4.5 A)
      real           :: r_ON,r23
      !real,parameter :: cosphic=0.866          ! 1.732/2; phiC=pi/6.
      !real(kind=rk),parameter :: delta_t=0.005 ! fs
      integer    :: begin_time,end_time,rat,i,j,nmovie,natoms,iatom,& 
                    imovie,m1,m2,m3,m4,i_O,i_H,i_N,&
                    i1,i2,ii,jj,i4,i_OW,i_O_Nitrate
      !real(kind=rk),allocatable,dimension (:,:,:,:) :: r23 !r12, r13, r23
      ! real(kind=rk),allocatable,dimension (:,:)     :: r_ON 
      ! real(kind=rk),allocatable,dimension (:,:,:,:)  :: cosphi, pm
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:)          :: ndx_O,ndx_OW,&
                                                    ndx_O_Nitrate,&  
                                                    ndx_H,ndx_N
      call system_clock(begin_time,rat) 

!==================
!read data in input
!==================
      !open(10,file='input_hbond_all_pair_0903')
      write(6,*)'What is the name of the syetem:' 
      read(5,*)filename   ! system name
      write(6,*)'What is the name of the trajectory file:' 
      read(5,*)pos_filename
      write(6,*)'What is the total number of steps:' 
      read(5,*)nmovie     !number of steps
      write(6,*)'What is the total number of atoms in the system:' 
      read(5,*)natoms     !number of atoms per molecules
      
      if (natoms >= nmax)then
        write(6,*)'!! natoms> nmax !!' 
        write(6,*)'!! stop !!'
        stop
      endif  

      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))

!=======================
!read in trajectory file 
!=======================
      open(10,file=trim(pos_filename))     
      do imovie=1,nmovie
         read(10,*)                  !Neglect data of this line
         read(10,*)                  !Neglect data of this line
         i=0
         ii=0
         jj=0
         do iatom= 1,natoms
            read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
            if (trim(atom_type(iatom)) .eq. 'O') then
                  i=i+1
            elseif(trim(atom_type(iatom)) .eq. 'H') then
                  ii=ii+1
            elseif(trim(atom_type(iatom)) .eq. 'N') then
                  jj=jj+1
            endif
         enddo
      enddo
      i_O=i
      i_H=ii
      i_N=jj
      !write(6,*)i_O,i_H,i_N
      close(10)
!=======================
!read in trajectory file 
!=======================

!      open(10,file='traj_pos.xyz')     
!      do imovie=1,nmovie
!         read(10,*)                  !Neglect data of this line
!         read(10,*)                  !Neglect data of this line
!         do iatom= 1,natoms
!            read (10,*)atom_type(iatom),x(iatom,imovie),& 
!                       y(iatom,imovie),z(iatom,imovie)
!         enddo
!      enddo
!      close(10)
!      write(6,*)'end of trajectory reading'
!============================

      allocate(ndx_O(i_O))    ! this should be put after i_O is defined
      allocate(ndx_H(i_H))    ! this should be put after i_H is defined
      allocate(ndx_N(i_N))    ! this should be put after i_H is defined
      i=0
      ii=0
      jj=0
      do iatom=1,natoms
          if (trim(atom_type(iatom)) .eq. 'N') then
                 jj=jj+1
                 ndx_N(jj)=iatom
          elseif (trim(atom_type(iatom)) .eq. 'O')then
                 i=i+1      
                 ndx_O(i)=iatom
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ii=ii+1
                 ndx_H(ii)=iatom
          else
          endif
      enddo 
     !TEST
     !write(6,*) jj
     !do i=1, i_N
     !   write(6,*) ndx_N(i)
     !enddo
      
      deallocate(atom_type)
!==================================================
! calculate the number of O atom in water molecules
!==================================================
      i_OW=0
      i_O_Nitrate=0
      do i1=1,i_O
          m1=ndx_O(i1)
          do i4=1,i_N
              m4=ndx_N(i4)
              do j=20000,20000
                  r_ON=((x(m4,j)-x(m1,j))**2+   &
                               (y(m4,j)-y(m1,j))**2+    &    
                               (z(m4,j)-z(m1,j))**2 )
                  if (r_ON>r_ON_min .and. r_ON <r_ON_max) then
                      !write(6,*)m4,m1,j
                      i_OW=i_OW+1
                  else
                      i_O_Nitrate=i_O_Nitrate+1    
                  endif
              enddo
          enddo
       enddo
      allocate(ndx_OW(i_OW))
      allocate(ndx_O_Nitrate(i_O_Nitrate))
      !write(6,*) i_OW

      ! index of O in water molecules
      ! index of O in NO3- ions
      i=0
      ii=0
      do i1=1,i_O
          m1=ndx_O(i1)
          do i4=1,i_N
              m4=ndx_N(i4)
              do j=20000,20000
                 r_ON=((x(m4,j)-x(m1,j))**2+ &
                             (y(m4,j)-y(m1,j))**2+    &    
                             (z(m4,j)-z(m1,j))**2 )
                 if (r_ON>r_ON_min .and. r_ON < r_ON_max) then
                     i=i+1
                     ndx_OW(i)=m1
                 else 
                     ii=ii+1
                     ndx_O_Nitrate(ii)=m1
                 endif
              enddo
          enddo
       enddo
!===========================key part====================================
!allocate(r23(1,i_OW,i_OW,i_H))
      open(20,file=trim(filename)//'_OW-OW-HW_list.dat')
          
      do i1=1, i_OW-1      ! No O atom can not be bonded to itself 
          m1=ndx_OW(i1)
          do i2=i1+1, i_OW 
              m2=ndx_OW(i2)
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j =20000, 20000
                  r23= (x(m2,j)-x(m3,j))**2+  &
                       (y(m2,j)-y(m3,j))**2+  &
                       (z(m2,j)-z(m3,j))**2
                       ! write(20,*) m2
                      if (r23<r_ohc .and.    &
                            r23>r_ohc_min ) then
                            write(20,*) m1,m2,m3
                      endif
                      if ((x(m1,j)-x(m3,j))**2+  &  ! Actually, we do not need to define the array r23
                         (y(m1,j)-y(m3,j))**2+  &
                         (z(m1,j)-z(m3,j))**2 <r_ohc &
                         .and.                       & 
                         (x(m1,j)-x(m3,j))**2+  &
                         (y(m1,j)-y(m3,j))**2+  &
                         (z(m1,j)-z(m3,j))**2 >r_ohc_min)then
                      write(20,*) m2,m1,m3
                      endif
                  enddo
              enddo
          enddo
      enddo
      close(20)
!--------------------------------
!print the list NitrateO-OH pairs
!--------------------------------
      open(20,file=trim(filename)//'_NitrateO-OW-HW_list.dat')
          
          do i1=1, i_O_Nitrate!i_O_Nitrate: Nr. of O atom in nitrate ions 
              m1=ndx_O_Nitrate(i1)
              do i2=i1+1, i_OW 
                  m2=ndx_OW(i2)
                  do ii=1,i_H
                      m3=ndx_H(ii)
                      do j =20000, 20000
                          r23= (x(m2,j)-x(m3,j))**2+  &
                               (y(m2,j)-y(m3,j))**2+  &
                               (z(m2,j)-z(m3,j))**2
                          ! write(20,*) m2
                          if (r23<r_ohc .and.    &
                              r23>r_ohc_min ) then
                              write(20,*) m1,m2,m3
                          endif
                      enddo
                  enddo
              enddo
          enddo
      close(20)

      deallocate(ndx_O,ndx_OW,ndx_O_Nitrate,ndx_H,&
                x,y,z)
!==============================================================
      call system_clock(end_time,rat) 
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
