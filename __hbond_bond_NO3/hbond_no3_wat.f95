!===========================
! 2015/10/27
! Author: Huang Gang
!===========================
!===========================================================================      
! The function 'hbond_bond_no3.f90' is different from 'hbond_shell2_hal.f90',
! because there are O atoms in NO3- ions, we can print '*_NitrateO-Ow-Hw_list
! .dat' and '*_OW-OW-HW_list.dat'.Therefore, we we split the whole trajectory
! into several pieces of trajectories, all the trajectories share same list
! files.      
!===========================================================================      
      program hbond_no3_wat
      implicit none
!=========================================      
! For specified k pairs of molecules，
! to calculate the average and the auto-
! correlation function of the H-Bond 
! population operator:
! <h>,<h_d>,<hh>, <hh_d>;      
! S_HB(t),S_d_HB(t),
! C_HB(t),C_d_HB(t).
! Here we list all the pairs of molecules 
! which may  form H bonds.      

! input filename 'input': 
     ! 1st line of input: sytem name
     ! 2nd line: posistion filename
     ! 3rd line: number of atoms
!========================
!parameters and variables
!========================
      character(LEN=200) :: filename             ! For specific filename to analyzing data
      character(LEN=200) :: pos_filename         ! For specific filename to analyzing data
      integer,parameter :: rk=4              
      real,parameter :: rooc=5**2                ! cutoff distance of rOO (Angstrom;to extract water molecules in NO3--shell )
      !real,parameter :: rohc=2.45**2            ! rOH (Angstrom)
      real(kind=rk) :: r_ohc                         ! rOH (1.2 A)
      real(kind=rk) :: r_ohc_min                     ! rOH lower bound(0.5 A)
      real(kind=rk) :: r_ON_min                      ! r_ON lower bound (3.0 A)
      real(kind=rk) :: r_ON_max                      ! r_ON upper bound (4.5 A)
      real(kind=rk) :: r_ON,r_OO,r23
      integer    :: begin_time,end_time,rat,i,j,nmovie,natoms,iatom,& 
                    imovie,m1,m2,m3,m4,i_O,i_H,i_N,&
                    i1,i2,ii,iii,i4,i_OW,i_O_Nitrate,i_OW_bond_Nitrate,&
                    i_OW_bond_wat
      real(kind=rk),allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:) :: ndx_O,ndx_OW,&
          ndx_O_Nitrate, ndx_H,ndx_N,ndx_OW_bond_Nitrate,ndx_OW_bond_wat
      call system_clock(begin_time,rat) 

!==================
!read data in input
!==================
      write(6,*)'What is the name of the syetem:' 
      read(5,*)filename   ! system name
      write(6,*)'What is the name of the trajectory file:' 
      read(5,*)pos_filename
      write(6,*)'What is the total number of atoms in the system:' 
      read(5,*)natoms     !number of atoms per molecules

      nmovie=1
      r_ohc=1.2
      r_ohc_min=0.5
      r_ON_min=3.0
      r_ON_max=4.5
      r_ohc=r_ohc**2
      r_ohc_min=r_ohc_min**2
      r_ON_min=r_ON_min**2
      r_ON_max=r_ON_max**2
      
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
         iii=0
         do iatom= 1,natoms
            read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
            if (trim(atom_type(iatom)) .eq. 'O') then
                  i=i+1
            elseif(trim(atom_type(iatom)) .eq. 'H') then
                  ii=ii+1
            elseif(trim(atom_type(iatom)) .eq. 'N') then
                  iii=iii+1
            endif
         enddo
      enddo
      i_O=i
      i_H=ii
      i_N=iii

      close(10)
!============================

      allocate(ndx_O(i_O))    ! this should be put after i_O is defined
      allocate(ndx_H(i_H))    ! this should be put after i_H is defined
      allocate(ndx_N(i_N))    ! this should be put after i_N is defined
      i=0
      ii=0
      iii=0
      do iatom=1,natoms
          if (trim(atom_type(iatom)) .eq. 'N') then
                 iii=iii+1
                 ndx_N(iii)=iatom
          elseif (trim(atom_type(iatom)) .eq. 'O')then
                 i=i+1      
                 ndx_O(i)=iatom
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ii=ii+1
                 ndx_H(ii)=iatom
          else
          endif
      enddo 
      
      deallocate(atom_type)
!==================================================
! calculate the number of O atom in water molecules
!==================================================
      do j=1,1

      i_OW=0
      i_O_Nitrate=0
      i=0
      ii=0
      do i1=1,i_O
          m1=ndx_O(i1)
          do i4=1,i_N
              m4=ndx_N(i4)
                  r_ON=((x(m4,j)-x(m1,j))**2+   &
                        (y(m4,j)-y(m1,j))**2+    &    
                        (z(m4,j)-z(m1,j))**2 )
                  if (r_ON >r_ON_max) then
                      i=i+1
                  endif
                  if (r_ON<r_ON_min) then
                      ii=ii+1    
                      !write(6,*) i_O_Nitrate
                  endif
          enddo
      enddo
      i_OW=i
      i_O_Nitrate=ii

      allocate(ndx_OW(i_OW))
      allocate(ndx_O_Nitrate(i_O_Nitrate))

      i=0
      ii=0
      do i1=1,i_O
          m1=ndx_O(i1)
          do i4=1,i_N
               m4=ndx_N(i4)
                 r_ON=((x(m4,j)-x(m1,j))**2+ &
                             (y(m4,j)-y(m1,j))**2+    &    
                             (z(m4,j)-z(m1,j))**2 )
                 if (r_ON > r_ON_max) then
                     i=i+1
                     ndx_OW(i)=m1
                 endif 
                 if (r_ON<r_ON_min) then
                     ii=ii+1
                     ndx_O_Nitrate(ii)=m1
                 endif
          enddo
       enddo
!==================================================
! calculate the number of O atom bonded to Nitrate 
!==================================================
      i_OW_bond_Nitrate=0
      i_OW_bond_wat=0
      i=0
      ii=0
      do i1=1,i_OW
          m1=ndx_OW(i1)
          do i4=1,i_O_Nitrate
              m4=ndx_O_Nitrate(i4)
                 r_OO=((x(m4,j)-x(m1,j))**2+ &
                    (y(m4,j)-y(m1,j))**2+    &    
                    (z(m4,j)-z(m1,j))**2 )
                 if (r_OO < rooc) then
                     i=i+1
                    ! write(6,*) m4,m1 
                 endif
                 if (r_OO > rooc) then
                     ii=ii+1
                 endif
          enddo
      enddo
      i_OW_bond_Nitrate=i
      i_OW_bond_wat=ii

      allocate(ndx_OW_bond_Nitrate(i_OW_bond_Nitrate)) 
      allocate(ndx_OW_bond_wat(i_OW_bond_wat)) 
      i=0
      ii=0
      do i1=1,i_OW
          m1=ndx_OW(i1)
          do i4=1,i_O_Nitrate
              m4=ndx_O_Nitrate(i4)
                 r_OO=((x(m4,j)-x(m1,j))**2+ &
                    (y(m4,j)-y(m1,j))**2+    &    
                    (z(m4,j)-z(m1,j))**2 )
                 if (r_OO < rooc) then
                     i=i+1
                     ndx_OW_bond_Nitrate(i)=m1
                 endif
                 if (r_OO > rooc) then 
                     ii=ii+1
                     ndx_OW_bond_wat(ii)=m1
                 endif
          enddo
       enddo

       enddo !j-loop
!==============================
!allocate(r23(1,i_OW,i_OW,i_H))
!==============================
      do j =1, 1
      !=============================
      !print the list OW-OW-OH pairs
      !=============================
      open(20,file=trim(filename)//'_OW-OW-HW_list.dat')
          
      do i1=1, i_OW_bond_wat-1      ! No O atom can not be bonded to itself 
          m1=ndx_OW_bond_wat(i1)
          do i2=i1+1, i_OW_bond_wat 
              m2=ndx_OW_bond_wat(i2)
              if (m2 .ne. m1) then
              do ii=1,i_H
                  m3=ndx_H(ii)
                  r23= (x(m2,j)-x(m3,j))**2+  &
                       (y(m2,j)-y(m3,j))**2+  &
                       (z(m2,j)-z(m3,j))**2
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
              endif
          enddo
      enddo
      close(20)
      !================================
      !print the list NitrateO-OH pairs
      !================================
      open(20,file=trim(filename)//'_NitrateO-OW-HW_list.dat')
          
          do i1=1, i_O_Nitrate!i_O_Nitrate: Nr. of O atom in nitrate ions 
              m1=ndx_O_Nitrate(i1)
              do i2=1, i_OW_bond_Nitrate 
                  m2=ndx_OW_bond_Nitrate(i2)
                  do ii=1,i_H
                      m3=ndx_H(ii)
                          r23= (x(m2,j)-x(m3,j))**2+  &
                               (y(m2,j)-y(m3,j))**2+  &
                               (z(m2,j)-z(m3,j))**2
                          if (r23<r_ohc .and.    &
                              r23>r_ohc_min ) then
                              write(20,*) m1,m2,m3
                          endif
                  enddo
              enddo
          enddo
      close(20)
      
      enddo!j-loop

      deallocate(ndx_O,ndx_OW,ndx_O_Nitrate,ndx_H,&
                ndx_OW_bond_Nitrate,ndx_OW_bond_wat,x,y,z)
!==============================================================
      call system_clock(end_time,rat) 
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
