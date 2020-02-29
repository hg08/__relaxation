!===========================
! 2015/11/19
! Author: Huang Gang
!===========================
!======================================================================
! In this function, we use Nitrate-N-OW distance to determine the water
! molecules in the NO3--shell. (To overcome repeated indices PROBLEM!)      
!======================================================================
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
      character(LEN=200) :: name_slide           ! To define filename for each traj. piece
      character(LEN=200) :: x1                   ! To define filename for each traj. piece
      integer,parameter :: rk=4              
      real(kind=rk) :: r_ohc                         ! rOH (1.2 A)
      real(kind=rk) :: r_ohc_min                     ! rOH lower bound(0.5 A)
      real(kind=rk) :: r_ON_min                      ! r_ON lower bound (3.0 A)
      real(kind=rk) :: r_ON_max                      ! r_ON upper bound (4.5 A)
      real(kind=rk) :: r_ON,r23
      integer    :: begin_time,end_time,rat,i,j,nmovie,natoms,iatom,& 
                    m1,m2,m3,m4,i_O,i_H,i_N,len_piece,j_slide,&
                    i1,i2,ii,iii,i4,i_OW,i_O_Nitrate,kmovie
      real(kind=rk),allocatable,dimension (:,:)           :: x,y,z
      integer,allocatable,dimension (:,:) :: ndx_OW_bond_Nitrate,& 
          ndx_OW_bond_wat
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:) :: ndx_O,ndx_OW,&
          ndx_O_Nitrate,ndx_H,ndx_N,i_OW_bond_Nitrate,i_OW_bond_wat,&
          kmo
      character(LEN=8) :: fmt  ! format descriptor
      call system_clock(begin_time,rat) 

!==================
!read data in input
!==================
      write(6,*)'What is the name of the syetem:' 
      read(5,*)filename   ! system name
      write(6,*)'What is the name of the trajectory file:' 
      read(5,*)pos_filename
      write(6,*)'What is the total number of movie steps:' 
      read(5,*)nmovie    
      write(6,*)'What is the total number of atoms in the system:' 
      read(5,*)natoms     !number of atoms per molecules
      write(6,*)'What is the length of each piece (unit: steps):' 
      read(5,*)len_piece    
      write(6,*)'The largest distance of Nit_N-OW distance r_ON_max(A):'
      read(5,*)r_ON_max
      
      kmovie=int(nmovie/len_piece)
      allocate(kmo(kmovie))
      do i =1, kmovie
          kmo(i)=len_piece*i !No. of movie steps
      enddo
      
      r_ohc=1.2
      r_ohc_min=0.5
      r_ON_min=3.0
      r_ohc=r_ohc**2
      r_ohc_min=r_ohc_min**2
      r_ON_min=r_ON_min**2
      r_ON_max=r_ON_max**2

      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      allocate(i_OW_bond_Nitrate(kmovie))
      allocate(i_OW_bond_wat(kmovie))

!=======================
!read in trajectory file 
!=======================
      open(10,file=trim(pos_filename))     

      do j_slide=1,kmovie
           !=========================================
           !skip (len_piece-1)*(natoms+2)*10/5 lines
           !=========================================
           do j=1, len_piece*(natoms+2)/2 - (natoms+2)
               read(10,*)
           enddo
           !=======
           !Reading
           !=======
           read(10,*)
           read(10,*)
           i=0
           ii=0
           iii=0
           do iatom =1, natoms
            read (10,*)atom_type(iatom),x(iatom,j_slide),& 
                       y(iatom,j_slide),z(iatom,j_slide)
            if (trim(atom_type(iatom)) .eq. 'O') then
                  i=i+1
            elseif(trim(atom_type(iatom)) .eq. 'H') then
                  ii=ii+1
            elseif(trim(atom_type(iatom)) .eq. 'N') then
                  iii=iii+1
            endif
         enddo
         write(6,*)x(iatom-1,j_slide), "x, FOR TESTING." ! FOR TESTING
           !=====================================
           !skip len_piece*(natoms+2)*10/5 lines
           !=====================================
           do j=1, len_piece*(natoms+2)/(2*10)
               read(10,*) !there are 20 lines on belowing
               read(10,*)
               read(10,*)
               read(10,*)
               read(10,*)
               read(10,*)
               read(10,*)
               read(10,*)
               read(10,*)
               read(10,*)
           enddo
      enddo
      i_O=i
      i_H=ii
      i_N=iii
      Write(6,*)'No. of O atoms  ','No. of NO3- ions'
      write(6,*)i_O," ",i_N
      close(10) 
!=======================
!Extracting ndx of O,H,N
!=======================
      allocate(ndx_O(i_O))    ! this should be put after i_O is defined
      allocate(ndx_H(i_H))    ! this should be put after i_H is defined
      allocate(ndx_N(i_N))    ! this should be put after i_N is defined
      i=0
      ii=0
      iii=0
      do iatom=1,natoms
          if(trim(atom_type(iatom)) .eq. 'O')then
                 i=i+1      
                 ndx_O(i)=iatom
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ii=ii+1
                 ndx_H(ii)=iatom
          elseif(trim(atom_type(iatom)) .eq. 'N') then
                 iii=iii+1
                 ndx_N(iii)=iatom
          else
          endif
      enddo 
      
      deallocate(atom_type)
!==================================================
! calculate the number of O atom in water molecules
!==================================================

      i_OW=0
      i_O_Nitrate=0
      i=0
      ii=0
      do i1=1,i_O
          m1=ndx_O(i1)
          do i4=1,i_N
              m4=ndx_N(i4)
              do j=1,1
                  r_ON=((x(m4,j)-x(m1,j))**2+   &
                        (y(m4,j)-y(m1,j))**2+   &    
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
               do j=1,1
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
       enddo
!==================================================
! calculate the number of O atom bonded to Nitrate 
!==================================================
      i_OW_bond_Nitrate(:)=0
      i_OW_bond_wat(:)=0

      do j_slide=1,kmovie

      i=0
      ii=0
      do i1=1,i_OW
          m1=ndx_OW(i1)
          do i4=1,i_N
              m4=ndx_N(i4)
              do j=j_slide,j_slide
                 r_ON=((x(m4,j)-x(m1,j))**2+ &
                    (y(m4,j)-y(m1,j))**2+    &    
                    (z(m4,j)-z(m1,j))**2 )
                 if (r_ON < r_ON_max .and.   &
                     r_ON > r_ON_min) then
                     i=i+1
                    ! write(6,*) m4,m1 
                 endif
                 if (r_ON > r_ON_max) then
                     ii=ii+1
                 endif
              enddo
          enddo
      enddo
      i_OW_bond_Nitrate(j_slide)=i ! Use "OW around Nitrate" is better.
      i_OW_bond_wat(j_slide)=ii
      enddo

      allocate(ndx_OW_bond_Nitrate(MAXVAL(i_OW_bond_Nitrate),kmovie)) 
      allocate(ndx_OW_bond_wat(MAXVAL(i_OW_bond_wat),kmovie)) 
      ndx_OW_bond_Nitrate(:,:)=0
      ndx_OW_bond_wat(:,:)=0
!=========================
!=========================
!=========================
      do j_slide=1,kmovie

      i=0
      ii=0
      do i1=1,i_OW
          m1=ndx_OW(i1)
          do i4=1,i_N
              m4=ndx_N(i4)
              do j=j_slide,j_slide
                 r_ON=((x(m4,j)-x(m1,j))**2+ &
                    (y(m4,j)-y(m1,j))**2+    &    
                    (z(m4,j)-z(m1,j))**2 )
                 if (r_ON < r_ON_max .and.   &
                     r_ON > r_ON_min) then
                     i=i+1
                     ndx_OW_bond_Nitrate(i,j_slide)=m1
                 endif
                 if (r_ON > r_ON_max) then
                     ii=ii+1
                     ndx_OW_bond_wat(ii,j_slide)=m1
                 endif
              enddo    
          enddo
       enddo

      !========
      !Printing
      !========
      fmt= '(I4.4)' !an integer of width 4 with zeros at the left 
      !====================================================
      !Print the O indics of O which bonded to Nitrate ions
      !====================================================
      write (x1,fmt) j_slide
      name_slide=trim(filename)//trim(x1)
      
      open(20,file=trim(name_slide)//'_bond_Nitrate_index.dat')
      do i1=1,i_OW_bond_Nitrate(j_slide)
           write(20,fmt='(1I4)',advance='no') &
               ndx_OW_bond_Nitrate(i1,j_slide)
      enddo 
      close(20)
      !=============================
      !print the list OW-OW-OH pairs
      !=============================
      open(20,file=trim(name_slide)//'_OW-OW-HW_list.dat')
          
      do i1=1, i_OW_bond_wat(j_slide)-1      ! No O atom can not be bonded to itself 
          m1=ndx_OW_bond_wat(i1,j_slide)
          do i2=i1+1, i_OW_bond_wat(j_slide) 
              m2=ndx_OW_bond_wat(i2,j_slide)
              if (m2 .ne. m1) then
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j=1,1
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
              enddo
              endif
          enddo
      enddo
      close(20)
      !================================
      !print the list NitrateO-OH pairs
      !================================
      open(20,file=trim(name_slide)//'_NitrateO-OW-HW_list.dat')
          
          do i1=1, i_O_Nitrate!i_O_Nitrate: Nr. of O atom in nitrate ions 
              m1=ndx_O_Nitrate(i1)
              do i2=1, i_OW_bond_Nitrate(j_slide) 
                  m2=ndx_OW_bond_Nitrate(i2,j_slide)
                  do ii=1,i_H
                      m3=ndx_H(ii)
                      do j=1,1
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
          enddo
      close(20)
      
      enddo !j_slide-loop

      deallocate(ndx_O,ndx_OW,ndx_O_Nitrate,ndx_H,&
                ndx_OW_bond_Nitrate,ndx_OW_bond_wat,kmo,x,y,z)
!==============================================================
      call system_clock(end_time,rat) 
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
