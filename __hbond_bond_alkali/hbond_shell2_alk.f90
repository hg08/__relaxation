!===========================
! 2015/10/27
! huang gang
!===========================
      program hbond_shell2_alk
      implicit none
!========================
!parameters and variables
!========================
      character(LEN=200) :: filename            ! specific filename to analyz data
      character(LEN=200) :: pos_filename        ! specific trajectory filename to analyz data
      character(LEN=200) :: name_slide          ! to define filename for each slide  
      character(LEN=200) :: x1                  ! to define filename for each slide  
      integer,parameter :: rk=4              
      real,parameter :: rooc=12.25             ! cutoff distance of rOO (5.4**2 )
      real,parameter :: rooc_min=1.0          ! O atom itself can not be included.
      real,parameter :: rohc=6.0025            ! rOH (3.5**2)
      real,parameter :: r_ohc=1.44             ! rOH (1.2**2)
      real           :: r,r23,r13,r_shell2
      integer    :: begin_time,end_time,rat,i,j,nmovie,kmovie,natoms,&
                    iatom,m1,m2,m3,i_H,i_alk,i1,i2,ii,i4,&
                    i_OW,j_slide,len_piece
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:)        :: ndx_OW,ndx_H,&
          ndx_alk,kmo,i_OW_shell2_alk
      integer,allocatable,dimension (:,:)      :: ndx_OW_shell2_alk

      character(len=8) :: fmt ! format descriptor

     

      call system_clock(begin_time,rat) 
!==================
!read data in input
!==================
      write(6,*)'What is the name of the system:'
      read(5,*)filename
      write(6,*)'What is the name of the trajectory file:'
      read(5,*)pos_filename
      write(6,*)'What is the total number of movie steps:'
      read(5,*)nmovie!number of atoms per molecules
      write(6,*)'What is the number of atoms in the system:'
      read(5,*)natoms!number of atoms per molecules
      write(6,*)'What is the lenght of each piece (Unit:steps):'
      read(5,*)len_piece!number of atoms per molecules
      write(6,*)'What is the distance of r_O-alk(1st and 2nd shell):'
      read(5,*)r_shell2

      kmovie=int(nmovie/len_piece)
      allocate(kmo(kmovie))
      do i = 1,kmovie
          kmo(i)=len_piece*i!number of movie steps:
      enddo
      r_shell2=r_shell2**2
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      allocate(i_OW_shell2_alk(kmovie))
!=======================
!read in trajectory file 
!=======================
      open(10,file=trim(pos_filename))     
      do j_slide=1, kmovie
          ! j_slide=(imovie+len_piece)/len_piece
          !=============================================
          !Skip len_piece*(natoms+2)/2 -(natoms+2) lines
          !=============================================
          do j=1,len_piece*(natoms+2)/2 -(natoms+2)
             read(10,*)
          enddo
          !=======
          !Reading
          !=======
          read(10,*)!Neglect data of this line
          read(10,*)!Neglect data of this line
          i=0
          ii=0
          i4=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,j_slide),& 
                       y(iatom,j_slide),z(iatom,j_slide)
              if (trim(atom_type(iatom)) .eq. 'O') then
                  i=i+1
              elseif(trim(atom_type(iatom)) .eq. 'H') then
                  ii=ii+1
              elseif (trim(atom_type(iatom)) .eq. 'Li' &
                 .or. trim(atom_type(iatom)) .eq. 'Na'&
                 .or. trim(atom_type(iatom)) .eq. 'K' ) then
                  i4=i4+1
              endif
          enddo
          write(6,*) x(iatom-1,j_slide)!FOR TESTING 
          do j=1,len_piece*(natoms+2)/(2*10)
            read(10,*)
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

      i_OW=i
      i_H=ii
      i_alk=i4
      write(6,*)'Number of Ow','    Number of alkali metal ions'
      write(6,*)i_OW,i_alk
      close(10)
!==================================
!Extracting ndx of alk cations
!and O (H) atoms in water molecules
!==================================
      allocate(ndx_OW(i_OW))! this should be put after i_OW is defined
      allocate(ndx_H(i_H))! this should be put after i_H is defined
      allocate(ndx_alk(i_alk))
      i=0
      ii=0
      i4=0
      do iatom=1,natoms
          if (trim(atom_type(iatom)) .eq. 'O')then
                 i=i+1      
                 ndx_OW(i)=iatom
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ii=ii+1
                 ndx_H(ii)=iatom
          elseif (trim(atom_type(iatom)) .eq. 'Li' &
                 .or. trim(atom_type(iatom)) .eq. 'Na'&
                 .or. trim(atom_type(iatom)) .eq. 'K' ) then
                 i4=i4+1
                 ndx_alk(i4)=iatom
          endif
      enddo 
     deallocate(atom_type)
!========================================
!Print index of O in the shell of I- ions
!========================================
      do j_slide=1,kmovie
      i=0
      do i1=1, i_OW! No O atom can not be bonded to itself 
          m1=ndx_OW(i1)
          do i2=1, i_alk 
              m2=ndx_alk(i2)
                  do j =j_slide,j_slide
                      r= (x(m2,j)-x(m1,j))**2+  &
                           (y(m2,j)-y(m1,j))**2+  &
                           (z(m2,j)-z(m1,j))**2
                  if (r<r_shell2) then
                    i=i+1 
                  endif
                  enddo
          enddo
      enddo
      i_OW_shell2_alk(j_slide)=i
      enddo
      allocate(ndx_OW_shell2_alk(MAXVAL(i_OW_shell2_alk),j_slide)) 
      ndx_OW_shell2_alk(:,:)=0 
      
      do j_slide=1,kmovie
      i=0
      do i1=1, i_OW
          m1=ndx_OW(i1)
          do i2=1, i_alk 
              m2=ndx_alk(i2)
              do j =j_slide,j_slide
                  r= (x(m2,j)-x(m1,j))**2+  &
                       (y(m2,j)-y(m1,j))**2+  &
                       (z(m2,j)-z(m1,j))**2
                  if (r<r_shell2) then
                    i=i+1 
                    ndx_OW_shell2_alk(i,j_slide)=m1
                  endif
              enddo
          enddo
      enddo

     !========
     !Printing
     !========

      fmt = '(I4.4)' ! an integer of width 4 with zeros at the left
     !===================================
     !Print the O indics in thei I- shell     
     !===================================
      write (x1,fmt) j_slide ! converting integer to string using a 'internal file'
      name_slide=trim(filename)//trim(x1)
      open(20,file=trim(name_slide)//'_shell2_index.dat')
      do i1=1, i_OW_shell2_alk(j_slide) ! No O atom can not be bonded to itself 
         write(20,fmt='(1I4)',advance='no')ndx_OW_shell2_alk(i1,j_slide)
      enddo
      close(20) 
     !=====================================
     !Print the O-O-H list in thei I- shell     
     !=====================================
      open(20,file=trim(name_slide)//'_shell2_list.dat')
      do i1=1, i_OW_shell2_alk(j_slide)-1! No O atom can not be bonded to itself 
          m1=ndx_OW_shell2_alk(i1,j_slide)
          do i2=i1+1,i_OW_shell2_alk(j_slide) 
              m2=ndx_OW_shell2_alk(i2,j_slide)
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j =1,1
                      r23= (x(m2,j)-x(m3,j))**2+  &
                           (y(m2,j)-y(m3,j))**2+  &
                           (z(m2,j)-z(m3,j))**2
                      r13= (x(m1,j)-x(m3,j))**2+  &
                           (y(m1,j)-y(m3,j))**2+  &
                           (z(m1,j)-z(m3,j))**2
                  if (r23<r_ohc .and. m1 .ne. m2) then
                      write(20,fmt='(3I4)')m1,m2,m3
                  endif
                  if (r13<r_ohc .and. m1 .ne. m2) then
                      write(20,fmt='(3I4)')m2,m1,m3
                  endif
                  enddo
              enddo
          enddo
      enddo
      close(20)

      enddo
      deallocate(ndx_OW,ndx_OW_shell2_alk,ndx_alk,&
                 ndx_H,x,y,z,i_OW_shell2_alk)
!==============================================================
      call system_clock(end_time,rat) 
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
