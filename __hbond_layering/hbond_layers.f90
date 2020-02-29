!============
! 2015/09/04
! huang gang
!============
!===============================================
      program hbond_interface
      implicit none
!========================
!parameters and variables
!========================
      character(LEN=20) :: filename            ! specific filename to analyz data
      character(LEN=20) :: pos_filename        ! specific trajectory filename to analyz data
      integer,parameter :: rk=4              
      real,parameter :: rooc=12.25             ! cutoff distance of rOO (5.4**2 )
      real,parameter :: rohc=6.0025            ! rOH (3.5**2)
      real,parameter :: r_ohc=1.44             ! rOH (1.2**2)
      real           :: r23,r13,l_1,l_2,l_3,l_4,l_5,l_6,l_7,l_8
      integer    :: begin_time,end_time,rat,i,j,nmovie,natoms,iatom,& 
                    imovie,m1,m2,m3,i_H,i_Li,i_Na,i_K,i_N,i_I,&
                    i_OW_l1,i_OW_l2,i_OW_l3,i_OW_l4,i_OW_l5,i_OW_l6,&
                    i_OW_l7,i_OW_l8,i1,i2,ii,iii,iiii,i5,i6,i7,i_OW
      logical    :: hbonded      
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:)          :: ndx_OW,ndx_H,&
          ndx_OW_l1,ndx_OW_l2,ndx_OW_l3,ndx_OW_l4,ndx_OW_l5,&
          ndx_OW_l6,ndx_OW_l7,ndx_OW_l8,ndx_Li
      call system_clock(begin_time,rat) 

      hbonded=.TRUE.
!==================
!read data in input
!==================
      write(6,*)'What is the name of the system:'
      read(5,*)filename
      write(6,*)'What is the name of the trajectory file:'
      read(5,*)pos_filename
      write(6,*)'What is the number of atoms in the system:'
      read(5,*)natoms!number of atoms per molecules
      write(6,*)'What is normal coordinate of 1st layer in the system:'
      read(5,*)l_1!number of atoms per molecules

      nmovie=1!number of movie steps: 1
      l_2=l_1+1
      l_3=l_1+2
      l_4=l_1+3
      l_5=l_1+4
      l_6=l_1+5
      l_7=l_1+6
      l_8=l_1+7
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
          iiii=0
          i5=0
          i6=0
          i7=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O') then
                  i=i+1
              elseif(trim(atom_type(iatom)) .eq. 'H') then
                  ii=ii+1
              elseif(trim(atom_type(iatom)) .eq. 'Li') then
                  iii=iii+1
              elseif(trim(atom_type(iatom)) .eq. 'Na') then
                  iiii=iiii+1
              elseif(trim(atom_type(iatom)) .eq. 'K') then
                  i5=i5+1
              elseif(trim(atom_type(iatom)) .eq. 'I') then
                  i6=i6+1
              elseif(trim(atom_type(iatom)) .eq. 'N') then
                  i7=i7+1
              endif
          enddo
      enddo
      i_OW=i
      i_H=ii
      i_Li=iii
      i_Na=iiii
      i_K=i5
      i_I=i6
      i_N=i7
      write(6,*)'Number of Ow','    Number of Hw'
      write(6,*)i_OW,i_H
      close(10)
!=======================
!read in trajectory file 
!=======================
      open(10,file=trim(pos_filename))     
      do imovie=1,nmovie
          read(10,*)        
          read(10,*)       
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O' .and.&
                  x(iatom,imovie) .lt. l_1) then
                  i=i+1
              endif
          enddo
      enddo
      i_OW_l1=i
      write(6,*)'Number of Ow at interface','    Number of Hw'
      write(6,*)i_OW_l1,i_H
      close(10)
      !=======
      !layer 2      
      !=======
      open(10,file=trim(pos_filename))     
      do imovie=1,nmovie
          read(10,*)        
          read(10,*)       
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O' .and.&
                  x(iatom,imovie) .lt. l_2) then
                  i=i+1
              endif
          enddo
      enddo
      i_OW_l2=i
      write(6,*)'Number of Ow at interface','    Number of Hw'
      write(6,*)i_OW_l2,i_H
      close(10)
      !=======
      !layer 3      
      !=======
      open(10,file=trim(pos_filename))     
      do imovie=1,nmovie
          read(10,*)        
          read(10,*)       
          i=0
          do iatom=1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O' .and.&
                  x(iatom,imovie) .lt. l_3) then
                  i=i+1
              endif
          enddo
      enddo
      i_OW_l3=i
      write(6,*)'Number of Ow at interface','    Number of Hw'
      write(6,*)i_OW_l3,i_H
      close(10)
      !=======
      !layer 4      
      !=======
      open(10,file=trim(pos_filename))     
      do imovie=1,nmovie
          read(10,*)        
          read(10,*)       
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O' .and.&
                  x(iatom,imovie) .lt. l_4) then
                  i=i+1
              endif
          enddo
      enddo
      i_OW_l4=i
      write(6,*)'Number of Ow at interface','    Number of Hw'
      write(6,*)i_OW_l4,i_H
      close(10)
      !=======
      !layer 5      
      !=======
      open(10,file=trim(pos_filename))     
      do imovie=1,nmovie
          read(10,*)        
          read(10,*)       
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O' .and.&
                  x(iatom,imovie) .lt. l_5) then
                  i=i+1
              endif
          enddo
      enddo
      i_OW_l5=i
      write(6,*)'Number of Ow at first 5 layers','    Number of Hw'
      write(6,*)i_OW_l5,i_H
      close(10)
      !=======
      !layer 6      
      !=======
      open(10,file=trim(pos_filename))     
      do imovie=1,nmovie
          read(10,*)        
          read(10,*)       
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O' .and.&
                  x(iatom,imovie) .lt. l_6) then
                  i=i+1
              endif
          enddo
      enddo
      i_OW_l6=i
      write(6,*)'Number of Ow at first 6 layers','    Number of Hw'
      write(6,*)i_OW_l6,i_H
      close(10)
      !=======
      !layer 7      
      !=======
      open(10,file=trim(pos_filename))     
      do imovie=1,nmovie
          read(10,*)        
          read(10,*)       
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O' .and.&
                  x(iatom,imovie) .lt. l_7) then
                  i=i+1
              endif
          enddo
      enddo
      i_OW_l7=i
      write(6,*)'Number of Ow at first 7 layers','    Number of Hw'
      write(6,*)i_OW_l7,i_H
      close(10)
      !=======
      !layer 8      
      !=======
      open(10,file=trim(pos_filename))     
      do imovie=1,nmovie
          read(10,*)        
          read(10,*)       
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O' .and.&
                  x(iatom,imovie) .lt. l_8) then
                  i=i+1
              endif
          enddo
      enddo
      i_OW_l8=i
      write(6,*)'Number of Ow at first 8 layers','    Number of Hw'
      write(6,*)i_OW_l8,i_H
      close(10)
!=======================================================
! calculate the number of O (H) atoms in water molecules
!=======================================================
      allocate(ndx_OW(i_OW))! this should be put after i_OW is defined
      allocate(ndx_H(i_H))! this should be put after i_H is defined
      allocate(ndx_Li(i_Li))! this should be put after i_H is defined
      allocate(ndx_OW_l1(i_OW_l1))
      allocate(ndx_OW_l2(i_OW_l2))
      allocate(ndx_OW_l3(i_OW_l3))
      allocate(ndx_OW_l4(i_OW_l4))
      allocate(ndx_OW_l5(i_OW_l5))
      allocate(ndx_OW_l6(i_OW_l6))
      allocate(ndx_OW_l7(i_OW_l7))
      allocate(ndx_OW_l8(i_OW_l8))
      i=0
      ii=0
      iii=0
      do iatom=1,natoms
          if (trim(atom_type(iatom)) .eq. 'O')then
                 i=i+1      
                 ndx_OW(i)=iatom
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ii=ii+1
                 ndx_H(ii)=iatom
          elseif(trim(atom_type(iatom)) .eq. 'Li') then
                 iii=iii+1
                 ndx_Li(iii)=iatom
          else
          endif
      enddo 
!=======================================================
! calculate the number of O (H) atoms in water molecules
! at layers of interface     
!=======================================================
      i=0
      do imovie=1,nmovie
          do iatom=1,natoms
              if (trim(atom_type(iatom)) .eq. 'O' &
                 .and. x(iatom,imovie) .lt. l_1)then
                     i=i+1      
                     ndx_OW_l1(i)=iatom
              endif
          enddo 
      enddo
      i=0
      do imovie=1,nmovie
          do iatom=1,natoms
              if (trim(atom_type(iatom)) .eq. 'O' &
                 .and. x(iatom,imovie) .lt. l_2)then
                     i=i+1      
                     ndx_OW_l2(i)=iatom
              endif
          enddo 
      enddo
      i=0
      do imovie=1,nmovie
          do iatom=1,natoms
              if (trim(atom_type(iatom)) .eq. 'O' &
                 .and. x(iatom,imovie) .lt. l_3)then
                     i=i+1      
                     ndx_OW_l3(i)=iatom
              endif
          enddo 
      enddo
      i=0
      do imovie=1,nmovie
          do iatom=1,natoms
              if (trim(atom_type(iatom)) .eq. 'O' &
                 .and. x(iatom,imovie) .lt. l_4)then
                     i=i+1      
                     ndx_OW_l4(i)=iatom
              endif
          enddo 
      enddo
      i=0
      do imovie=1,nmovie
          do iatom=1,natoms
              if (trim(atom_type(iatom)) .eq. 'O' &
                 .and. x(iatom,imovie) .lt. l_5)then
                     i=i+1      
                     ndx_OW_l5(i)=iatom
              endif
          enddo 
      enddo
      i=0
      do imovie=1,nmovie
          do iatom=1,natoms
              if (trim(atom_type(iatom)) .eq. 'O' &
                 .and. x(iatom,imovie) .lt. l_6)then
                     i=i+1      
                     ndx_OW_l6(i)=iatom
              endif
          enddo 
      enddo
      i=0
      do imovie=1,nmovie
          do iatom=1,natoms
              if (trim(atom_type(iatom)) .eq. 'O' &
                 .and. x(iatom,imovie) .lt. l_7)then
                     i=i+1      
                     ndx_OW_l7(i)=iatom
              endif
          enddo 
      enddo
      i=0
      do imovie=1,nmovie
          do iatom=1,natoms
              if (trim(atom_type(iatom)) .eq. 'O' &
                 .and. x(iatom,imovie) .lt. l_8)then
                     i=i+1      
                     ndx_OW_l8(i)=iatom
              endif
          enddo 
      enddo
      deallocate(atom_type)
!==============
!Print the list      
!==============
      open(20,file=trim(filename)//'_list.dat')
      do i1=1, i_OW-1! No O atom can not be bonded to itself 
          m1=ndx_OW(i1)
          do i2=i1+1, i_OW 
              m2=ndx_OW(i2)
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j =1,1
                      r23= (x(m2,j)-x(m3,j))**2+  &
                           (y(m2,j)-y(m3,j))**2+  &
                           (z(m2,j)-z(m3,j))**2
                      r13= (x(m1,j)-x(m3,j))**2+  &
                           (y(m1,j)-y(m3,j))**2+  &
                           (z(m1,j)-z(m3,j))**2
                  if (r23<r_ohc) then
                      write(20,fmt='(3I4)')m1,m2,m3
                  endif
                  if (r13<r_ohc) then
                      write(20,fmt='(3I4)')m2,m1,m3
                  endif
                  enddo
              enddo
          enddo
      enddo
      close(20)
!====================================
!Print the list (layers of interface)     
!====================================
      open(20,file=trim(filename)//'_l1_list.dat')
      do i1=1, i_OW_l1-1! No O atom can not be bonded to itself 
          m1=ndx_OW_l1(i1)
          do i2=i1+1, i_OW_l1 
              m2=ndx_OW_l1(i2)
              !write(20,fmt='(2I4)',advance='no')m1,m2
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j =1,1
                      r23= (x(m2,j)-x(m3,j))**2+  &
                           (y(m2,j)-y(m3,j))**2+  &
                           (z(m2,j)-z(m3,j))**2
                      r13= (x(m1,j)-x(m3,j))**2+  &
                           (y(m1,j)-y(m3,j))**2+  &
                           (z(m1,j)-z(m3,j))**2
                  if (r23<r_ohc) then
                      write(20,fmt='(3I4)')m1,m2,m3
                  endif
                  if (r13<r_ohc) then
                      write(20,fmt='(3I4)')m2,m1,m3
                  endif
                  enddo
              enddo
          enddo
      enddo
      close(20)
      open(20,file=trim(filename)//'_l2_list.dat')
      do i1=1, i_OW_l2-1! No O atom can not be bonded to itself 
          m1=ndx_OW_l2(i1)
          do i2=i1+1, i_OW_l2 
              m2=ndx_OW_l2(i2)
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j =1,1
                      r23= (x(m2,j)-x(m3,j))**2+  &
                           (y(m2,j)-y(m3,j))**2+  &
                           (z(m2,j)-z(m3,j))**2
                      r13= (x(m1,j)-x(m3,j))**2+  &
                           (y(m1,j)-y(m3,j))**2+  &
                           (z(m1,j)-z(m3,j))**2
                  if (r23<r_ohc) then
                      write(20,fmt='(3I4)')m1,m2,m3
                  endif
                  if (r13<r_ohc) then
                      write(20,fmt='(3I4)')m2,m1,m3
                  endif
                  enddo
              enddo
          enddo
      enddo
      close(20)
      open(20,file=trim(filename)//'_l3_list.dat')
      do i1=1, i_OW_l3-1! No O atom can not be bonded to itself 
          m1=ndx_OW_l3(i1)
          do i2=i1+1, i_OW_l3
              m2=ndx_OW_l3(i2)
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j =1,1
                      r23= (x(m2,j)-x(m3,j))**2+  &
                           (y(m2,j)-y(m3,j))**2+  &
                           (z(m2,j)-z(m3,j))**2
                      r13= (x(m1,j)-x(m3,j))**2+  &
                           (y(m1,j)-y(m3,j))**2+  &
                           (z(m1,j)-z(m3,j))**2
                  if (r23<r_ohc) then
                      write(20,fmt='(3I4)')m1,m2,m3
                  endif
                  if (r13<r_ohc) then
                      write(20,fmt='(3I4)')m2,m1,m3
                  endif
                  enddo
              enddo
          enddo
      enddo
      close(20)
      open(20,file=trim(filename)//'_l4_list.dat')
      do i1=1, i_OW_l4-1! No O atom can not be bonded to itself 
          m1=ndx_OW_l4(i1)
          do i2=i1+1, i_OW_l4
              m2=ndx_OW_l4(i2)
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j =1,1
                      r23= (x(m2,j)-x(m3,j))**2+  &
                           (y(m2,j)-y(m3,j))**2+  &
                           (z(m2,j)-z(m3,j))**2
                      r13= (x(m1,j)-x(m3,j))**2+  &
                           (y(m1,j)-y(m3,j))**2+  &
                           (z(m1,j)-z(m3,j))**2
                  if (r23<r_ohc) then
                      write(20,fmt='(3I4)')m1,m2,m3
                  endif
                  if (r13<r_ohc) then
                      write(20,fmt='(3I4)')m2,m1,m3
                  endif
                  enddo
              enddo
          enddo
      enddo
      close(20)
      open(20,file=trim(filename)//'_l5_list.dat')
      do i1=1, i_OW_l5-1! No O atom can not be bonded to itself 
          m1=ndx_OW_l5(i1)
          do i2=i1+1, i_OW_l5
              m2=ndx_OW_l5(i2)
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j =1,1
                      r23= (x(m2,j)-x(m3,j))**2+  &
                           (y(m2,j)-y(m3,j))**2+  &
                           (z(m2,j)-z(m3,j))**2
                      r13= (x(m1,j)-x(m3,j))**2+  &
                           (y(m1,j)-y(m3,j))**2+  &
                           (z(m1,j)-z(m3,j))**2
                  if (r23<r_ohc) then
                      write(20,fmt='(3I4)')m1,m2,m3
                  endif
                  if (r13<r_ohc) then
                      write(20,fmt='(3I4)')m2,m1,m3
                  endif
                  enddo
              enddo
          enddo
      enddo
      close(20)
      open(20,file=trim(filename)//'_l6_list.dat')
      do i1=1, i_OW_l6-1! No O atom can not be bonded to itself 
          m1=ndx_OW_l6(i1)
          do i2=i1+1, i_OW_l6
              m2=ndx_OW_l6(i2)
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j =1,1
                      r23= (x(m2,j)-x(m3,j))**2+  &
                           (y(m2,j)-y(m3,j))**2+  &
                           (z(m2,j)-z(m3,j))**2
                      r13= (x(m1,j)-x(m3,j))**2+  &
                           (y(m1,j)-y(m3,j))**2+  &
                           (z(m1,j)-z(m3,j))**2
                  if (r23<r_ohc) then
                      write(20,fmt='(3I4)')m1,m2,m3
                  endif
                  if (r13<r_ohc) then
                      write(20,fmt='(3I4)')m2,m1,m3
                  endif
                  enddo
              enddo
          enddo
      enddo
      close(20)
      open(20,file=trim(filename)//'_l7_list.dat')
      do i1=1, i_OW_l7-1! No O atom can not be bonded to itself 
          m1=ndx_OW_l7(i1)
          do i2=i1+1, i_OW_l7
              m2=ndx_OW_l7(i2)
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j =1,1
                      r23= (x(m2,j)-x(m3,j))**2+  &
                           (y(m2,j)-y(m3,j))**2+  &
                           (z(m2,j)-z(m3,j))**2
                      r13= (x(m1,j)-x(m3,j))**2+  &
                           (y(m1,j)-y(m3,j))**2+  &
                           (z(m1,j)-z(m3,j))**2
                  if (r23<r_ohc) then
                      write(20,fmt='(3I4)')m1,m2,m3
                  endif
                  if (r13<r_ohc) then
                      write(20,fmt='(3I4)')m2,m1,m3
                  endif
                  enddo
              enddo
          enddo
      enddo
      close(20)
      open(20,file=trim(filename)//'_l8_list.dat')
      do i1=1, i_OW_l8-1! No O atom can not be bonded to itself 
          m1=ndx_OW_l8(i1)
          do i2=i1+1, i_OW_l8
              m2=ndx_OW_l8(i2)
              do ii=1,i_H
                  m3=ndx_H(ii)
                  do j =1,1
                      r23= (x(m2,j)-x(m3,j))**2+  &
                           (y(m2,j)-y(m3,j))**2+  &
                           (z(m2,j)-z(m3,j))**2
                      r13= (x(m1,j)-x(m3,j))**2+  &
                           (y(m1,j)-y(m3,j))**2+  &
                           (z(m1,j)-z(m3,j))**2
                  if (r23<r_ohc) then
                      write(20,fmt='(3I4)')m1,m2,m3
                  endif
                  if (r13<r_ohc) then
                      write(20,fmt='(3I4)')m2,m1,m3
                  endif
                  enddo
              enddo
          enddo
      enddo
      close(20)
      deallocate(ndx_OW,ndx_OW_l1,ndx_OW_l2,ndx_OW_l3,&
                 ndx_OW_l4,ndx_OW_l5,ndx_OW_l6,ndx_OW_l7,&
                 ndx_OW_l8,ndx_H,x,y,z)
!==============================================================
      call system_clock(end_time,rat) 
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
