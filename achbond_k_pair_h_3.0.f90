!============
! 2015/09/1
! huang gang
!============

!=========================================      
!For specified k pairs of water molecules，
!to calculate the average and the auto-
!correlation function of the H-Bond 
!population operator:
! <h>;      
! C_HB(t).

     ! input file: input
     ! the first line of input: nmovie
     ! the 2nd line: natoms
     ! the 3rd line: number of pairs of
     ! molecules, or npair 
     ! the 4th-(npair-3)th: the index 
     ! of O,O,H atoms which will be H-
     ! bonded. 
!=========================================
      
      program achbond_k_pair
      implicit none

!========================
!parameters and variables
!========================
      character(LEN=20) :: filename            ! For specific filename to analyzing data
      integer,parameter :: rk=4              
      integer,parameter :: nmax=775            ! max number of atoms
      integer,parameter :: nmovie_max=300000   ! max number of movie
      real,parameter :: rate=0.4               ! cutoff rate of total correlation 
                                               ! (determine the length of the correlation data)
      real,parameter :: rooc=29.16             ! cutoff distance of rOO (5.4**2 )
      real,parameter :: rohc=12.25             ! rOH (3.5**2)
      real,parameter :: cosphic=0.866          ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: delta_t=0.0005 ! fs
      real           :: ave_h
      integer        :: i,j,k,nmovie,natoms,iatom,& 
                        imovie,npair,m1,m2,m3,mt
      real(kind=rk),allocatable,dimension (:,:) :: r12, r13, r23
      real(kind=rk),allocatable,dimension (:,:) :: h
      real(kind=rk),allocatable,dimension (:,:) :: cosphi, pm
      real,allocatable,dimension (:,:)         :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:)   :: ndx_1, ndx_2, ndx_3
      real,allocatable,dimension (:)      :: hb      ! For calcualate the average HB population
      real(kind=rk),allocatable,dimension (:) :: corr_h
      real(kind=rk),allocatable,dimension (:) :: scalar_h      ! Auxiliary variables
      
!==================
!read data in input
!==================
      open(10,file='achbond_k_pair_input')
      read(10,*)filename
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
      allocate(cosphi(nmovie,npair))
      allocate(pm(nmovie,npair))
      allocate(r12(nmovie,npair))
      allocate(r13(nmovie,npair))
      allocate(r23(nmovie,npair))
      allocate(h(nmovie,npair))
      allocate(hb(npair))            !Average H-bonded population 

!=======================
!read in trajectory file 
!=======================
      open(10,file='traj_pos.xyz')     
      do imovie=1,nmovie
         read(10,*)                  !Neglect data of this line
         read(10,*)                  !Neglect data of this line
         do iatom= 1,natoms
            read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
         enddo
      enddo
      close(10)
      write(6,*)'end of trajectory reading'
!===========================key part====================================
      write(6,*)'For test:','ndx_1(1)=',ndx_1(1),&
                ', ndx_2(1)=',ndx_2(1)

      do k=1,npair
          hb(k)=0
          m1=ndx_1(k)
          m2=ndx_2(k)
          m3=ndx_3(k)

          do j =1, nmovie
              h(j,k)=0
              r13(j,k)= (x(m1,j)-x(m3,j))**2+  &
                  (y(m1,j)-y(m3,j))**2+  &
                  (z(m1,j)-z(m3,j))**2      !r:squra of distances
              r12(j,k)= (x(m1,j)-x(m2,j))**2+  &
                  (y(m1,j)-y(m2,j))**2+  &
                  (z(m1,j)-z(m2,j))**2
              r23(j,k)= (x(m2,j)-x(m3,j))**2+  &
                  (y(m2,j)-y(m3,j))**2+  &
                  (z(m2,j)-z(m3,j))**2
              pm(j,k)= (x(m3,j)-x(m2,j))*   &
                  (x(m1,j)-x(m2,j))+   & 
                  (y(m3,j)-y(m2,j))*   & 
                  (y(m1,j)-y(m2,j))+   & 
                  (z(m3,j)-z(m2,j))*   &
                  (z(m1,j)-z(m2,j)) 
              cosphi(j,k)= pm(j,k)/(sqrt(r23(j,k)*r12(j,k)))            !pm: point multiplication.
              if (r13(j,k) .lt. rohc .and. r12(j,k).lt.rooc          & 
                  .and. cosphi(j,k).gt. cosphic) then    
                  h(j,k)=1.0 
                  hb(k)=hb(k)+h(j,k)                          
              endif
          enddo   
          hb(k)=hb(k)/nmovie
      enddo
      deallocate (x,y,z,atom_type)
!================
!Write the result
!================
      
      !For testing
      open(10,file=trim(filename)//'_nhb_test.dat')
      write(10,*)'#pair(k) ','Ave HBond population(h_k) ','hb_d_k'
                  
      do k=1,npair
         write(10,*) k, hb(k)
      enddo
      write(6,*)'written in '//trim(filename)//'_nhb_test.dat'
      close(10)

      !Write the result
      open(10,file=trim(filename)//'_nhb.dat')

      do k=1,npair
          write(10,*)'#pair(k) ','Ave HBond population(h_k) '
                
          write(10,*) k, hb(k)
          write (10,*)'#Step ','time(ps)','Population(h_k(t))  '      

          do j =1,nmovie
              write(10,*)j,j*delta_t,h(j,k)
          enddo
      enddo

      write(6,*)'written in '//trim(filename)//'_nhb.dat'
      close(10)
      deallocate (cosphi,pm,                                          & 
                 r12,r13,r23,ndx_1,ndx_2,ndx_3)

!==================================
!Calculate autocorrelation function
!==================================
      !initialization   
      allocate(corr_h(nmovie))
      allocate(scalar_h(npair))
      do i=1, nmovie
          corr_h(i)=0
      enddo
      ave_h=0.d0               ! ave_h are average hb population over all pairs

      ! calculate the ave_h,
      do k=1, npair
          ave_h=ave_h+hb(k)
      enddo
      ave_h=ave_h/npair
      write(6,*) ave_h
      
      ! calculate <h(0)h(t)>/<h> 
      ! Notice here <> is not average over
      ! different pairs of water molecules,
      ! but average over the time steps.
      do mt=0,nmovie-1     ! time interval
          do k=1, npair
              scalar_h(k)=0.d0
              do j=1, nmovie-mt-1
                  scalar_h(k)=scalar_h(k)+h(j,k)*h(j+mt,k)      ! k: the k_th pair of water molecules
              enddo
              scalar_h(k)=scalar_h(k)/(nmovie-mt)
              corr_h(mt+1)=corr_h(mt+1)+scalar_h(k)
          enddo  
          corr_h(mt+1)=corr_h(mt+1)/(npair*ave_h)
      enddo    

!=====================
!Write the correlation
!C_HB(t)      
!=====================
      open(10,file=trim(filename)//'_acf_h.dat')
          do i=1,int(nmovie*rate)
              write(10,*)i-1,corr_h(i)
          enddo
          write(6,*)'written in '//trim(filename)//'_acf_h.dat'
      close(10)

!==========================
!Write the correlation
!ln(C_HB(t))     
!==========================
      open(10,file=trim(filename)//'_ln_acf_h.dat')
          do i=1,int(nmovie*rate)
              write(10,*)i-1,log(corr_h(i))
          enddo
          write(6,*)'written in '//trim(filename)//'_ln_acf_h.dat'
      close(10)

      deallocate (h,hb,scalar_h,corr_h)
!==============================================================
      STOP
      END
