! 2015.08.26
! autocorrelation of velocities 
! new defination of autocorrelation of velocities is given in this function to make sure the autocorrelation of velocities decay from 1 to almost zero.
! Since the later part of the autocorrelation of velocities is not accurate, only the first half part of the autocorrelation of velocities is use to calculate the VDOS. This make the VDOS much neat and clear. 
! In this function,  we consider Hydrogens and Oxygen in a water molecule separately, because we are about to find that the VDOS of two Hydrogens is proportional to the VDOS of the Oxygen!
! No  mass-weighted VDOS here, since the absolute values do not have pysical meannings.
! 含 if () then 语句
! For the autocorrelation function: e-folding decay time (Te) , degree of freeedom (N*) ,and rate of degree of freedom(N*/N) are calculated.
!*****************************************************************************************************************************    
      program velocities
      implicit none 

      integer, parameter  :: rk=8
      real(kind=rk)       :: velx1,vely1,velz1,velx2,vely2,velz2,&
                             scalar,scalar_norm,te
      character(LEN=30)  :: vel_finame
      character(LEN=3)   :: atom_type(nmax)
      integer,parameter   :: nmax=50 ! max number of atoms
      integer,parameter   :: nmovie_max = 300000 ! max number of movie steps
      integer             :: i, ii,imovie,iatom,t1,nts,nhmaxt,natoms 
      real(kind=rk), parameter :: r_e=0.36788 ! reverse of e=2.718281828.
      real(kind=rk), parameter :: Delta_t=0.5 ! the time interval (unit: fs)  
      real(kind=rk), allocatable, dimension(:,:) :: vx,vy,vz
      real(kind=rk), allocatable, dimension(:,:) :: x,y,z
      real(kind=rk), allocatable, dimension(:,:) :: corr, norm 
      real(kind=rk), allocatable, dimension(:,:) :: vcmx,vcmy,vcmz
      real(kind=rk), allocatable, dimension(:,:) :: amass 
      !dimension :: x(nmax,nmovie_max),y(nmax,nmovie_max), z(nmax,nmovie_max)
      !dimension :: corr(nmovie_max), norm(nmovie_max)
      !dimension :: amass(nmax)
      !dimension :: vcmx(nmovie_max),vcmy(nmovie_max),vcmz(nmovie_max)


! data in input 
      write(6,*) 'What is the name of velocity file:' 
      read(5,*) vel_filename
      write(6,*) 'What is the total steps of trajectory:' 
      read(5,*)nmovie          ! number of movie steps
      write(6,*) 'What is the total number of atoms:' 
      read(5,*)natoms          ! number of atoms per  molecule

      if(nmovie.ge.nmovie_max)then
         write(6,*)'!!! nmovie > nmovie_max !!!'
         write(6,*)'!!! stop !!!'
         stop
      endif
      if(natoms.ge.nmax)then
         write(6,*)'!!! natoms > nmax !!!'
         write(6,*)'!!! stop !!!'
         stop
      endif
      write(6,*)'number of movie steps:',nmovie
      write(6,*)'number of atoms in molecule:',natoms
      write(6,*)


! read in TRAJECTORY/VELOCITY file from CP2K 
      open(10,file=trim(vel_filename))

      do imovie = 1,nmovie
         write(6,*)'imovie',imovie
         read(10,*)
         read(10,*)
      allocate(vx(natoms,nmovie))
      allocate(vy(natoms,nmovie))
      allocate(vz(natoms,nmovie))
      allocate(atom_type(natoms))
         do iatom = 1,natoms
            read(10,*)atom_type(iatom), &
                      vx(iatom,imovie),vy(iatom,imovie),vz(iatom,imovie)
         enddo

      enddo
      close(10)
      write(6,*)'end of TRAJECTORY reading'

      allocate(amass(natoms))
! remove ensemble translation from velocities of atoms 
      do iatom = 1,natoms
         if(atom_type(iatom).eq.'N')amass(iatom) = 14.d0        
         if(atom_type(iatom).eq.'C')amass(iatom) = 12.d0        
         if(atom_type(iatom).eq.'H')amass(iatom) = 1.d0        
         if(atom_type(iatom).eq.'O')amass(iatom) = 16.d0        
         if(atom_type(iatom).eq.'K')amass(iatom) = 39.d0        
         write(6,*)'iatom =',iatom,'atom_type ',atom_type(iatom),&
                   'masse=',amass(iatom)
      enddo
      deallocate(atom_type)

      amass_tot = 0.d0
      do ii = 1,natoms
         amass_tot = amass_tot + amass(ii)
      enddo
      write(6,*)'masse tot', amass_tot

      do imovie = 1,nmovie
         vcmx(imovie) = 0.d0
         vcmy(imovie) = 0.d0
         vcmz(imovie) = 0.d0
         do iatom = 1,natoms
            vcmx(imovie) = vcmx(imovie) + amass(iatom)*vx(iatom,imovie)
            vcmy(imovie) = vcmy(imovie) + amass(iatom)*vy(iatom,imovie)
            vcmz(imovie) = vcmz(imovie) + amass(iatom)*vz(iatom,imovie)
         enddo
         vcmx(imovie) = vcmx(imovie)/amass_tot
         vcmy(imovie) = vcmy(imovie)/amass_tot
         vcmz(imovie) = vcmz(imovie)/amass_tot
      enddo

      do imovie = 1,nmovie
         do iatom = 1,natoms
            vx(iatom,imovie) = vx(iatom,imovie) - vcmx(imovie)
            vy(iatom,imovie) = vy(iatom,imovie) - vcmy(imovie)
            vz(iatom,imovie) = vz(iatom,imovie) - vcmz(imovie)
         enddo
      enddo
 23   continue

! _______________________________
! Total VDOS (all atoms taken into account in the sum)      
      
      do i=1,nmovie_max
         corr(i) = 0.d0
      enddo

      maxt = nmovie
      hmaxt = nmovie/2
      nhmaxt = int(hmaxt)      ! To get accurate VDOS, we only consider the first half of the correlation function.

      do mt = 0,maxt-1
         do nts = 1,maxt-mt-1
            scalar = 0.d0
            scalar_norm = 0.d0

            do iatom=1,natoms
               scalar = scalar +                   & 
                (vx(iatom,nts)*vx(iatom,nts+mt) +  &
                vy(iatom,nts)*vy(iatom,nts+mt) +   &
                vz(iatom,nts)*vz(iatom,nts+mt))    
               scalar_norm = scalar_norm+(vx(iatom,nts)*vx(iatom,nts)+ &
                vy(iatom,nts)*vy(iatom,nts) +      &
                vz(iatom,nts)*vz(iatom,nts)) ! scalar product of velocities
            enddo
            corr(mt+1) = corr(mt+1) + scalar
         enddo
         corr(mt+1) = corr(mt+1)/dble(maxt-mt)
      enddo

      open(10,file='correlation_allvels.dat')
      do i = 1,nhmaxt
         write(10,*)i-1,corr(i)
      enddo
      write(6,*)'correlation results written in file correlation'
      close(10)
! _______________________________
! Order of atoms :  N(1) O(2) O(3) O(4) (NO3-), O(5) H(6) H(7) (Wat1), O(8) H(9) H(10) (Wat2), O(11) H(12) H(13) (Wat3)

      do i=1,nmovie_max
         corr(i) = 0.d0
      enddo

      maxt = nmovie
      do mt = 0,maxt-1
         do nts = 1,maxt-mt-1
            scalar = 0.d0

            do iatom=1,1 ! N of NO3-
               scalar = scalar +                  & 
                vx(iatom,nts)*vx(iatom,nts+mt) +  &
                vy(iatom,nts)*vy(iatom,nts+mt) +  &
                vz(iatom,nts)*vz(iatom,nts+mt)  
            enddo
            corr(mt+1) = corr(mt+1) + scalar
         enddo
         corr(mt+1) = corr(mt+1)/dble(maxt-mt)
      enddo

      do mt =0, maxt-1     
      if (corr(mt+1) < r_e) then
          t1=mt  ! t1 is a real type varible 
          Te=t1*Delta_t
      exit 
      endif
      enddo 
 
      open (10, file= 'degree_freedom_N.dat')
        write(10,*)  'e-folding decay time(Te):', &
                    te,'(fs)'   !"," should not be negnected. 
        write(10,*)  'degree of freeedom(N*):',maxt/(2*t1)
        write(10,*)  'rate of degree of freedom(N*/N):',1/(2*t1)  
      close(10)

      open(10,file= 'correlation_N.dat')
      do i = 1,nhmaxt
         write(10,*)i-1,corr(i)
      enddo
      close(10)
! _______________________________
      END
