!2015/10/17
!========================================
!For ONE pair of molecules，we can obtain
!C_HB(t):
!   C_{HB}(t)= <h(0)h(t)>/<h>
!===================================================
!We can implementate these correlation functions by
!calculate the average over N pairs of molecules. 
!===================================================      
! For any time step, instead of the original time 
! step.      
!===================================================     
! THIS FUNCTION CONSIDERED THE PBC, THEREFORE, THE 
! PARAMETERS a,b and c ARE REQUIRED.
!===================================================     
! input file: input
     ! name of system
     ! name of trajectory
     ! name of list
     ! nmo
     ! nat
     ! number of pairs of molecules(np) 
     ! the new time step, or ns
     ! size: la,lb,lc
!=======================================
      program hbacf_C_function
      implicit none
!==========
!parameters
!==========
      character(LEN=200) :: filename ,pos_filename,list_filename       
      integer,parameter :: rk=4              
      real(kind=rk),parameter :: rate=0.80        ! condition for cutting off autocorrelation functions
      real,parameter :: rooc=12.25                ! cutoff distance of rOO (3.5**2 )
      real,parameter :: rohc=6.0025                ! rOH (2.45**2)
      real,parameter :: cosphic=0.866             ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: h_min=0.5 ! condition for the existence of a h-bond for a step
      real(kind=rk),parameter :: hb_min=0.5 ! condition for the existence of h-bond for a pair of water molecules
      real(kind=rk)           :: r12,r13,r23,cosphi,pm,qj,&
                                 tot_hb,delta_t,&
                                 r12_1,r13_1,r23_1,cosphi_1,pm_1,&
                                 r12_2,r13_2,r23_2,cosphi_2,pm_2,&
                                 r12_3,r13_3,r23_3,cosphi_3,pm_3,&
                                 r12_4,r13_4,r23_4,cosphi_4,pm_4,&
                                 r12_5,r13_5,r23_5,cosphi_5,pm_5,&
                                 r12_6,r13_6,r23_6,cosphi_6,pm_6,&
                                 r12_11,r13_11,r23_11,cosphi_11,pm_11,&
                                 r12_22,r13_22,r23_22,cosphi_22,pm_22,&
                                 r12_33,r13_33,r23_33,cosphi_33,pm_33,&
                                 r12_44,r13_44,r23_44,cosphi_44,pm_44,&
                                 r12_55,r13_55,r23_55,cosphi_55,pm_55,&
                                 r12_66,r13_66,r23_66,cosphi_66,pm_66,&
                                 la,lb,lc
      integer :: begin_time,end_time,rat,&
                 i,j,k,jj,nmo,nat,iatom,& 
                 imovie,np,m1,m2,m3,mt,&
                 nqj,tot_nhb,n_bonded_pairs,ns 
      real(kind=rk),allocatable,dimension (:)    :: h,hb,corr_h
      real,allocatable,dimension (:,:)           :: x,y,z,x1,y1,z1,&
                                                    x11,y11,z11
      character(LEN=3)  :: atom_type  ! we are not interested in the atom type in this calculation,
                                      ! thus I do not use allocatable array for it.
      integer,allocatable,dimension(:)           :: ndx_1,ndx_2,ndx_3,&
          nhb_exist
      real(kind=rk)  :: scalar 
      logical,allocatable,dimension (:)  :: hb_exist
      call system_clock(begin_time,rat)
!==================
!read data in input
!==================
      write(6,*)'What is the name of the system:'
      read(5,*)filename
      write(6,*)'What is the name of the trajecotry file:'
      read(5,*)pos_filename     
      write(6,*)'What is the name of thei list file:'
      read(5,*)list_filename     
      write(6,*)'What is the total steps of the trajecotry:'
      read(5,*)nmo     !number of movie steps
      write(6,*)'What is the total number of atoms in the system:'
      read(5,*)nat     !number of atoms per mole.
      write(6,*)'What is the total number of water pairs:'
      read(5,*)np      !number of pairs   
      write(6,*)'What is the time step for calculating CORRELATION:' 
      read(5,*)ns! [ns*0.0005] ps is the new time step for calculating correl func.
      write(6,*)'SIZE OF EACH DIMENSION (a,b,c) (unit: Angstrom):' 
      read(5,*)la,lb,lc
      
      list_filename=trim(list_filename)
      allocate(ndx_1(np))          
      allocate(ndx_2(np))          
      allocate(ndx_3(np))          
      open(10,file=list_filename)     
      do k=1,np
          read(10,*)ndx_1(k),ndx_2(k),ndx_3(k)
      enddo
      close(10)

      delta_t=ns*0.0005 ! unit: ps
      nmo=nmo/ns ! Length of the correl. function 
      allocate(x(nat,nmo))
      allocate(y(nat,nmo))
      allocate(z(nat,nmo)) 
      allocate(x1(nat,nmo)) ! FOR PBC along x axis
      allocate(y1(nat,nmo)) ! FOR PBC along x axis
      allocate(z1(nat,nmo)) ! FOR PBC along x axis
      allocate(x11(nat,nmo)) ! FOR PBC along y axis
      allocate(y11(nat,nmo)) ! FOR PBC along y axis
      allocate(z11(nat,nmo)) ! FOR PBC along y axis
      allocate(h(nmo))
      allocate(hb(np)) ! Average H-bonded population 
      allocate(nhb_exist(np))
!=======================
!read in trajectory file 
!=======================
      open(10,file=trim(pos_filename))     
      do imovie=1,nmo
         read(10,*)!Neglect data of this line
         read(10,*)                 
         do iatom= 1,nat
             read (10,*)atom_type,x(iatom,imovie),y(iatom,imovie),&
                        z(iatom,imovie)
             x1(iatom,imovie)=x(iatom,imovie)+la
             y1(iatom,imovie)=y(iatom,imovie)
             z1(iatom,imovie)=z(iatom,imovie)
             x11(iatom,imovie)=x(iatom,imovie)
             y11(iatom,imovie)=y(iatom,imovie)+lb
             z11(iatom,imovie)=z(iatom,imovie)
         enddo

         do i=1, (nat+2)*(ns-1)
             read(10,*)
         enddo

      enddo
      close(10)
      write(6,*) 'end of trajectory reading'
!====================================
!Calculate autocorrelation function
!====================================
! calculate <h(0)h(t)>/<h>  
! Notice here <> is not average over
! different pairs of water molecules,
! and over all starting time points i
! with h(i)=1.
!====================================      
      allocate(corr_h(nmo))
      allocate(hb_exist(nmo))
      !do i=1, nmo
      corr_h(:)=0.0
      !enddo
      tot_hb=0.0
      tot_nhb=0
      !do k=1,np
      hb(:)=0.0
      nhb_exist(:)=0
      !enddo
!=============
!The main loop
!=============      
      do k=1,np
        qj=0
        nqj=0
        m1=ndx_1(k)
        m2=ndx_2(k)
        m3=ndx_3(k)
        !Calculate h(j)
        do jj =1, nmo
          h(jj)=0.0
          hb_exist(jj)=.False.
          r13= (x(m1,jj)-x(m3,jj))**2+       &
                    (y(m1,jj)-y(m3,jj))**2+  &
                    (z(m1,jj)-z(m3,jj))**2!r12,r13,r23:square of distances
          r12= (x(m1,jj)-x(m2,jj))**2+       &
                    (y(m1,jj)-y(m2,jj))**2+  &
                    (z(m1,jj)-z(m2,jj))**2
          r23= (x(m2,jj)-x(m3,jj))**2+       &
                    (y(m2,jj)-y(m3,jj))**2+  &
                    (z(m2,jj)-z(m3,jj))**2
          pm= (x(m3,jj)-x(m2,jj))*           &
                   (x(m1,jj)-x(m2,jj))+      & 
                   (y(m3,jj)-y(m2,jj))*      & 
                   (y(m1,jj)-y(m2,jj))+      & 
                   (z(m3,jj)-z(m2,jj))*      &
                   (z(m1,jj)-z(m2,jj))!pm: point multiplication. 
          cosphi= pm/(sqrt(r23*r12))
          !========================================
          ! x-aixs-mirror-original cell interaction
          !========================================
          r13_1= (x1(m1,jj)-x1(m3,jj))**2+       &
                    (y1(m1,jj)-y1(m3,jj))**2+  &
                    (z1(m1,jj)-z1(m3,jj))**2!r12,r13,r23:square of distances
          r12_1= (x1(m1,jj)-x(m2,jj))**2+       &
                    (y1(m1,jj)-y(m2,jj))**2+  &
                    (z1(m1,jj)-z(m2,jj))**2
          r23_1= (x(m2,jj)-x1(m3,jj))**2+       &
                    (y(m2,jj)-y1(m3,jj))**2+  &
                    (z(m2,jj)-z1(m3,jj))**2
          pm_1= (x1(m3,jj)-x(m2,jj))*           &
                   (x1(m1,jj)-x(m2,jj))+      & 
                   (y1(m3,jj)-y(m2,jj))*      & 
                   (y1(m1,jj)-y(m2,jj))+      & 
                   (z1(m3,jj)-z(m2,jj))*      &
                   (z1(m1,jj)-z(m2,jj))!pm: point multiplication. 
          cosphi_1= pm_1/(sqrt(r23_1*r12_1))
          r13_2= (x1(m1,jj)-x(m3,jj))**2+       &
                    (y1(m1,jj)-y(m3,jj))**2+  &
                    (z1(m1,jj)-z(m3,jj))**2!r12,r13,r23:square of distances
          r12_2= (x1(m1,jj)-x(m2,jj))**2+       &
                    (y1(m1,jj)-y(m2,jj))**2+  &
                    (z1(m1,jj)-z(m2,jj))**2
          r23_2= (x(m2,jj)-x(m3,jj))**2+       &
                    (y(m2,jj)-y(m3,jj))**2+  &
                    (z(m2,jj)-z(m3,jj))**2
          pm_2= (x(m3,jj)-x(m2,jj))*           &
                   (x1(m1,jj)-x(m2,jj))+      & 
                   (y(m3,jj)-y(m2,jj))*      & 
                   (y1(m1,jj)-y(m2,jj))+      & 
                   (z(m3,jj)-z(m2,jj))*      &
                   (z1(m1,jj)-z(m2,jj))!pm: point multiplication. 
          cosphi_2= pm_2/(sqrt(r23_2*r12_2))
          r13_3= (x(m1,jj)-x(m3,jj))**2+       &
                    (y(m1,jj)-y(m3,jj))**2+  &
                    (z(m1,jj)-z(m3,jj))**2!r12,r13,r23:square of distances
          r12_3= (x(m1,jj)-x1(m2,jj))**2+       &
                    (y(m1,jj)-y1(m2,jj))**2+  &
                    (z(m1,jj)-z1(m2,jj))**2
          r23_3= (x1(m2,jj)-x(m3,jj))**2+       &
                    (y1(m2,jj)-y(m3,jj))**2+  &
                    (z1(m2,jj)-z(m3,jj))**2
          pm_3= (x(m3,jj)-x1(m2,jj))*           &
                   (x(m1,jj)-x1(m2,jj))+      & 
                   (y(m3,jj)-y1(m2,jj))*      & 
                   (y(m1,jj)-y1(m2,jj))+      & 
                   (z(m3,jj)-z1(m2,jj))*      &
                   (z(m1,jj)-z1(m2,jj))!pm: point multiplication. 
          cosphi_3= pm_3/(sqrt(r23_3*r12_3))
          r13_4= (x(m1,jj)-x1(m3,jj))**2+       &
                    (y(m1,jj)-y1(m3,jj))**2+  &
                    (z(m1,jj)-z1(m3,jj))**2!r12,r13,r23:square of distances
          r12_4= (x(m1,jj)-x1(m2,jj))**2+       &
                    (y(m1,jj)-y1(m2,jj))**2+  &
                    (z(m1,jj)-z1(m2,jj))**2
          r23_4= (x1(m2,jj)-x1(m3,jj))**2+       &
                    (y1(m2,jj)-y1(m3,jj))**2+  &
                    (z1(m2,jj)-z1(m3,jj))**2
          pm_4= (x1(m3,jj)-x1(m2,jj))*           &
                   (x(m1,jj)-x1(m2,jj))+      & 
                   (y1(m3,jj)-y1(m2,jj))*      & 
                   (y(m1,jj)-y1(m2,jj))+      & 
                   (z1(m3,jj)-z1(m2,jj))*      &
                   (z(m1,jj)-z1(m2,jj))!pm: point multiplication. 
          cosphi_4= pm_4/(sqrt(r23_4*r12_4))
          r13_5= (x(m1,jj)-x1(m3,jj))**2+       &
                    (y(m1,jj)-y1(m3,jj))**2+  &
                    (z(m1,jj)-z1(m3,jj))**2!r12,r13,r23:square of distances
          r12_5= (x(m1,jj)-x(m2,jj))**2+       &
                    (y(m1,jj)-y(m2,jj))**2+  &
                    (z(m1,jj)-z(m2,jj))**2
          r23_5= (x(m2,jj)-x1(m3,jj))**2+       &
                    (y(m2,jj)-y1(m3,jj))**2+  &
                    (z(m2,jj)-z1(m3,jj))**2
          pm_5= (x1(m3,jj)-x(m2,jj))*           &
                   (x(m1,jj)-x(m2,jj))+      & 
                   (y1(m3,jj)-y(m2,jj))*      & 
                   (y(m1,jj)-y(m2,jj))+      & 
                   (z1(m3,jj)-z(m2,jj))*      &
                   (z(m1,jj)-z(m2,jj))!pm: point multiplication. 
          cosphi_5= pm_5/(sqrt(r23_5*r12_5))
          r13_6= (x1(m1,jj)-x(m3,jj))**2+       &
                    (y1(m1,jj)-y(m3,jj))**2+  &
                    (z1(m1,jj)-z(m3,jj))**2!r12,r13,r23:square of distances
          r12_6= (x1(m1,jj)-x1(m2,jj))**2+       &
                    (y1(m1,jj)-y1(m2,jj))**2+  &
                    (z1(m1,jj)-z1(m2,jj))**2
          r23_6= (x1(m2,jj)-x(m3,jj))**2+       &
                    (y1(m2,jj)-y(m3,jj))**2+  &
                    (z1(m2,jj)-z(m3,jj))**2
          pm_6= (x(m3,jj)-x1(m2,jj))*           &
                   (x1(m1,jj)-x1(m2,jj))+      & 
                   (y(m3,jj)-y1(m2,jj))*      & 
                   (y1(m1,jj)-y1(m2,jj))+      & 
                   (z(m3,jj)-z1(m2,jj))*      &
                   (z1(m1,jj)-z1(m2,jj))!pm: point multiplication. 
          cosphi_6= pm_6/(sqrt(r23_6*r12_6))
          !========================================
          ! y-aixs-mirror-original cell interaction
          !========================================
          r13_11= (x11(m1,jj)-x11(m3,jj))**2+       &
                    (y11(m1,jj)-y11(m3,jj))**2+  &
                    (z11(m1,jj)-z11(m3,jj))**2!r12,r13,r23:square of distances
          r12_11= (x11(m1,jj)-x(m2,jj))**2+       &
                    (y11(m1,jj)-y(m2,jj))**2+  &
                    (z11(m1,jj)-z(m2,jj))**2
          r23_11= (x(m2,jj)-x11(m3,jj))**2+       &
                    (y(m2,jj)-y11(m3,jj))**2+  &
                    (z(m2,jj)-z11(m3,jj))**2
          pm_11= (x11(m3,jj)-x(m2,jj))*           &
                   (x11(m1,jj)-x(m2,jj))+      & 
                   (y11(m3,jj)-y(m2,jj))*      & 
                   (y11(m1,jj)-y(m2,jj))+      & 
                   (z11(m3,jj)-z(m2,jj))*      &
                   (z11(m1,jj)-z(m2,jj))!pm: point multiplication. 
          cosphi_11= pm_11/(sqrt(r23_11*r12_11))
          r13_22= (x11(m1,jj)-x(m3,jj))**2+       &
                    (y11(m1,jj)-y(m3,jj))**2+  &
                    (z11(m1,jj)-z(m3,jj))**2!r12,r13,r23:square of distances
          r12_22= (x11(m1,jj)-x(m2,jj))**2+       &
                    (y11(m1,jj)-y(m2,jj))**2+  &
                    (z11(m1,jj)-z(m2,jj))**2
          r23_22= (x(m2,jj)-x(m3,jj))**2+       &
                    (y(m2,jj)-y(m3,jj))**2+  &
                    (z(m2,jj)-z(m3,jj))**2
          pm_22= (x(m3,jj)-x(m2,jj))*           &
                   (x11(m1,jj)-x(m2,jj))+      & 
                   (y(m3,jj)-y(m2,jj))*      & 
                   (y11(m1,jj)-y(m2,jj))+      & 
                   (z(m3,jj)-z(m2,jj))*      &
                   (z11(m1,jj)-z(m2,jj))!pm: point multiplication. 
          cosphi_22= pm_22/(sqrt(r23_22*r12_22))
          r13_33= (x(m1,jj)-x(m3,jj))**2+       &
                    (y(m1,jj)-y(m3,jj))**2+  &
                    (z(m1,jj)-z(m3,jj))**2!r12,r13,r23:square of distances
          r12_33= (x(m1,jj)-x11(m2,jj))**2+       &
                    (y(m1,jj)-y11(m2,jj))**2+  &
                    (z(m1,jj)-z11(m2,jj))**2
          r23_33= (x11(m2,jj)-x(m3,jj))**2+       &
                    (y11(m2,jj)-y(m3,jj))**2+  &
                    (z11(m2,jj)-z(m3,jj))**2
          pm_33= (x(m3,jj)-x11(m2,jj))*           &
                   (x(m1,jj)-x11(m2,jj))+      & 
                   (y(m3,jj)-y11(m2,jj))*      & 
                   (y(m1,jj)-y11(m2,jj))+      & 
                   (z(m3,jj)-z11(m2,jj))*      &
                   (z(m1,jj)-z11(m2,jj))!pm: point multiplication. 
          cosphi_33= pm_33/(sqrt(r23_33*r12_33))
          r13_44= (x(m1,jj)-x11(m3,jj))**2+       &
                    (y(m1,jj)-y11(m3,jj))**2+  &
                    (z(m1,jj)-z11(m3,jj))**2!r12,r13,r23:square of distances
          r12_44= (x(m1,jj)-x11(m2,jj))**2+       &
                    (y(m1,jj)-y11(m2,jj))**2+  &
                    (z(m1,jj)-z11(m2,jj))**2
          r23_44= (x11(m2,jj)-x11(m3,jj))**2+       &
                    (y11(m2,jj)-y11(m3,jj))**2+  &
                    (z11(m2,jj)-z11(m3,jj))**2
          pm_44= (x11(m3,jj)-x11(m2,jj))*           &
                   (x(m1,jj)-x11(m2,jj))+      & 
                   (y11(m3,jj)-y11(m2,jj))*      & 
                   (y(m1,jj)-y11(m2,jj))+      & 
                   (z11(m3,jj)-z11(m2,jj))*      &
                   (z(m1,jj)-z11(m2,jj))!pm: point multiplication. 
          cosphi_44= pm_44/(sqrt(r23_44*r12_44))
          r13_55= (x(m1,jj)-x11(m3,jj))**2+       &
                    (y(m1,jj)-y11(m3,jj))**2+  &
                    (z(m1,jj)-z11(m3,jj))**2!r12,r13,r23:square of distances
          r12_55= (x(m1,jj)-x(m2,jj))**2+       &
                    (y(m1,jj)-y(m2,jj))**2+  &
                    (z(m1,jj)-z(m2,jj))**2
          r23_55= (x(m2,jj)-x11(m3,jj))**2+       &
                    (y(m2,jj)-y11(m3,jj))**2+  &
                    (z(m2,jj)-z11(m3,jj))**2
          pm_55= (x1(m3,jj)-x(m2,jj))*           &
                   (x(m1,jj)-x(m2,jj))+      & 
                   (y11(m3,jj)-y(m2,jj))*      & 
                   (y(m1,jj)-y(m2,jj))+      & 
                   (z11(m3,jj)-z(m2,jj))*      &
                   (z(m1,jj)-z(m2,jj))!pm: point multiplication. 
          cosphi_55= pm_55/(sqrt(r23_55*r12_55))
          r13_66= (x11(m1,jj)-x(m3,jj))**2+       &
                    (y11(m1,jj)-y(m3,jj))**2+  &
                    (z11(m1,jj)-z(m3,jj))**2!r12,r13,r23:square of distances
          r12_66= (x11(m1,jj)-x11(m2,jj))**2+       &
                    (y11(m1,jj)-y11(m2,jj))**2+  &
                    (z11(m1,jj)-z11(m2,jj))**2
          r23_66= (x11(m2,jj)-x(m3,jj))**2+       &
                    (y11(m2,jj)-y(m3,jj))**2+  &
                    (z11(m2,jj)-z(m3,jj))**2
          pm_66= (x(m3,jj)-x11(m2,jj))*           &
                   (x11(m1,jj)-x11(m2,jj))+      & 
                   (y(m3,jj)-y11(m2,jj))*      & 
                   (y11(m1,jj)-y11(m2,jj))+      & 
                   (z(m3,jj)-z11(m2,jj))*      &
                   (z11(m1,jj)-z11(m2,jj))!pm: point multiplication. 
          cosphi_66= pm_66/(sqrt(r23_66*r12_66))
          if (                                       & 
             (r13 .lt. rohc .and. r12 .lt. rooc      & 
             .and. cosphi .gt. cosphic) .or.         &
             (r13_44 .lt. rohc .and. r12_44 .lt. rooc      & 
             .and. cosphi_44 .gt. cosphic) .or.         &
             (r13_33 .lt. rohc .and. r12_33 .lt. rooc      & 
             .and. cosphi_33 .gt. cosphic) .or.         &
             (r13_22 .lt. rohc .and. r12_22 .lt. rooc      & 
             .and. cosphi_22 .gt. cosphic) .or.         &
             (r13_11 .lt. rohc .and. r12_11 .lt. rooc      & 
             .and. cosphi_11 .gt. cosphic) .or.         &
             (r13_5 .lt. rohc .and. r12_5 .lt. rooc      & 
             .and. cosphi_5 .gt. cosphic) .or.         &
             (r13_4 .lt. rohc .and. r12_4 .lt. rooc      & 
             .and. cosphi_4 .gt. cosphic) .or.         &
             (r13_3 .lt. rohc .and. r12_3 .lt. rooc      & 
             .and. cosphi_3 .gt. cosphic) .or.         &
             (r13_2 .lt. rohc .and. r12_2 .lt. rooc      & 
             .and. cosphi_2 .gt. cosphic) .or.         &
             (r13_1 .lt. rohc .and. r12_1 .lt. rooc  & 
             .and. cosphi_1 .gt. cosphic)            & 
             )then    
              h(jj)=1.0 
              hb_exist(jj)=.True.
              qj=qj+h(jj)!To calculate ave population of hbond over all starting points for one pair of water molecules.                          
              nqj=nqj+1
          endif
        enddo   
        !qj=qj/nmo!Ave of hb over all starting points for each pair 
        hb(k)=qj 
        nhb_exist(k)=nqj
        tot_hb=tot_hb+hb(k)
        tot_nhb=tot_nhb+nhb_exist(k)
        !==================================
        !Calcualte the correlation function
        !==================================
        if (hb(k)>hb_min) then
            do mt=0,nmo-1! time interval
                scalar=0.d0
                do j=1,nmo-mt-1
                    scalar=scalar+h(j)*h(j+mt)  
                    enddo
                corr_h(mt+1)=corr_h(mt+1)+scalar! sum_C_k(t)
            enddo
        endif
      enddo!End of k-loop 
      deallocate(hb_exist,nhb_exist)
      !=========================================
      !Calculate the number of ever bonded pairs
      !=========================================
      n_bonded_pairs=0 
      do k=1,np
          if (hb(k)>hb_min) then
              n_bonded_pairs=n_bonded_pairs+1      
          endif
      enddo
      !========================
      !Normalization of C_HB(t)
      !========================
      do mt=0,nmo-1! time interval
          corr_h(mt+1)=corr_h(mt+1)/tot_nhb   
      enddo
      deallocate(x,y,z,ndx_1,ndx_2,ndx_3)          
!======================
!Write the correlation
!C_HB(t)     
!======================
      open(10,file=trim(filename)//'_pbc_hbacf_h.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,corr_h(i)
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_hbacf_h.dat'
      close(10)
!===========================
!Write the correlation
!ln(C_HB(t))     
!===========================
      open(10,file=trim(filename)//'_pbc_hbacf_ln_h.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,log(corr_h(i))
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_hbacf_ln_h.dat'
      close(10)
!===========
! Print <h>      
!===========      
      open(10,file=trim(filename)//'_pbc_ave_h.dat')
        write(10,*) 'Ave. No. bonds:',tot_hb/nmo
        write(10,*) '<h>:',(tot_hb/nmo)/np
        write(6,*)'written in '//trim(filename)//&
                  '_ave_h.dat'
      close(10)
      deallocate (h,corr_h,hb)
!=======================
!Print ending time point
!=======================      
      call system_clock(end_time,rat)
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
