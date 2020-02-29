!2015/09/13
! modified on 2015/09/25
!=====================================
!对指定的 1 对原子，可求其氢键布居算符之
!平均值和其自关联函数S_HB(t),S_d_HB(t),
!C_HB(t),C_d_HB(t)。
!      
! C_{HB}(t)= <h(0)h(t)>/<h>
! C^(d)_{HB}(t)= <h(0)h^(d)(t)>/<h>
! 用此方法可得 N 对水分子之氢键布居算符. 
     ! input file: input
     ! name of system
     ! name of trajectory
     ! name of list
     ! nmo
     ! nat
     ! number of pairs of molecules(np) 
!=====================================
      program hbacf2
      implicit none
!==========
!parameters
!==========
      character(LEN=30) :: filename ,pos_filename,list_filename         
      integer,parameter :: rk=4              
      integer,parameter :: nmax=2001               ! max number of atoms
      integer,parameter :: nmo_max=300000         ! max number of movie
      real(kind=rk),parameter :: rate=0.50        ! condition for cutting off autocorrelation functions
      real,parameter :: rooc=29.16                ! cutoff distance of rOO (5.4**2 )
      real,parameter :: rohc=12.25                ! rOH (3.5**2)
      real,parameter :: cosphic=0.866             ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: hb_min=0.0000001 ! condition for the existence of h-bond
      real(kind=rk),parameter :: delta_t=0.01   ! ps
      real(kind=rk)           :: r12,r13,r23,cosphi,pm,qj,qj_d,&
                                 tot_hb
      integer :: begin_time,end_time,rat,&
                        i,j,k,jj,nmo,nat,iatom,& 
                        imovie,np,m1,m2,m3,mt 
      real(kind=rk),allocatable,dimension (:)    :: h,h_d
      real(kind=rk),allocatable,dimension (:)    :: hb,hb_d
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3)  :: atom_type  ! we are not interested in the atom type in this calculation,
                                      ! thus I do not use allocatable array for it.
      integer,allocatable,dimension(:)           :: ndx_1, ndx_2, ndx_3
      real(kind=rk),allocatable,dimension (:)    :: corr_h, corr_h_d
      real(kind=rk)  :: scalar, scalar_d 
      call system_clock(begin_time,rat)
!==================
!read data in input
!==================

      !open(10,file='input')
      write(6,*)'What is the name of the system:'
      read(5,*)filename
      write(6,*)'What is the name of the trajecotry file:'
      read(5,*)pos_filename     
      write(6,*)'What is the name of the list file:'
      read(5,*)list_filename     
      write(6,*)'What is the total steps of the trajecotry:'
      read(5,*)nmo     !number of movie steps
      write(6,*)'What is the total number of atoms in the system:'
      read(5,*)nat     !number of atoms per mole.
      write(6,*)'What is the total number of water pairs:'
      read(5,*)np      !number of pairs   

      allocate(ndx_1(np))          
      allocate(ndx_2(np))          
      allocate(ndx_3(np))          

      if (nmo >= nmo_max)then
        write(6,*)'!! nmo> nmo_max !!'
        write(6,*)'!! stop !!'
        stop
      endif
      
      if (nat >= nmax)then
        write(6,*)'!! nat> nmax !!' 
        write(6,*)'!! stop !!'
        stop
      endif  

      list_filename=trim(list_filename)
      open(10,file=list_filename)     
      do k=1,np
          read(10,*)ndx_1(k),ndx_2(k),ndx_3(k)
      enddo
      close(10)

      allocate(x(nat,nmo))
      allocate(y(nat,nmo))
      allocate(z(nat,nmo))
      allocate(h(nmo))
      allocate(h_d(nmo))
      !allocate(hh(nmo))
      !allocate(hh_d(nmo))
      allocate(hb(np))         !Average H-bonded population 
      allocate(hb_d(np))

!=======================
!read in trajectory file 
!=======================
      open(10,file=trim(pos_filename))     
      do imovie=1,nmo
         read(10,*)                 !Neglect data of this line
         read(10,*)                 !Neglect data of this line
         do iatom= 1,nat
             read (10,*)atom_type,x(iatom,imovie),y(iatom,imovie),&
                        z(iatom,imovie)
         enddo
      enddo
      close(10)
      write(6,*) 'end of trajectory reading'
!=================关键部分=================
      allocate(corr_h(nmo))
      allocate(corr_h_d(nmo))
!==================================
!Calculate autocorrelation function
!==================================

      ! calculate <h(0)h(t)>/<h> and 
      ! calculate <h^d(0)h^d(t)>/<h^d> 
      ! Notice here <> is not average over
      ! different pairs of water molecules,
      ! but average over the time steps.
      do i=1, nmo
        corr_h(i)=0.0
        corr_h_d(i)=0.0
      enddo
      
      tot_hb=0.0
      do k=1,np
          hb(k)=0.0
          hb_d(k)=0.0
      enddo

      do k=1,np
      !calculate h(j)
        qj=0
        qj_d=0
        m1=ndx_1(k)
        m2=ndx_2(k)
        m3=ndx_3(k)
        
        do jj =1, nmo
          h(jj)=0
          h_d(jj)=0 
          r13= (x(m1,jj)-x(m3,jj))**2+       &
                    (y(m1,jj)-y(m3,jj))**2+  &
                    (z(m1,jj)-z(m3,jj))**2      !r:squra of distances
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
                   (z(m1,jj)-z(m2,jj)) 
          cosphi= pm/(sqrt(r23*r12))            !pm: point multiplication.
          if (r13 .lt. rohc .and. r12 .lt. rooc   & 
             .and. cosphi .gt. cosphic) then    
              h(jj)=1.0 
              qj=qj+h(jj)                          
          endif
          if (r13 .lt. rohc ) then                           
              h_d(jj)=1.0
              qj_d=qj_d+h_d(jj)                          
          endif
        enddo   
        qj=qj/nmo                            ! ave of hb for each pair 
        qj_d=qj_d/nmo
        hb(k)=qj
        hb_d(k)=qj_d
        tot_hb=tot_hb+hb(k)

        do mt=0,nmo-1     ! time interval
            if(hb(k)>hb_min) then
                scalar=0.d0
                do j=1, nmo-mt-1
                    scalar=scalar+h(j)*h(j+mt)     ! 1: the first pair of water molecules
                enddo
                scalar=scalar/(nmo-mt-1)         ! C_k(t)
                corr_h(mt+1)=corr_h(mt+1)+scalar  ! sum_C_k(t)
            endif
        enddo
        do mt=0,nmo-1                    ! time interval
            if(hb_d(k)>hb_min) then
                scalar_d=0.d0
                do j=1, nmo-mt-1
                    scalar_d=scalar_d+h(j)*h_d(j+mt)       ! 1: the first pair of water molecules
                enddo
                scalar_d=scalar_d/(nmo-mt-1)           ! C_k(t)
                corr_h_d(mt+1)=corr_h_d(mt+1)+scalar_d  ! sum_C_k(t)
            endif
        enddo

      enddo                                       ! k loop 
      tot_hb=tot_hb/np

      do mt=0,nmo-1     ! time interval
          corr_h(mt+1)=corr_h(mt+1)/(np*tot_hb)  
          corr_h_d(mt+1)=corr_h_d(mt+1)/(np*tot_hb)  
      enddo

      deallocate(x,y,z)
      deallocate(ndx_1,ndx_2,ndx_3)          
  
!===============================
!Write the result of HB dynamics
!===============================
!   open(10,file=trim(filename)//'_achbond_nhb.dat')

!      do k=1,np
!          if (hb(k)>hb_min) then
!              write(10,*)'#pair(k) ','Ave HBonded population(h_k) ',&
!                         'hb_d_k'
!              write(10,*) k,hb,hb_d
!              write(10,*)'#Step  ','time (ps)','Population(h_k(t))  ',&
!                         'h_d_k(t)'
!              do j =1,nmo
!                  write(10,*)j*delta_t,h(j),h_d(j) 
!              enddo
!          endif
!      enddo
!      write(6,*)'written in nhb_profile'
!      close(10)

!======================
!Write the correlation
!C_HB(t) and C_HB_d(t)     
!======================
      open(10,file=trim(filename)//'_Chandra_achbond_h.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,corr_h(i)
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_Chandra_achbond_h.dat'
      close(10)

      open(10,file=trim(filename)//'_Chandra_achbond_h_d.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,corr_h_d(i)
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_Chandra_achbond_h_d.dat'
      close(10)

!==========================
!Write the correlation
!ln(C_HB(t)) and lnC_HB_d(t)     
!==========================
      open(10,file=trim(filename)//'_Chandra_achbond_ln_h.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,log(corr_h(i))
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_Chandra_achbond_ln_h.dat'
      close(10)

      open(10,file=trim(filename)//'_Chandra_achbond_ln_h_d.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,log(corr_h_d(i))
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_Chandra_achbond_ln_h_d.dat'
      close(10)
!==========================
! print <h>      
      open(10,file=trim(filename)//'_Chandra_average_h.dat')
        write(10,*) '<h>:',tot_hb
        write(6,*)'written in '//trim(filename)//&
                  '_Chandra_average_h.dat'
      close(10)

      deallocate (h,h_d,corr_h,corr_h_d)
!==============================================================
      call system_clock(end_time,rat)
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END

