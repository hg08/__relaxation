!2015/09/17
!============================================
!对指定的1对原子，可求其氢键布居算符之
!平均值和其自关联函数S_HB(t),C_HB(t)。
! 用此方法，可得 N 对水分子之氢键布居随
! 时间之改变：k(t)。 
     ! input file: 
     ! system name 
     ! trajectory name
     ! list name
     ! nmo
     ! nat
     ! number of config. of molecules, or np 
!============================================
      program rfachbond_1pair
      implicit none

!==========
!parameters
!==========
      character(LEN=30) :: filename ,pos_filename,list_filename         
      integer,parameter :: rk=4              
      integer,parameter :: nmax=775               ! max number of atoms
      integer,parameter :: nmo_max=300000         ! max number of movie
      real(kind=rk),parameter :: rate=0.50        ! condition for cutting off autocorrelation functions
      real,parameter :: rooc=29.16                ! cutoff distance of rOO (5.4**2 )
      real,parameter :: rohc=12.25                ! rOH (3.5**2)
      real,parameter :: cosphic=0.866             ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: hb_min=0.0000001 ! condition for the existence of h-bond
      real(kind=rk)  :: delta_t                   ! fs
      real(kind=rk)           :: r12,r13,r23,cosphi,pm,qj,&
                                 tot_hb
      integer :: begin_time,end_time,rat,&
                        i,j,k,jj,nmo,nat,iatom,& 
                        imovie,np,m1,m2,m3,mt 
      real(kind=rk),allocatable,dimension (:)    :: h,dh
      real(kind=rk),allocatable,dimension (:)    :: hb
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3)  :: atom_type  ! we are not interested in the atom type in this calculation,
                                      ! thus I do not use allocatable array for it.
      integer,allocatable,dimension(:)           :: ndx_1, ndx_2, ndx_3
      real(kind=rk),allocatable,dimension (:)    :: corr_h
      real(kind=rk)  :: scalar
      call system_clock(begin_time,rat)
!==================
!read data in input
!==================

      !open(10,file='input')
      write(6,*)'What is the timestep of the trajectory:'
      read(5,*) delta_t
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
      allocate(dh(nmo))
      !allocate(hh(nmo))
      allocate(hb(np))         !Average H-bonded population 

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
!===========关键部分===============
      allocate(corr_h(nmo))
!==================================
!Calculate autocorrelation function
!==================================

      ! calculate <h(0)h(t)>/<h>  
      ! Notice here <> is not average over
      ! different pairs of water molecules,
      ! but average over the time steps.

      do i=1, nmo-1
        corr_h(i)=0.0
      enddo
      
      tot_hb=0.0
      do k=1,np
          hb(k)=0.0
      enddo

      do k=1,np
      !calculate h(j)
        qj=0
        m1=ndx_1(k)
        m2=ndx_2(k)
        m3=ndx_3(k)
        
        do jj =1, nmo
          h(jj)=0
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
        enddo   
        do jj=1,nmo
            if (jj<nmo .and. jj>0)then
                dh(jj)=h(jj)-h(jj+1)             ! -1*(increament) of h 
            else
                dh(jj)=h(jj-1)-h(jj)               
            endif    
        enddo  
        qj=qj/nmo                           ! ave of hb for each pair 
        hb(k)=qj
        tot_hb=tot_hb+hb(k)

        do mt=0,nmo-1                       ! time interval
            if(hb(k)>hb_min) then
                scalar=0.d0
                do j=1, nmo-mt
                    scalar=scalar+dh(j)*(1-h(j+mt))     ! 1: the first pair of water molecules
                enddo
                scalar=scalar/(nmo-mt)                ! C_k(t)
                corr_h(mt+1)=corr_h(mt+1)+scalar        ! sum_C_k(t)
            endif
        enddo

      enddo                                       ! k loop 
      tot_hb=tot_hb/np

      do mt=0,nmo-1                               ! time interval
          corr_h(mt+1)=corr_h(mt+1)/(np*tot_hb*delta_t)  
      enddo
      write(6,*) corr_h(1),corr_h(2),corr_h(3)
      deallocate(ndx_1,ndx_2,ndx_3,hb)          
  
!======================
!calculate k_HB(t)  
         ! Notice that we start from 'i=2', instead of 'i=1' !
         ! We should not include the first term! Since it related to forming of HB! 
         ! This is not consist with our assumption that the HB is already formed!
!======================
      open(10,file=trim(filename)//'_rfachb_h.dat')
        do i=2,int(nmo*rate)                
            write(10,*)(i-1)*delta_t,corr_h(i)
        enddo
        write(6,*)'written in rfachb_h'
      close(10)

!======================      
!write ln(k_HB(t))  
      open(10,file=trim(filename)//'_rfachb_ln_h.dat')
        do i=2,int(nmo*rate)                
            write(10,*)(i-1)*delta_t,log(corr_h(i))
        enddo
        write(6,*)'written in rfachb_ln_h'
      close(10)

!=============
!test the traj
!=============
!      open(10,file='x_pos1.dat')
!        do i=1,int(nmo)
!            write(10,*)i,x(10,i),x(11,i)
!        enddo
!        write(6,*)'written in x_pos1.dat'
!      close(10)
!      open(10,file='x_pos2.dat')
!        do i=1,int(nmo)
!            write(10,*)i,x(12,i),x(13,i)
!        enddo
!        write(6,*)'written in x_pos2.dat'
!      close(10)
!
!      deallocate(x,y,z)
!      deallocate (h,dh,corr_h)
!==============================================================
      call system_clock(end_time,rat)
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
