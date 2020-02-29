!2015/09/9
!==============
!Calculate k(t) 
!==============
     ! input file: 
     ! system name 
     ! trajectory name
     ! list name
     ! nmo
     ! nat
     ! number of config. of molecules, or np 
!============================================
      program rfhb_kpair
              use hb1_module
      implicit none
!==========
!parameters
!==========
      character(LEN=30) :: filename ,pos_filename,list_filename         
      real(kind=rk),parameter :: rate=0.80        ! condition for cutting off autocorrelation functions
      real(kind=rk),parameter :: delta_t=0.0005   ! fs
      integer :: begin_time,end_time,rat,&
                        i,j,k,jj,nmo,nat,iatom,& 
                        imovie,np,m1,m2,m3,mt 
      real(kind=rk),allocatable,dimension (:)    :: dh
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3)  :: atom_type  ! we are not interested in the atom type in this calculation,
                                      ! thus I do not use allocatable array for it.
      integer,allocatable,dimension(:)           :: ndx_1, ndx_2, ndx_3
      real(kind=rk)  :: scalar
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
      read(5,*)nmo!number of movie steps
      write(6,*)'What is the total number of atoms in the system:'
      read(5,*)nat!number of atoms per mole.
      write(6,*)'What is the total number of water pairs:'
      read(5,*)np!number of pairs   

      allocate(ndx_1(np))          
      allocate(ndx_2(np))          
      allocate(ndx_3(np))          
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
      allocate(hb(np))!Average H-bonded population 
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
         enddo
      enddo
      close(10)
      write(6,*) 'end of trajectory reading'
!==================================
!Calculate autocorrelation function
!==================================
! calculate <h(0)h(t)>/<h>  
! Notice here <> is not average over
! different pairs of water molecules,
! but average over the time steps.
      allocate(corr_h(nmo))
      !do i=1, nmo-1
      corr_h(:)=0.0
      !enddo
      tot_hb=0.0
      !do k=1,np
      hb(:)=0.0
      !enddo

      do k=1,np
        qj=0
        m1=ndx_1(k)
        m2=ndx_2(k)
        m3=ndx_3(k)
        !calculate h(j)
        do jj =1, nmo
          h(jj)=0
          r13= (x(m1,jj)-x(m3,jj))**2+       &
                    (y(m1,jj)-y(m3,jj))**2+  &
                    (z(m1,jj)-z(m3,jj))**2!r:squra of distances
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
          cosphi= pm/(sqrt(r23*r12))!pm: point multiplication.
          if (r13 .lt. rohc .and. r12 .lt. rooc   & 
             .and. cosphi .gt. cosphic) then    
              h(jj)=1.0 
              qj=qj+h(jj)                          
          endif
        enddo   
        do jj=1,nmo
            if (jj<nmo-1 .and. jj>0)then
                dh(jj)=3*h(jj)-4*h(jj+1)+h(jj+2)!threepoint formula; -1*(increament) of h
                dh(jj)=dh(jj)/2
            else
                dh(jj)=h(jj-1)-h(jj)               
            endif    
        enddo  
        qj=qj/nmo! ave of hb for each pair 
        hb(k)=qj
        tot_hb=tot_hb+hb(k)
        do mt=0,nmo-1! time interval
            if(hb(k)>hb_min) then
                scalar=0.d0
                do j=1, nmo-mt
                    scalar=scalar+dh(j)*(1-h(j+mt))! 1: the first pair of water molecules
                enddo
                scalar=scalar/(nmo-mt)                ! C_k(t)
                corr_h(mt+1)=corr_h(mt+1)+scalar        ! sum_C_k(t)
            endif
        enddo
      enddo! k loop 
      tot_hb=tot_hb/np
      do mt=0,nmo-1! time interval
          corr_h(mt+1)=corr_h(mt+1)/(np*tot_hb*delta_t)  
      enddo
      write(6,*) corr_h(1),corr_h(2),corr_h(3)
      deallocate(ndx_1,ndx_2,ndx_3,hb)          
!========================================================================
!calculate k(t)  
!Notice that we start from 'i=2', instead of 'i=1' !
!We should not include the first term! Since it related to forming of HB! 
!This is not consistant with our assumption that the HB is already formed!
!========================================================================
      open(10,file=trim(filename)//'_rfhb_h.dat')
        do i=2,int(nmo*rate)                
            write(10,*)(i-1)*delta_t,corr_h(i)
        enddo
        write(6,*)'written in rfhb_h'
      close(10)
!==============
!write ln(k(t))  
!==============
      open(10,file=trim(filename)//'_rfhb_ln_h.dat')
        do i=2,int(nmo*rate)                
            write(10,*)(i-1)*delta_t,log(corr_h(i))
        enddo
        write(6,*)'written in rfhb_ln_h'
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
