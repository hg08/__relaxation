     !2020/2/15
     !===================================================
     !For ONE pair of molecules，we can obtain C_HB(t):
     !   C_{HB}(t)= <h(0)h(t)>/<h>
     !1)We also implementate these correlation functions by
     !calculate the average over N pairs of molecules. 
     !2)For any time step, instead of the original time 
     ! step.      
     !3)Newly added function of this code: We considered 
     ! the PBC in an easier way.
     !===================================================      
     ! input file: input
     ! time step
     ! name of system
     ! name of trajectory
     ! name of list
     ! nmo
     ! nat
     ! number of pairs of molecules(np) 
     ! the new time step, or ns
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
      real,parameter :: rohc=6.0025               ! rOH (2.45**2)
      real,parameter :: cosphic=0.866             ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: h_min=0.5 ! condition for the existence of a h-bond for a step
      real(kind=rk),parameter :: hb_min=0.5 ! condition for the existence of h-bond for a pair of water molecules


      real(kind=rk)           :: r12,r13,r23,cosphi,pm,qj,&
                                 tot_hb,delta_t,delta_t0 
      integer :: begin_time,end_time,rat,&
                 i,j,k,jj,nmo,nat,iatom,& 
                 imovie,np,m1,m2,m3,mt,&
                 nqj,tot_nhb,n_bonded_pairs,ns 
      real(kind=rk),allocatable,dimension (:)    :: h,hb,corr_h
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3)  :: atom_type    ! we are not interested in the atom type in this calculation,
                                        ! thus I do not use allocatable array for it.
      integer,allocatable,dimension(:)           :: ndx_1,ndx_2,ndx_3,&
          nhb_exist
      real(kind=rk)  :: scalar 
      logical,allocatable,dimension (:)  :: hb_exist
      
      ! For PBC
      real(kind=rk),parameter :: a=31.00        ! condition for cutting off autocorrelation functions
      real(kind=rk),parameter :: b=15.64       ! condition for cutting off autocorrelation functions
      real(kind=rk),parameter :: c=15.64       ! condition for cutting off autocorrelation functions
      REAL(KIND=rk) :: distance2,diff_axis

      call system_clock(begin_time,rat)
     !==================
     !read data in input
     !==================
      write(6,*)'What is the time step in the traj. file (ps):'
      read(5,*)delta_t0
      write(6,*)'What is the name of the system:'
      read(5,*)filename
      write(6,*)'What is the name of the trajecotry file:'
      read(5,*)pos_filename     
      write(6,*)'What is the name of thei list file:'
      read(5,*)list_filename     
      write(6,*)'What is the total steps of the trajecotry:'
      read(5,*)nmo    ! number of movie steps
      write(6,*)'What is the total number of atoms in the system:'
      read(5,*)nat    ! number of atoms per mole.
      write(6,*)'What is the total number of water pairs:'
      read(5,*)np     ! number of pairs   
      write(6,*)'What is the time step for calculating CORRELATION:' 
      read(5,*)ns     ! [ns*0.0005] ps is the new time step for calculating correl func.

      list_filename=trim(list_filename)
      allocate(ndx_1(np))          
      allocate(ndx_2(np))          
      allocate(ndx_3(np))          
      open(10,file=list_filename)     
      do k=1,np
          read(10,*)ndx_1(k),ndx_2(k),ndx_3(k)
      enddo
      close(10)

      delta_t=ns*delta_t0    ! unit: ps
      nmo=nmo/ns    ! Length of the correl. function 
      allocate(x(nat,nmo))
      allocate(y(nat,nmo))
      allocate(z(nat,nmo))
      allocate(h(nmo))
      allocate(hb(np))    ! Average H-bonded population 
      allocate(nhb_exist(np))
     !=======================
     !read in trajectory file 
     !=======================
      open(10,file=trim(pos_filename))     
      do imovie=1,nmo
         read(10,*)    !Neglect data of this line
         read(10,*)                 
         do iatom= 1,nat
             read (10,*)atom_type,x(iatom,imovie),y(iatom,imovie),&
                        z(iatom,imovie)
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
      ! loop
      corr_h(:)=0.0
      
      tot_hb=0.0
      tot_nhb=0
      
      ! loop
      hb(:)=0.0
      nhb_exist(:)=0
     !=============
     !The main loop
     !=============      
      do k=1,np
        qj=0
        nqj=0
        m1=ndx_1(k)
        m2=ndx_2(k)
        m3=ndx_3(k)
        ! Calculate h(j)
        do jj =1, nmo
          h(jj)=0.0
          hb_exist(jj)=.False.
          r13 = distance2(x(m1,jj),y(m1,jj),z(m1,jj), &
                          x(m3,jj),y(m3,jj),z(m3,jj),a,b,c)! r12,r13,r23:square of distances 
          r12 = distance2(x(m1,jj),y(m1,jj),z(m1,jj), &
                          x(m2,jj),y(m2,jj),z(m2,jj),a,b,c)! r12,r13,r23:square of distances 
          r23 = distance2(x(m2,jj),y(m2,jj),z(m2,jj), &
                          x(m3,jj),y(m3,jj),z(m3,jj),a,b,c)! r12,r13,r23:square of distances 
          pm= diff_axis(x(m3,jj),x(m2,jj),a)*           &
                   diff_axis(x(m1,jj),x(m2,jj),a)+      & 
                   diff_axis(y(m3,jj),y(m2,jj),b)*      & 
                   diff_axis(y(m1,jj),y(m2,jj),b)+      & 
                   diff_axis(z(m3,jj),z(m2,jj),c)*      &
                   diff_axis(z(m1,jj),z(m2,jj),c)       ! pm: point multiplication. 
          cosphi= pm/(sqrt(r23*r12))
          if (r13 .lt. rohc .and. r12 .lt. rooc  & 
             .and. cosphi .gt. cosphic) then    
              h(jj)=1.0 
              hb_exist(jj)=.True.
              qj=qj+h(jj)    ! To calculate ave population of HB over all starting points for one pair of water molecules.                          
              nqj=nqj+1
          endif
        enddo   
        hb(k)=qj 
        nhb_exist(k)=nqj
        tot_hb=tot_hb+hb(k)
        tot_nhb=tot_nhb+nhb_exist(k)
        !==========================================
        !Calcualte the correlation function C_HB(t)
        !==========================================
        if (hb(k)>hb_min) then
            do mt=0,nmo-1    ! time interval
                scalar=0.d0
                do j=1,nmo-mt-1
                    scalar=scalar+h(j)*h(j+mt)  
                    enddo
                corr_h(mt+1)=corr_h(mt+1)+scalar    ! sum_C_k(t)
            enddo
        endif
      enddo    ! End of k-loop 
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
      open(10,file=trim(filename)//'_hbacf_h.dat')
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
      open(10,file=trim(filename)//'_hbacf_log_h.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,log(corr_h(i))
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_hbacf_log_h.dat'
      close(10)
     !===========
     ! Print <h>      
     !===========      
      open(10,file=trim(filename)//'_ave_h.dat')
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
      END PROGRAM

      !==============
      ! The functions 
      !==============
      REAL(KIND=4) FUNCTION direct_distance2(u1,v1,w1,u2,v2,w2)
          INTEGER, PARAMETER :: rk=4
          real(kind=rk) :: u1,v1,w1,u2,v2,w2
          direct_distance2 = (u2-u1)**2 + (v2-v1)**2 + (w2-w1)**2
      END FUNCTION direct_distance2

      REAL(KIND=4) FUNCTION distance2(u1,v1,w1,u2,v2,w2,a,b,c)
          INTEGER, PARAMETER :: rk=4
          real(kind=rk) :: u1,v1,w1,u2,v2,w2,a,b,c
          logical :: A1,A2,A3,B1,B2,B3,C1,C2,C3
          A1 = (abs(u1-u2) > a/2 .AND. u1 > u2)
          A2 = (abs(u1-u2) > a/2 .AND. u2 > u1)
          A3 = (abs(u1-u2) < a/2)
          B1 = (abs(v1-v2) > b/2 .AND. v1 > v2)
          B2 = (abs(v1-v2) > b/2 .AND. v2 > v1)
          B3 = (abs(v1-v2) < b/2)
          C1 = (abs(w1-w2) > c/2 .AND. w1 > w2)
          C2 = (abs(w1-w2) > c/2 .AND. w2 > w1)
          C3 = (abs(w1-w2) < c/2)
          if (A3 .and. B3 .and. C3) then
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B3 .and. C3) then
              u2 = u2-a
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B3 .and. C3) then
              u1 = u1-a
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B1 .and. C3) then
              v1 = v1 - b
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B3 .and. C2) then
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B2 .and. C3) then
              v2 = v2 - b
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B3 .and. C1) then
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B1 .and. C3) then
              u1 = u1 - a
              v1 = v1 - b
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B2 .and. C3) then
              u1 = u1 - a
              v2 = v2 - b
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B3 .and. C1) then
              u1 = u1 - a
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B3 .and. C2) then
              u1 = u1 - a
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B1 .and. C3) then
              u2 = u2 - a
              v1 = v1 - b
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B2 .and. C3) then
              u2 = u2 - a
              v2 = v2 - b
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B3 .and. C1) then
              u2 = u2 - a
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B3 .and. C2) then
              u2 = u2 - a
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B1 .and. C1) then
              v1 = v1 - b
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B1 .and. C2) then
              v1 = v1 - b
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B2 .and. C1) then
              v2 = v2 - b
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A3 .and. B2 .and. C2) then
              v2 = v2 - b
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B1 .and. C1) then
              u1 = u1 - a
              v1 = v1 - b
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B1 .and. C2) then
              u1 = u1 - a
              v1 = v1 - b
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B2 .and. C1) then
              u1 = u1 - a
              v2 = v2 - b
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A1 .and. B2 .and. C2) then
              u1 = u1 - a
              v2 = v2 - b
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B1 .and. C1) then
              u2 = u2 - a
              v1 = v1 - b
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B1 .and. C2) then
              u2 = u2 - a
              v1 = v1 - b
              w2 = w2 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B2 .and. C1) then
              u2 = u2 - a
              v2 = v2 - b
              w1 = w1 - c
              distance2= direct_distance2(u1,v1,w1,u2,v2,w2) 
          elseif (A2 .and. B2 .and. C2) then
              u2 = u2 - a
              v2 = v2 - b
              w2 = w2 - c
              distance2 = direct_distance2(u1,v1,w1,u2,v2,w2) 
          endif
      END FUNCTION distance2

      REAL(KIND=4) FUNCTION diff_axis(u1,u2,a)
          logical :: A1,A2,A3
          INTEGER, PARAMETER :: rk=4
          real(kind=rk) :: u1, u2, a
          A1 = (abs(u1-u2) > a/2 .AND. u1 > u2)
          A2 = (abs(u1-u2) > a/2 .AND. u2 > u1)
          A3 = (abs(u1-u2) < a/2)
          if (A3) then
              diff_axis= u1 - u2
          elseif (A1) then
              u1 = u1 - a
              diff_axis= u1-u2
          elseif (A2) then
              u2 = u2 - a
              diff_axis= u1 - u2
          endif
      END FUNCTION diff_axis
