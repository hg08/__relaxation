!2017/06/17
!==========================================
!Calculate k_in(t) and its intrgration n(t)
!==========================================
     ! input file: 
     ! system name 
     ! trajectory name
     ! list name
     ! nmo
     ! nat
     ! number of config. of molecules, or np 
!============================================
      program rfhb_in_kpair
      implicit none
!==========
!parameters
!==========
      character(LEN=30) :: filename ,pos_filename,list_filename         
      integer,parameter :: rk=4              
      real(kind=rk),parameter :: rate=0.80        ! condition for cutting off autocorrelation functions
      real,parameter :: rooc=12.25                ! cutoff distance of rOO (3.5**2 )
      real,parameter :: rohc=6.0025               ! rOH (2.45**2)
      real,parameter :: cosphic=0.866             ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: hb_min=0.0000001 ! condition for the existence of h-bond
      real(kind=rk),parameter :: delta_t=0.0005   ! ps
      real(kind=rk)           :: r12,r13,r23,cosphi,pm,qj,&
                                 tot_hb
      integer :: begin_time,end_time,rat,&
                        i,j,k,jj,nmo,nat,iatom,& 
                        imovie,np,m1,m2,m3,mt 
      real(kind=rk),allocatable,dimension (:)    :: h,h_d,dh,prob_in_h_d
      real(kind=rk),allocatable,dimension (:)    :: hb
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3)  :: atom_type  ! we are not interested in the atom type in this calculation,
                                      ! thus I do not use allocatable array for it.
      integer,allocatable,dimension(:)           :: ndx_1, ndx_2, ndx_3
      real(kind=rk),allocatable,dimension (:)    :: corr_h
      real(kind=rk)  :: scalar

      ! For PBC
      real(kind=rk),parameter :: a=31.00       ! condition for cutting off autocorrelation functions
      real(kind=rk),parameter :: b=15.64       ! condition for cutting off autocorrelation functions
      real(kind=rk),parameter :: c=15.64       ! condition for cutting off autocorrelation functions
      REAL(KIND=rk) :: distance2,diff_axis
      REAL(KIND=rk) :: direct_distance2,direct_diff_axis

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
      allocate(h_d(nmo))
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
      allocate(prob_in_h_d(nmo))
      !do i=1, nmo-1
      corr_h(:)=0.0
      prob_in_h_d(:)=0.0
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
          r13 = direct_distance2(x(m1,jj),y(m1,jj),z(m1,jj), &
                          x(m3,jj),y(m3,jj),z(m3,jj),a,b,c)! r12,r13,r23:square of distances 
          r12 = direct_distance2(x(m1,jj),y(m1,jj),z(m1,jj), &
                          x(m2,jj),y(m2,jj),z(m2,jj),a,b,c)! r12,r13,r23:square of distances 
          r23 = direct_distance2(x(m2,jj),y(m2,jj),z(m2,jj), &
                          x(m3,jj),y(m3,jj),z(m3,jj),a,b,c)! r12,r13,r23:square of distances 
          pm= direct_diff_axis(x(m3,jj),x(m2,jj),a)*           &
                   direct_diff_axis(x(m1,jj),x(m2,jj),a)+      & 
                   direct_diff_axis(y(m3,jj),y(m2,jj),b)*      & 
                   direct_diff_axis(y(m1,jj),y(m2,jj),b)+      & 
                   direct_diff_axis(z(m3,jj),z(m2,jj),c)*      &
                   direct_diff_axis(z(m1,jj),z(m2,jj),c)       ! pm: point multiplication. 
          cosphi= pm/(sqrt(r23*r12))!pm: point multiplication.
          if (r13 .lt. rohc .and. r12 .lt. rooc   & 
             .and. cosphi .gt. cosphic) then    
              h(jj)=1.0 
              qj=qj+h(jj)                          
          endif
          if (r12 .lt. rooc) then    
              h_d(jj)=1.0 
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
                    scalar=scalar+dh(j)*(1-h(j+mt))*h_d(j+mt)! 1: the first pair of water molecules
                enddo
                scalar=scalar/(nmo-mt)                ! C_k(t)
                corr_h(mt+1)=corr_h(mt+1)+scalar        ! sum_C_k(t)
            endif
        enddo
      enddo! k loop 
      tot_hb=tot_hb/np
      do mt=0,nmo-1! time interval
          corr_h(mt+1)=corr_h(mt+1)/(np*tot_hb*delta_t)  

          scalar=0.d0  ! Calculate n(t)=\int_0^t dt' k_in(t'), donoted
          ! by prob_in_h_d
          do jj=1, nmo-1
              scalar=scalar+corr_h(jj)*delta_t
          enddo
          prob_in_h_d(mt+1)=scalar
      enddo
      write(6,*) corr_h(1),corr_h(2),corr_h(3)
      deallocate(ndx_1,ndx_2,ndx_3,hb)          
!========================================================================
!calculate k_in(t)  
!Notice that we start from 'i=2', instead of 'i=1' !
!We should not include the first term! Since it related to forming of HB! 
!This is not consist with our assumption that the HB is already formed!
!========================================================================
      open(10,file=trim(filename)//'_rfhb_in_h--n.dat')
        do i=2,int(nmo*rate)                
            write(10,*)'# t (ps)    ', '    h(t)    ','    n(t)    '
            write(10,*)(i-1)*delta_t,corr_h(i), prob_in_h_d(i)
        enddo
        write(6,*)'written in '//trim(filename)//'_rfhb_in_h--n.dat'
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
      deallocate (h,dh,corr_h,h_d)
      call system_clock(end_time,rat)
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
!=====================
!write log for k_in(t)
!=====================
      open(10,file=trim(filename)//'_k_in(t).log')
      write(10, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      write(6,*)'written in rfhb.log'
      close(10)
      END PROGRAM

      !==============
      ! The functions 
      !==============
      REAL(KIND=4) FUNCTION direct_distance2(u1,v1,w1,u2,v2,w2,a,b,c)
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
      REAL(KIND=4) FUNCTION direct_diff_axis(u1,u2,a)
          INTEGER, PARAMETER :: rk=4
          real(kind=rk) :: u1, u2, a
              direct_diff_axis= u1 - u2
      END FUNCTION direct_diff_axis
