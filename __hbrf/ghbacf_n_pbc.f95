 program ghbacf_n_pbc
 USE tools_0
 USE tools_1
 USE tools_2
 implicit none
!2015/09/9
!2020-2-16: k(t) is devided by dleta_t=(ns*delta_t0)
!==============
!Calculate k(t) 
!==============
  ! input file:
  ! time step
  ! system name 
  ! trajectory name
  ! list name
  ! nmo_start
  ! nmo_end
  ! nat
  ! number of config. of molecules, or np 
  ! new time step
!==========
!parameters
!==========
character(LEN=200) :: filename,pos_filename,list_filename      
integer,parameter :: rk=4   
INTEGER,PARAMETER :: n_neighbors = 5 ! used for do window averaging for k(t)
real(kind=rk),parameter :: rate=0.80 ! for cutting off autocorrelation functions
real(kind=rk),parameter :: rooc=12.25 ! cutoff distance square of rOO(3.5 A )
real(kind=rk),parameter :: rohc=6.0025 ! square of rOH (2.45 A)
real(KIND=rk),parameter :: cosphic=0.866             ! 1.732/2; phiC=pi/6.
real(kind=rk),parameter :: hb_min=0.00001 ! condition for the existence of h-bond
real(kind=rk)  :: r12,r13,r23,cosphi,pm,qj,tot_hb,delta_t0,delta_t
integer :: begin_time,end_time,rat,i,j,k,jj,&
           nmo,nat,np,m1,m2,m3,mt,ns,nmo_start,nmo_end 
real(kind=rk),allocatable,dimension (:)    :: h
real(kind=rk),allocatable,dimension (:)    :: h_d 
real(kind=rk),allocatable,dimension (:)    :: hb
real,allocatable,dimension (:,:)           :: x,y,z
character(LEN=3) :: atom_type  
integer,allocatable,dimension(:) :: ndx_1, ndx_2, ndx_3
real(kind=rk),allocatable,dimension(:) :: corr_n ! the correlation function n(t)
real(kind=rk) :: scalar_n ! for n(t)

! For PBC
real(kind=rk),parameter :: a=31.00       ! BOX size a 
real(kind=rk),parameter :: b=15.64       ! BOX size b
real(kind=rk),parameter :: c=15.64       ! BOX size c
REAL(KIND=rk) :: distance2,diff_axis

call system_clock(begin_time,rat)

!==================
!read data in input
!==================
write(6,*)'What is the timestep (ps):'
read(5,*)delta_t0
write(6,*)'What is the name of the system:'
read(5,*)filename
write(6,*)'What is the name of the trajecotry file:'
read(5,*)pos_filename     
write(6,*)'What is the name of thei list file:'
read(5,*)list_filename     
write(6,*)'What is the initial step of the trajecotry:'
read(5,*)nmo_start !number of the first movie steps
write(6,*)'What is the end step of the trajecotry:'
read(5,*)nmo_end !number of the last movie steps
write(6,*)'What is the total number of atoms in the system:'
read(5,*)nat!number of atoms per mole.
write(6,*)'What is the total number of water pairs:'
read(5,*)np!number of pairs   
write(6,*)'What is the time step for calculating the correlation:'
read(5,*)ns! [ns*0.0005] ps is the new time step for calculating correl,if delta_t0=0.0005 ps.

allocate(ndx_1(np))          
allocate(ndx_2(np))          
allocate(ndx_3(np))          
! Read list file
CALL read_index_file(list_filename,np,ndx_1,ndx_2,ndx_3)

delta_t=ns*delta_t0
nmo=nmo_end-nmo_start
nmo=nmo/ns
allocate(x(nat,nmo))
allocate(y(nat,nmo))
allocate(z(nat,nmo))
allocate(h(nmo))
allocate(h_d(nmo))
allocate(hb(np))!Average H-bonded population 

! Initialization
h(:) = 0
h_d(:) = 0
hb(:) = 0

!=======================
!read in trajectory file 
!=======================
! call en explicit interface
CALL read_trajectory(pos_filename,nat,nmo,nmo_start,ns,atom_type,x,y,z)

!==================================
!Calculate autocorrelation function
!==================================
! calculate <h(0)(1-h(t))[h^d](t)>/<h>  
! Notice here <> is not average over
! different pairs of water molecules,
! but average over the time steps.
allocate(corr_n(nmo))
tot_hb=0.0
corr_n=0.0

do k=1,np
  qj=0
  m1=ndx_1(k)
  m2=ndx_2(k)
  m3=ndx_3(k)
  !calculate h(j)
  do jj =1, nmo
    h(jj)=0
    r13 = distance2(x(m1,jj),y(m1,jj),z(m1,jj), &
      x(m3,jj),y(m3,jj),z(m3,jj),a,b,c) ! r12,r13,r23:square of distances 
    r12 = distance2(x(m1,jj),y(m1,jj),z(m1,jj), &
          x(m2,jj),y(m2,jj),z(m2,jj),a,b,c) ! r12,r13,r23:square of distances 
    r23 = distance2(x(m2,jj),y(m2,jj),z(m2,jj), &
          x(m3,jj),y(m3,jj),z(m3,jj),a,b,c) ! r12,r13,r23:square of distances 
    pm= diff_axis(x(m3,jj),x(m2,jj),a)*        &
          diff_axis(x(m1,jj),x(m2,jj),a)+      & 
          diff_axis(y(m3,jj),y(m2,jj),b)*      & 
          diff_axis(y(m1,jj),y(m2,jj),b)+      & 
          diff_axis(z(m3,jj),z(m2,jj),c)*      &
          diff_axis(z(m1,jj),z(m2,jj),c)       ! pm: point multiplication.
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

  !do jj=1,nmo
  !    if (jj >2) then
  !        dh(jj)=(3*h(jj)-4*h(jj-1)+h(jj-2))/delta_t !threepoint formula; -1*(increament) of h
  !        !dh(jj)=2*dh(jj)
  !    else
  !        dh(jj)=(h(jj+1)-h(jj))/delta_t
  !    endif
  !enddo  
  qj=qj/nmo ! ave of hb for each pair 
  hb(k)=qj
  tot_hb=tot_hb+hb(k)
  do mt=0,nmo-1 ! time interval
      if(hb(k)>hb_min) then
          scalar_n=0.d0
          do j=1, nmo-mt
              !scalar=scalar-dh(j)*(1-h(j+mt))*(h_d(j+mt)) ! 1: the first pair of water molecules
              scalar_n=scalar_n + h(j)*(1-h(j+mt))*h_d(j+mt)   ! 1: For calculating n(t)
          enddo
          !scalar=scalar/(nmo-mt) ! k_k(t)
          scalar_n=scalar_n/(nmo-mt)
          !corr(mt+1)=corr(mt+1)+scalar ! sum_k_k(t)
          corr_n(mt+1)=corr_n(mt+1)+scalar_n
      endif
  enddo
enddo ! k loop 

!=====================
!Nromalization of k(t)
!=====================
do mt=0,nmo-1! time interval
    !corr(mt+1)=corr(mt+1)/tot_hb
    corr_n(mt+1)=corr_n(mt+1)/tot_hb
enddo
!=============
!Smoothing n(t)
!=============
CALL moving_window_average(corr_n,n_neighbors)

!========
! TESTING
!========
!write(6,*) corr(1),corr(2),corr(3)

!==================
! Deallocate arrays
!==================
deallocate(ndx_1,ndx_2,ndx_3,hb)          
  
!===============
! write out n(t)
!===============
open(10,file=trim(filename)//'_n.dat')
    do i=1,int(nmo*rate)                
        write(10,*)(i-1)*delta_t,corr_n(i)
    enddo
    write(6,*)'written in'//trim(filename)//'_n.dat'
close(10)
!====================
! write out log[n(t)]
!====================
open(10,file=trim(filename)//'_log_n.dat')
    do i=1,int(nmo*rate)                
        write(10,*)(i-1)*delta_t,log(corr_n(i))
    enddo
    write(6,*)'written in'//trim(filename)//'_log_n.dat'
close(10)
!==============================================================
deallocate(x,y,z)
deallocate (h,corr_n,h_d)

! Count time
call system_clock(end_time,rat)
write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 

END PROGRAM ghbacf_n_pbc 

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
