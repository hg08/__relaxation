module hb1_module
      implicit none
      integer,parameter :: rk = 4
      real,public,parameter :: rooc=12.25                ! cutoff distance of rOO (3.5**2 )
      real,public,parameter :: rohc=6.0025                ! rOH (2.45**2)
      real,public,parameter :: cosphic=0.866             ! 1.732/2; phiC=pi/6.
      real(kind=rk),public,parameter :: h_min=0.5 ! condition for the existence of a h-bond for a step
      real(kind=rk),public,parameter :: hb_min=0.5 ! condition for the existence of h-bond for a pair of water molecules
      real(kind=rk),allocatable,dimension (:)    :: h,hb,corr_h
      real(kind=rk),allocatable,dimension (:)    :: hh,corr_hh
      real(kind=rk),public :: r12,r13,r23,cosphi,pm,qj, tot_hb
end module hb1_module
