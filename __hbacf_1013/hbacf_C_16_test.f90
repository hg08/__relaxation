!2015/09/13
!=====================================
!对指定的 1 对原子，可求其氢键布居算符之
!平均值和其自关联函数S_HB(t),S_d_HB(t),
!C_HB(t),C_d_HB(t)。
!      
! C_{HB}(t)= <h(0)h(t)>/<h>
! C^(d)_{HB}(t)= <h(0)h^(d)(t)>/<h>
!
! 用此方法，可得 N 对水分子之氢键布居算
! 符. 
     ! input file: input
     ! name of system
     ! name of trajectory
     ! name of list
     ! nmo
     ! nat
     ! number of pairs of molecules(np) 
!=====================================
      program bhbaf
      implicit none
!==========
!parameters
!==========
      character(LEN=30) :: filename ,pos_filename,list_filename         
      integer,parameter :: rk=4              
      integer,parameter :: nmax=775               ! max number of atoms
      integer,parameter :: nmo_max=300000         ! max number of movie
      real(kind=rk),parameter :: rate=0.50        ! condition for cutting off autocorrelation functions
      !real,parameter :: rooc=29.16               ! cutoff distance of rOO (5.4**2 )
      real,parameter :: roh_min=1.25              ! rOH (1.5**2) in O-H
      real,parameter :: rooc=12.25                ! cutoff distance of rOO (3.5**2 )
      real,parameter :: rohc=6.0025               ! rOH (2.45**2)
      real,parameter :: cosphic=0.866             ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: hb_min=0.0000001 ! condition for the existence of h-bond
      real(kind=rk),parameter :: delta_t=0.0005   ! ps
      real(kind=rk)           :: r12,r13,r23,cosphi,pm,qj,&
                                 tot_hb
      integer :: begin_time,end_time,rat,&
                        i,j,k,jj,nmo,nat,iatom,& 
                        imovie,np,n_count,m1,m2,m3,mt
      real(kind=rk),allocatable,dimension (:)    :: h
      real(kind=rk),allocatable,dimension (:)    :: hb
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3)  :: atom_type  ! we are not interested in the atom type in this calculation,
                                      ! thus I do not use allocatable array for it.
      integer,allocatable,dimension(:)           ::x_1,x_2,x_3,&
                        ndx_1, ndx_2, ndx_3
      real(kind=rk),allocatable,dimension (:)    :: corr_h
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
      read(5,*)nmo     !number of movie steps
      write(6,*)'What is the total number of atoms in the system:'
      read(5,*)nat     !number of atoms per mole.
      write(6,*)'What is the total number of water pairs:'
      read(5,*)np      !number of pairs   

      allocate(x_1(np))          
      allocate(x_2(np))          
      allocate(x_3(np))          
      list_filename=trim(list_filename)
      open(10,file=list_filename)     
      do k=1,np
          read(10,*)x_1(k),x_2(k),x_3(k)
      enddo
      close(10)

      allocate(x(nat,nmo))
      allocate(y(nat,nmo))
      allocate(z(nat,nmo))
      allocate(h(nmo))
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

!======================================
!Simplifying list file
!count the new value of np as n_count
!write new list file with n_count pairs
!======================================
      n_count=0   
      do imovie=1,1
         do k= 1,np
         if(                                              &
             x(x_3(k),imovie)**2-x(x_2(k),imovie)**2+ &
             y(x_3(k),imovie)**2-y(x_2(k),imovie)**2+ &
             z(x_3(k),imovie)**2-z(x_2(k),imovie)**2 <roh_min &
            )then 
             n_count=n_count+1
         endif
         enddo
      enddo
!=======================      
      allocate(ndx_1(np),ndx_2(np),ndx_3(np))
      open(10,file='new_'//trim(list_filename))  
      do imovie=1,1
         do k= 1,np
         if(                                              &
             x(x_3(k),imovie)**2-x(x_2(k),imovie)**2+ &
             y(x_3(k),imovie)**2-y(x_2(k),imovie)**2+ &
             z(x_3(k),imovie)**2-z(x_2(k),imovie)**2 <roh_min &
            )then 
             write(10,*) x_1(k),x_2(k),x_3(k)
         endif
         enddo
      enddo
      close(10)

      deallocate(x_1,x_2,x_3)

!=======================
!Read in new list file      
!=======================
      list_filename='new_'//trim(list_filename)
      open(10,file=list_filename)    
      np=n_count
      do k=1,np
        read(10,*) ndx_1(k),ndx_2(k),ndx_3(k)
      enddo
      close(10)

!=========================================
!Allocate hb and hb_d with new value of np      
!=========================================
      allocate(hb(np))!Average H-bonded population 

!==================================
!Calculate autocorrelation function
!==================================
      ! calculate <h(0)h(t)>/<h> and 
      ! calculate <h^d(0)h^d(t)>/<h^d> 
      ! Notice here <> is not average over
      ! different pairs of water molecules,
      ! but average over the time steps.
      allocate(corr_h(nmo))
      do i=1, nmo
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
        qj=qj/nmo! ave of hb for each pair 
        hb(k)=qj
        tot_hb=tot_hb+hb(k)

        do mt=0,nmo-1! time interval
            if(hb(k)>hb_min) then
                scalar=0.d0
                do j=1, nmo-mt-1
                    scalar=scalar+h(j)*h(j+mt)     ! 1: the first pair of water molecules
                enddo
                scalar=scalar/(nmo-mt-1)         ! C_k(t)
                corr_h(mt+1)=corr_h(mt+1)+scalar  ! sum_C_k(t)
            endif
        enddo
      enddo! k loop 
      tot_hb=tot_hb/np

      do mt=0,nmo-1! time interval
          corr_h(mt+1)=corr_h(mt+1)/(np*tot_hb)  
      enddo
      deallocate(x,y,z)
      deallocate(ndx_1,ndx_2,ndx_3)          
!======================
!Write the correlation
!C_HB(t) and C_HB_d(t)     
!======================
      open(10,file=trim(filename)//'_achbond_h.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,corr_h(i)
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_achbond_h.dat'
      close(10)
!==========================
!Write the correlation
!ln(C_HB(t)) and lnC_HB_d(t)     
!==========================
      open(10,file=trim(filename)//'_achbond_ln_h.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,log(corr_h(i))
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_achbond_ln_h.dat'
      close(10)
!=========
!print <h>      
!=========
      open(10,file=trim(filename)//'_average_h.dat')
        write(10,*) '<h>:',tot_hb
        write(6,*)'written in '//trim(filename)//&
                  '_average_h.dat'
      close(10)
      deallocate (h,corr_h)
!==============================================================
      call system_clock(end_time,rat)
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
