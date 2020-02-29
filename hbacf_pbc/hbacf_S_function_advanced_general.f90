!2015/10/21
!huang gang      
!===================================================
!For N pair of molecules，we can obtain
!S_HB(t):
!   S_{HB}(t)= <h(0)H(t)>/<h>
!===================================================      
! With this function one can do calculation with an
! time step, instead of 0.5 ps.
!===================================================      
! input file: input
     ! name of system
     ! name of trajectory
     ! name of list
     ! nmo
     ! nat
     ! number of pairs of molecules(np) 
     ! the new time step, or ns
!===================================================
      program hbacf_S_function
      implicit none
!==========
!parameters
!==========
      character(LEN=200) :: filename ,pos_filename,list_filename       
      integer,parameter :: rk=4              
      real(kind=rk),parameter :: rate=0.80        ! Condition for cutting off autocorrelation functions
      real,parameter :: rooc=12.25                ! Cutoff distance of rOO (5.4**2 )
      real,parameter :: rohc=6.0025                ! rOH (3.5**2)
      real,parameter :: cosphic=0.866             ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: h_min=0.5        ! Condition for h=1.0
      real(kind=rk),parameter :: hb_min=0.5 ! Condition for the existence of h-bond
      real(kind=rk)           :: r12,r13,r23,cosphi,pm,qj,&
                                 tot_hb,hh,tau_hb,delta_t,delta_t0
      integer :: begin_time,end_time,rat,&
                 i,j,k,jj,nmo,nat,iatom,& 
                 imovie,np,m1,m2,m3,mt,m,nqj,tot_nhb,&
                 n_bonded_pairs,ns 
      real(kind=rk),allocatable,dimension (:)    :: h,hb,corr_hh
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3)  :: atom_type!We are not interested in the atom type in this calculation,
                                    !thus I do not use allocatable array for it.
      integer,allocatable,dimension(:)     :: ndx_1, ndx_2, ndx_3,&
           nhb_start
      real(kind=rk) :: scalar_hh 
      logical,allocatable,dimension (:)  :: hb_start
      !logical :: exist!For TESTING
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
      write(6,*)'What is the name of the list file:'
      read(5,*)list_filename     
      write(6,*)'What is the total steps of the trajecotry:'
      read(5,*)nmo     !number of movie steps
      write(6,*)'What is the total number of atoms in the system:'
      read(5,*)nat     !number of atoms per mole.
      write(6,*)'What is the total number of water pairs:'
      read(5,*)np      !number of pairs   
      write(6,*)'What is the time step for calculating CORRELATION:' 
      read(5,*)ns! [ns*0.0005] ps is the new time step for calculation

      allocate(ndx_1(np))          
      allocate(ndx_2(np))          
      allocate(ndx_3(np))          
      list_filename=trim(list_filename)
      open(10,file=list_filename)     
      do k=1,np
          read(10,*)ndx_1(k),ndx_2(k),ndx_3(k)
      enddo
      close(10)

      delta_t=ns*delta_t0 !unit: ps
      nmo=nmo/ns !Length of the correl.function 
      allocate(x(nat,nmo))
      allocate(y(nat,nmo))
      allocate(z(nat,nmo))
      allocate(h(nmo))
      allocate(hb(np))!Average H-bonded population 
      allocate(nhb_start(np)) 
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

         do i=1, (nat+2)*(ns-1) !impilimetation of changing time step  
             read(10,*)
         enddo

      enddo
      close(10)
      write(6,*) 'end of trajectory reading'

      allocate(corr_hh(nmo))
      allocate(hb_start(nmo))
!===================================
!Calculate autocorrelation function
!===================================
! here calculate <h(0)h(t)>/<h>. 
! Notice here <> is not average over
! different pairs of water molecules,
! but over all starting time points.
!===================================      
      tot_hb=0.0
      tot_nhb=0
      corr_hh(:)=0.0
      !do i=1, nmo
      !  corr_hh(i)=0.0
      !enddo
      hb(:)=0.0
      nhb_start(:)=0
      !do k=1,np
      !    hb(k)=0.0
      !    nhb_start(k)=0
      !enddo

      do k=1, np
        qj=0
        nqj=0
        m1=ndx_1(k)
        m2=ndx_2(k)
        m3=ndx_3(k)
        !Calculate h(j)
        do jj =1, nmo
            h(jj)=0
            hb_start(jj)=.False.
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
            if (r13 .lt. rohc .and. r12 .lt. rooc  & 
               .and. cosphi .gt. cosphic) then    
                h(jj)=1.0
                if (jj .eq. 1) then
                    hb_start(jj)=.True. 
                    qj=qj+h(jj)                          
                    nqj=nqj+1 
                else 
                    if(h(jj-1)<h_min) then
                        hb_start(jj)=.True. 
                        qj=qj+h(jj)                          
                        nqj=nqj+1 
                    endif
                endif    
            endif
        enddo! jj-loop   
        hb(k)=qj!Ave of hb over all starting points for each pair is not necessary in calculation of S(t) 
        nhb_start(k)=nqj
        tot_hb=tot_hb+hb(k)
        tot_nhb=tot_nhb+nhb_start(k)
        !==================================
        !Calcualte the correlation function
        !==================================
        if (hb(k)>hb_min) then
            do mt=0,nmo-1!Time interval
                scalar_hh=0.d0
                do j=1, nmo-mt-1
                hh=0.0
                    !====================================================================
                    !The definition/calculation of H(t): hh(j+mt)=h(j)*h(j+1)*...*h(j+mt)
                    if (hb_start(j)) then
                        do m=0,mt
                            if (h(j+m)>hb_min) then
                                hh=1.0  
                            else
                                hh=0.0
                                exit
                            endif    
                        enddo
                    endif
                    !End of definition/calculation of H(t).
                    !====================================================================  
                scalar_hh=scalar_hh+h(j)*hh        
                enddo
            corr_hh(mt+1)=corr_hh(mt+1)+scalar_hh! sum_S_k(t)
            enddo
        endif
      enddo!end of k-loop
      deallocate(hb_start)
      !=========================================
      !Calculate the number of ever bonded pairs
      !=========================================
      n_bonded_pairs=0 
      do k=1,np
          if (hb(k)>hb_min) then
              n_bonded_pairs=n_bonded_pairs+1      
          endif
      enddo
      !====
      !TEST
      !====
      !inquire(file="test.txt", exist=exist)
      !if (exist) then
      !    open(12, file="test.txt", status="old", position="append", action="write")
      !else
      !    open(12, file="test.txt", status="new",position="append", action="write")
      !end if
      !do k=1,np
      !    write(12, *) hb(k)
      !enddo
      !close(12)  
      !=====================
      !Normalization of S(t)
      !=====================
      corr_hh(:)=corr_hh(:)/tot_nhb  
      !do mt=0,nmo-1!time interval
      tau_hb=sum(corr_hh)*delta_t 
      !enddo
      deallocate(x,y,z,ndx_1,ndx_2,ndx_3)          
!=====================
!Write the correlation
!S_HB(t)  
!=====================
      open(10,file=trim(filename)//'_hbacf_hh.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,corr_hh(i)
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_hbacf_hh.dat'
      close(10)
!======================
!Write the correlation
!ln(S_HB(t))     
!======================
      open(10,file=trim(filename)//'_hbacf_ln_hh.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,log(corr_hh(i))
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_hbacf_ln_hh.dat'
      close(10)
!================================
!Print No. bonds starting at t=0.
!================================
      open(10,file=trim(filename)//'_No._bonds_starting_at_t=0.dat')
        write(10,*) 'No._bonds_starting_at_t=0:', tot_nhb
        write(10,*) 'H_bond lifetime(ps):',tau_hb
        write(6,*)'written in '//trim(filename)//&
                  '_No._bonds_starting_at_t=0.dat'
      close(10)
      deallocate (h,corr_hh)
!=====================
!Print the ending time
!=====================      
      call system_clock(end_time,rat)
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
