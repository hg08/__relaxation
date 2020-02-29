!2015/10/17
!===================================================
!For ONE pair of molecules，we can obtain
!C_HB(t):
!      
!   C_{HB}(t)= <h(0)h(t)>/<h>
!and      
!Cc_HB(t):
!   Cc_{HB}(t)= <h(0)H(t)>/<h>
!===================================================
!We can implementate these correlation functions by
!calculate the average over N pairs of molecules. 
!===================================================      
! input file: input
     ! name of system
     ! name of trajectory
     ! name of list
     ! nmo
     ! nat
     ! number of pairs of molecules(np) 
!===================================================
      program hbacf_C_and_S_function
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
      real(kind=rk),parameter :: hb_min=0.0000001 ! condition for the existence of h-bond
      real(kind=rk),parameter :: delta_t=0.0005   ! fs
      real(kind=rk)           :: r12,r13,r23,cosphi,pm,qj,&
                                 tot_hb,ave_h,hh
      integer :: begin_time,end_time,rat,&
                 i,j,k,jj,nmo,nat,iatom,& 
                 imovie,np,m1,m2,m3,mt,m 
      real(kind=rk),allocatable,dimension (:)    :: h,hb
      real,allocatable,dimension (:,:)           :: x,y,z
      character(LEN=3)  :: atom_type  ! we are not interested in the atom type in this calculation,
                                      ! thus I do not use allocatable array for it.
      integer,allocatable,dimension(:)           :: ndx_1, ndx_2, ndx_3
      real(kind=rk),allocatable,dimension (:)    :: corr_h,corr_hh
      real(kind=rk)  :: scalar, scalar_hh 
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

      allocate(corr_h(nmo))
      allocate(corr_hh(nmo))
!===================================
!Calculate autocorrelation function
!===================================
! calculate <h(0)h(t)>/<h> and 
! Notice here <> is not average over
! different pairs of water molecules,
! but over all starting time points.
!===================================      
      corr_h(:)=0.0
      corr_hh(:)=0.0
      tot_hb=0.0
      !do k=1,np
      hb(:)=0.0
      !enddo

      do k=1,np
        qj=0
        m1=ndx_1(k)
        m2=ndx_2(k)
        m3=ndx_3(k)
        !Calculate h(j)
        do jj =1, nmo
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
              qj=qj+h(jj)                          
          else                            
              h(jj)=0.0
          endif
        enddo   
        qj=qj/nmo!Ave of hb over all starting points for each pair 
        hb(k)=qj
        tot_hb=tot_hb+hb(k)
!==================================
!Calcualte the correlation function
!==================================
        do mt=0,nmo-1! time interval
            if (hb(k)>hb_min) then
                scalar=0.d0
                scalar_hh=0.d0
                do j=1,nmo-mt-1
                    !--------------------------------------------------------
                    !the definition of H(t): hh(j+mt)=h(j)*h(j+1)*...*h(j+mt)
                    hh=1
                    do m=0,mt
                        if (h(j+m)>hb_min) then
                            !hh=hh*h(j+m) 
                            hh=1
                        else
                            hh=0
                            exit
                        endif    
                    enddo
                    !--------------------------------------------------------  
                    scalar=scalar+h(j)*h(j+mt)  
                    scalar_hh=scalar_hh+h(j)*hh        
                    enddo
                scalar=scalar/(nmo-mt-1)            ! C_k(t)
                scalar_hh=scalar_hh/(nmo-mt-1)      ! Cc_k(t)
                corr_h(mt+1)=corr_h(mt+1)+scalar    ! sum_C_k(t)
                corr_hh(mt+1)=corr_hh(mt+1)+scalar_hh    ! sum_Cc_k(t)
            endif
        enddo
      enddo!End of k-loop 
      tot_hb=tot_hb/np
      ave_h=4*tot_hb

      do mt=0,nmo-1! time interval
          corr_h(mt+1)=corr_h(mt+1)/(np*tot_hb)  
          corr_hh(mt+1)=corr_hh(mt+1)/(np*tot_hb)  
      enddo

      deallocate(x,y,z)
      deallocate(ndx_1,ndx_2,ndx_3)          
!======================
!Write the correlation
!C_HB(t) and S_HB(t)     
!======================
      open(10,file=trim(filename)//'_hbacf_h.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,corr_h(i)
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_hbacf_h.dat'
      close(10)

      open(10,file=trim(filename)//'_hbacf_hh_Cc.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,corr_hh(i)
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_hbacf_hh_Cc.dat'
      close(10)
!===========================
!Write the correlation
!ln(C_HB(t)) and ln(Cc_HB(t))     
!===========================
      open(10,file=trim(filename)//'_hbacf_ln_h.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,log(corr_h(i))
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_hbacf_ln_h.dat'
      close(10)

      open(10,file=trim(filename)//'_hbacf_ln_hh_Cc.dat')
        do i=1,int(nmo*rate)
            write(10,*)(i-1)*delta_t,log(corr_hh(i))
        enddo
        write(6,*)'written in '//trim(filename)//&
                  '_hbacf_ln_hh_Cc.dat'
      close(10)
!===========
! print <h>      
!===========      
      open(10,file=trim(filename)//'_average_h.dat')
        write(10,*) '<h>:',ave_h
        write(6,*)'written in '//trim(filename)//&
                  '_average_h.dat'
      close(10)
      deallocate (h,corr_h,corr_hh)
!=======================
!Print ending time point
!=======================      
      call system_clock(end_time,rat)
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
