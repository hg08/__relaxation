!2015/10/17
!===================================================
!hb lifetime distribution
!===================================================      
! input file: input
     ! name of system
     ! name of trajectory
     ! name of list
     ! nmo
     ! nat
     ! number of pairs of molecules(np) 
!===================================================
      program hbacf_lifetime_dist
      implicit none
!==========
!parameters
!==========
      character(LEN=200) :: filename ,pos_filename,list_filename       
      integer,parameter :: rk=4              
      real(kind=rk),parameter :: rate=0.50        ! Condition for cutting off autocorrelation functions
      real,parameter :: rooc=12.25                ! Cutoff distance of rOO (5.4**2 )
      real,parameter :: rohc=6.0025                ! rOH (3.5**2)
      real,parameter :: cosphic=0.866             ! 1.732/2; phiC=pi/6.
      real(kind=rk),parameter :: h_min=0.5        ! Condition for h=1.0
      real(kind=rk),parameter :: hb_min=0.0000001 ! Condition for the existence of h-bond
      real(kind=rk),parameter :: delta_t=0.0005   ! fs
      real(kind=rk)           :: r12,r13,r23,cosphi,pm,qj,&
                                 tot_hb,hh
      integer :: begin_time,end_time,rat,&
                 i,j,k,jj,nmo,nat,iatom,& 
                 imovie,np,m1,m2,m3,mt,m,l,nsubarray 
      real(kind=rk),allocatable,dimension (:)    ::h,hb,corr_hh
      real,allocatable,dimension (:,:)           :: x,y,z,span
      character(LEN=3)  :: atom_type!We are not interested in the atom type in this calculation,
                                    !thus I do not use allocatable array for it.
      integer,allocatable,dimension(:)           :: ndx_1, ndx_2, ndx_3
      real(kind=rk)  :: scalar_hh 
      logical        :: insubarray
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

      allocate(corr_hh(nmo))
!===================================
!Calculate autocorrelation function
!===================================
! here calculate <h(0)h(t)>/<h>. 
! Notice here <> is not average over
! different pairs of water molecules,
! but over all starting time points.
!===================================      
      do i=1, nmo
        corr_hh(i)=0.0
      enddo
      tot_hb=0.0
      do k=1,np
          hb(k)=0.0
      enddo

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
!        do mt=0,nmo-1!Time interval
            if (hb(k)>hb_min) then
!                scalar_hh=0.d0
!                do j=1, nmo-mt-1
                    !====================================================================
                    !The definition/calculation of H(t): hh(j+mt)=h(j)*h(j+1)*...*h(j+mt)
                    nsubarray=0
                    insubarray=.False.
                    !do m=0,mt
                        if (h(j)>h_min) then
                            !hh=hh*h(j+m) ,or
                            if(.not.insubarray) nsubarray = nsubarray+1
                            insubarray=.True.
                            !hh=1.0
                        else
                            insubarray=.False.
                            !hh=0
                        endif    
                    allocate(span(nmo,nsubarray))
        enddo
      enddo!end of k-loop 

                    span=0
                    insubarray=.False.
                    l=1
                    jj=1
                    do j=1,nmo
                        if(h(j)>h_min) then
                                insubarray=.True.
                                span(jj,l)=h(j)
                                jj=jj+1
                        else
                                insubarray=.False.
                        endif
                        if(.not. insubarray) jj=1
                            if (j .gt.1) then
                                if((h(j-1)>h_min) .and. (h(j)<h_min)) then
                                        l=l+1
                                endif
                            endif
                    enddo

                    !Write the output in the required format:
                    open(10,file=trim(filename)//'_achbond.log',&
                        status='old')
                    do i=1, nsubarray
                        write(20,'(a1,i1,a2)',advance="no")'span',i,'=('
                        do j=1,nmo
                            if (j .lt.nmo) then
                                write(20,'(i1,a1)',advance="no"),&
                                    span(j,i),','
                            else
                                write(20,'(i1)',advance="no") span(j,i)
                            endif    
                        enddo
                        write(20,'(a1)')')'
                     enddo        
                     close(10)
            endif
        enddo
      enddo!end of k-loop 
      tot_hb=tot_hb/np
      deallocate(x,y,z,ndx_1,ndx_2,ndx_3,span)          
!======================
! print <h>      
!======================      
      open(10,file=trim(filename)//'_average_h.dat')
        write(10,*) '<h>:',tot_hb
        write(6,*)'written in '//trim(filename)//&
                  '_average_h.dat'
      close(10)
      deallocate (h,corr_hh)
!=====================
!Print the ending time
!=====================      
      call system_clock(end_time,rat)
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END

