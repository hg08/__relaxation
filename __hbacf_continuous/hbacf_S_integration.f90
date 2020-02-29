!2015/1/1
!huang gang      
!===================================================
!For N pair of molecules，we can obtain integration 
!of S_HB(t):
!   S_{HB}(t)= <h(0)H(t)>/<h>
!===================================================      
! input file: input
     ! name of system
     ! name of correlation file 
     ! nmo
!===================================================
      program hbacf_S_function
      implicit none
!==========
!parameters
!==========
      character(LEN=200) :: filename ,correl_filename       
      integer,parameter :: rk=4              
      real(kind=rk),parameter :: delta_t=0.0005   ! ps
      real(kind=rk)           :: imo,tau_hb,tau_hb1
      integer :: begin_time,end_time,rat,&
                 nmo,imovie
      real(kind=rk),allocatable,dimension (:)    :: corr_hh
      call system_clock(begin_time,rat)
!==================
!read data in input
!==================
      write(6,*)'What is the name of the system:'
      read(5,*)filename
      write(6,*)'What is the name of the correlation file S(t):'
      read(5,*)correl_filename     
      write(6,*)'What is the total steps of the trajecotry:'
      read(5,*)nmo     !number of movie steps

      allocate(corr_hh(nmo))
!==================
!read in correl file 
!==================
      open(10,file=trim(correl_filename))     
      do imovie=1,nmo
             read (10,*)imo,corr_hh(nmo)
      enddo
      close(10)
      write(6,*) 'end of correlation reading'
      tau_hb=sum(corr_hh)*delta_t
      tau_hb1=sum(corr_hh)
!===================================
!Write result 
!===================================
      open(10,file=trim(filename)//'_hb_lifetime.dat')
        write(6,*) 'Life time of Hbonds', tau_hb1
        write(10,*) 'Life time of Hbonds', tau_hb
        write(6,*)'written in '//trim(filename)//&
                  '_hb_lifetime.dat'
      close(10)
      deallocate (corr_hh)
!=====================
!Print the ending time
!=====================      
      call system_clock(end_time,rat)
      write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 
      END
