!This program has to extract the position of central atoms according to their lifetime
!This program has to be used after orientation.f95

!Information about the main tables
!=========
!tab_nb_dt
!=========
  !tab_nb_dt(1,:) contains the minimal lifetime taken into account
  !tab_nb_dt(2,:) contains the maximal lifetime taken into account
  !tab_nb_dt(3,:) contains the position where starts a new lifetime in tabz
  !tab_nb_dt(4,:) is a counter to know what is the next index where will be writen the data
  !tab_nb_dt(5,:) is a counter to know what is the number of molecules with the good lifetime
  !tab_nb_dt(6,:) contains the position where starts a new lifetime in tabz when the values are over limit2
!========
!tab_life
!========
  !tab_life(1,:) contains the line number (initially the number of the atom)
  !tab_life(2:3,:) contains the number of the ligands
  !tab_life(4:5,:) contains the starting and the final steps where we have a molecule
  !tab_life(6,:) contains the lifetime of the molecule
!====
!tabz
!====
  !tabz is a 1D table which will contain all the z positions which will be sorted according to their lifetime. tab_nb_dt(3,:) allows us to know where starts and stops a lifetime. tab_nb_dt(4,i) allows us to know where is the next value for a specific lifetime.
  !Two others solutions were considered instead of this 1D tabz:
  !1)Read several tme the z-file and take only the information about a specific lifetime (maybe too many "rewind" to do)
  !2)Build tabz as a 2D table, but this table can be too big to be allocated: number_of_lifetime_to_consider*maximal_number_of_data_stored_for_1_lifetime_window

program main

  implicit none

  integer i,j,k,l,m,istep,start_step,end_step,dt_sel,nb_dt,max_life,nb_center,end_file,end_line
  integer, dimension(:), allocatable:: tab_center
  integer, dimension(:,:), allocatable ::tab_life,tab_nb_dt
  double precision cell,limit1,limit2,max_mol,tmp_real,dt,inc1_dt,inc2_dt
  double precision, dimension(:), allocatable :: tabz
  character char
  character(len=200) filename,tmp

  !=================================
  !Reading the lifetime file (1/2)
  !Line number: to allocate tab_life
  !=================================
  write(*,*)"Filename of the lifetime file"
  read(*,*)filename
  open(10,file=filename)

  read(10,*)
  read(10,*)
  read(10,*)start_step,end_step
  
  max_life=-1!The the last line (blank line) of 'filename' has to be removed
  end_file=0
  do while(end_file==0)
     read(10,*,iostat=end_file)
     max_life=max_life+1
  enddo
  allocate(tab_life(6,max_life))

  !===============================
  !Reading the lifetime file (2/2)
  !Filling tab_life
  !===============================
  rewind(10)
  read(10,*)  
  read(10,*)
  read(10,*)
  do i=1,max_life
     read(10,*)tab_life(1:5,i)
  enddo
  tab_life(6,:)=tab_life(5,:)-tab_life(4,:)+1
  close(10)


  !====================================================
  !Reading of the basic information of the z-file (1/2)
  !Allocation of tab_center
  !====================================================
  write(*,*)"Filename of the z-projection of the central atom"
  read(*,*)filename
  open(10,file=filename)
  !I just read the first line (comment) and the beginning of the second line which contains a comment
  read(10,*)
  char=""  
  do while (char/=":")
     read(10,'(a1)',advance='no')char
  enddo
  !Number of central atoms
  end_line=0
  nb_center=-1!nb_center will be counted one extra time, so we have to substract 1
  do while(end_line==0)
     read(10,'(i6)',iostat=end_line,advance='no')i
     nb_center=nb_center+1
  enddo
  allocate(tab_center(nb_center))

  !====================================================
  !Reading of the basic information of the z-file (2/2)
  !Filling of tab_center
  !Rewriting of tab_life(1,i)
  !====================================================
  rewind(10)
  !I just read the first line (comment) and the beginning of the second line which contains a comment
  read(10,*)
  char=""  
  do while (char/=":")
     read(10,'(a1)',advance='no')char
  enddo
  !Now we know the size of the table, we record the value of the central atoms
  do i=1,nb_center-1
     read(10,'(i6)',advance='no')tab_center(i)
  enddo
  read(10,'(i6)')tab_center(nb_center)
  !We can convert the tab_life(1,:) from the number of the atom to its column in the z-file
  !tab_life is sorted 1st about the final time, 2nd about the number of the central atom
  do i=1,max_life
     do j=1,nb_center
        if(tab_center(j)==tab_life(1,i))tab_life(1,i)=j
     enddo
  enddo

  !======================
  !Additional information
  !======================
  write(*,*)"Cell parameter along the normal axis (Ang)"
  read(*,*)cell
  write(*,*)"Center of the vacuum/solid slab (Ang)"
  read(*,*)limit1
  limit2=limit1+cell/2!It is the center of the slab which is studied
  write(*,*)"Filename of the output"
  read(*,*)filename

  !===============================================
  !Division of the molecules lifetime into Delta t
  !===============================================
  write(*,*)"For the lifetime windows, do you want:"
  write(*,*)"   -a constant timestep (1)"
  write(*,*)"   -a constant ratio 'timestep' over 'average value' (2)"
  !To have a constant ratio, we need to have Mi+1=M1^(i+1)/M0^i
  !Where M0 is the starting minimum for the lifetime (fixed at 1),
  !M1 is the minimum lifetime of the second lifetime window ('the maximum lifetime of the first window'+1)
  !Mi+1 is the minimum lifetime of the i+1 th lifetime window ('the maximum lifetime of the i th window'+1)
  read(*,'(i10)')dt_sel

  if(dt_sel==1)then!If constant time step
     write(*,*)"Step division (Should be an integer (a rounding is done). Default: 1)"
     tmp=""
     read(*,'(a)')tmp
     if(tmp=="")then
        dt=1
     else
        read(tmp,*)dt
     endif
     !Allocation of tab_nb_dt
     nb_dt=ceiling(maxval(tab_life(6,:))/anint(dt))
     do i=1,nb_dt
        j=(i-1)*nint(dt)+1
        k=i*nint(dt)+1
        if(.not.any((tab_life(6,:)>=j).and.(tab_life(6,:)<k)))then
           nb_dt=nb_dt-1
        endif
     enddo
     
  else!If constant ratio
     write(*,*)"The first lifetime division starts at 1."
     write(*,*)"What is the beggining of the second lifetime windows (Real accepted. Default: 1.01):"
     tmp=""
     read(*,'(a)')tmp
     if(tmp=="")then
        dt=1.01
     else
        read(tmp,*)dt
     endif

     !Allocation of tab_nb_dt
     nb_dt=floor(log(dble(maxval(tab_life(6,:))))/log(dt))+1
     do i=1,nb_dt
        !WARNING!
        !WARNING!What is the good increment
        inc1_dt=dt**(i-1)
        inc2_dt=dt**i
        !For small i and dt, sometimes there is no integer between j and k. 
        !Therefore, it is useless to test if there is a lifetime between them.
        if(ceiling(inc1_dt)<ceiling(inc2_dt))then
           if(.not.any((tab_life(6,:)>=inc1_dt).and.(tab_life(6,:)<inc2_dt)))then
              nb_dt=nb_dt-1
           endif
        endif
     enddo
  endif

  allocate(tab_nb_dt(6,nb_dt))
  
  if(dt_sel==1)then!If constant time step
     !Writing the boundaries of tab_nb_dt
     tab_nb_dt(3,:)=0
     j=1
     do i=1,ceiling(maxval(tab_life(6,:))/anint(dt))
        l=(i-1)*nint(dt)+1!minimal boundary
        m=i*nint(dt)+1!maximal boundary
        do k=1,max_life
           if((tab_life(6,k)>=l).and.(tab_life(6,k)<m))then
              tab_nb_dt(3,j)=tab_nb_dt(3,j)+tab_life(6,k)          
           endif
        enddo
        if(tab_nb_dt(3,j)>0)then
           tab_nb_dt(1,j)=l
           tab_nb_dt(2,j)=m-1
           j=j+1
        endif
     enddo

  else!If constant ratio

     !Writing the boundaries of tab_nb_dt
     tab_nb_dt(3,:)=0
     j=1
     do i=1,floor(log(dble(maxval(tab_life(6,:))))/log(dt))+1
        inc1_dt=dt**(i-1)!minimal boundary
        inc2_dt=dt**i!maximal boundary
        if(ceiling(inc1_dt)<ceiling(inc2_dt))then
           do k=1,max_life
              if((tab_life(6,k)>=inc1_dt).and.(tab_life(6,k)<inc2_dt))then
                 tab_nb_dt(3,j)=tab_nb_dt(3,j)+tab_life(6,k)          
              endif
           enddo

           if(tab_nb_dt(3,j)>0)then
              tab_nb_dt(1,j)=ceiling(inc1_dt)
              tab_nb_dt(2,j)=ceiling(inc2_dt)-1
              j=j+1
           endif
        endif
     enddo

  endif

  !Rewriting the tab_nb_dt(3,i)
  j=tab_nb_dt(3,1)
  tab_nb_dt(3,1)=1
  do i=2,nb_dt
     k=tab_nb_dt(3,i)
     tab_nb_dt(3,i)=tab_nb_dt(3,i-1)+j
     j=k
  enddo

  !allocation of tabz
  allocate(tabz(tab_nb_dt(3,nb_dt)+j-1))


  !========================================================================================
  !Reading of the main part of the z-file
  !The z-file (so tab_life) is first sorted according to the destruction step (5th collumn)
  !And then according to the creation step (4th collumn)...
  !========================================================================================

  !Getting of the full information (I do not take into account the molecules already formed when the file starts, neither the molecules not yet destroyed after the end of the file)
  tab_nb_dt(4,:)=tab_nb_dt(3,:)!Initialization of the counter
  tab_nb_dt(5,:)=0

  read(10,'(i10)',advance='no',iostat=end_file)istep
  do while(end_file==0)
 
     do i=1,nb_center
        read(10,'(a10)',advance='no')tmp!Sometimes the value is an integer, sometimes it is x, so I read a character chain
        if(tmp(10:10)/="x")then
           do j=1,max_life
              if(tab_life(1,j)==i)then!tab_life is ordered according to its ending step, so the first j with the good central atom is the one which has to be studied.
                 if((tab_life(4,j)/=start_step).or.(tab_life(5,j)/=end_step))then!The molecules created/destroyed before/after the beggining/end of the z-file are not taken into account
                    do k=1,nb_dt
                       if((tab_nb_dt(1,k)<=tab_life(6,j)).and.(tab_nb_dt(2,k)>=tab_life(6,j)))then!We search where in tabz we have to write the position (tab_nb_dt(4,k))
                          read(tmp,*)tabz(tab_nb_dt(4,k))!Conversion from a chain to a real
                          if(tabz(tab_nb_dt(4,k))<limit1)tabz(tab_nb_dt(4,k))=tabz(tab_nb_dt(4,k))+cell!Centering the data
                          tab_nb_dt(4,k)=tab_nb_dt(4,k)+1
                          exit
                       endif
                    enddo
                 endif
                 
                 !At the end of the life of the molecule
                 !The occurence counter is incremented (if start/end of the molecule are "in" the z-file)
                 !and the lines of tab_life are removed if istep is too big
                 if(tab_life(5,j)==istep)then
                    !Counter
                    if((tab_life(4,j)/=start_step).or.(tab_life(5,j)/=end_step))then
                       do k=1,nb_dt
                          if((tab_life(6,j)>=tab_nb_dt(1,k)).and.(tab_life(6,j)<=tab_nb_dt(2,k)))then
                             tab_nb_dt(5,k)=tab_nb_dt(5,k)+1
                          endif
                       enddo
                    endif

                    !Line removing
                    do k=j,max_life-1
                       tab_life(:,k)=tab_life(:,k+1)
                    enddo
                    max_life=max_life-1
                 endif

                 exit

              endif
           enddo!loop about tab_life
        endif
     enddo

     if(mod(istep,1000)==0)write(*,*)"Step",istep,"treated."
     read(10,'(i10)')
     read(10,'(i10)',advance='no',iostat=end_file)istep

  end do!while loop reading the file
  close(10)

  
  !============
  !Sorting tabz
  !============
  write(*,*)
  write(*,*)"Data recorded. Sorting the data."

  do i=1,nb_dt
     do j=tab_nb_dt(3,i)+1,tab_nb_dt(4,i)-1
        tmp_real=tabz(j)
        do k=j-1,tab_nb_dt(3,i),-1
           if(tmp_real<tabz(k))then
              tabz(k+1)=tabz(k)
           else
              exit
           endif
        enddo
        tabz(k+1)=tmp_real
     enddo
     if(mod(i,50)==0)write(*,*)i," over",nb_dt," lifetime windows treated."
  enddo

  !======================================
  !Spliting data according to the symetry
  !======================================
  write(*,*)
  write(*,*)"Data sorted. Spliting the data."
  do i=1,nb_dt
     tab_nb_dt(6,i)=tab_nb_dt(4,i)
     do j=tab_nb_dt(3,i),tab_nb_dt(4,i)-1
        if(tabz(j)>limit2)then
           tab_nb_dt(6,i)=j
           exit
        endif
     enddo
  enddo


  !================
  !Writing the data
  !================
  open(10,file=filename)
  write(10,*)"#Projection on the axis of the molecules with this lifetime"
  write(10,*)"#Two lines have to be read together to have the full slab description (1 line represent one half of the slab)"
  write(10,*)"# #",limit2," is the limit between the first and the second half"!If there is only one sharp, awk does not read limit2
  write(10,'(a)',advance='no')" #Lifetime init (step) / Lifetime end (step) / "
  write(10,'(a)')"Number of molecules concerned (max=1) / Occurence on the side i (max=1) / z1 z2 z3 ...:"

  max_mol=dble(maxval(tab_nb_dt(5,:)))

  do i=1,nb_dt
     if(tab_nb_dt(5,i)>0)then
        !Two lines are writen: one for both side of limit2
        write(10,'(2i10,2f10.6)',advance='no')tab_nb_dt(1,i),tab_nb_dt(2,i),&
             tab_nb_dt(5,i)/max_mol,dble(tab_nb_dt(6,i)-tab_nb_dt(3,i))/(tab_nb_dt(4,i)-tab_nb_dt(3,i))
        do j=tab_nb_dt(3,i),tab_nb_dt(6,i)-1
           write(10,'(f10.5)',advance='no')tabz(j)
        enddo
        write(10,*)
        write(10,'(2i10,2f10.6)',advance='no')tab_nb_dt(1,i),tab_nb_dt(2,i),&
             tab_nb_dt(5,i)/max_mol,dble(tab_nb_dt(4,i)-tab_nb_dt(6,i))/(tab_nb_dt(4,i)-tab_nb_dt(3,i))
        do j=tab_nb_dt(6,i),tab_nb_dt(4,i)-1
           write(10,'(f10.5)',advance='no')tabz(j)
        enddo
        write(10,*)

     endif
  enddo

  close(10)
  deallocate(tab_life,tab_center,tab_nb_dt,tabz)

end program main




