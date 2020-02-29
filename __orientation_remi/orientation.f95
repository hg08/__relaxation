module derived_type
  implicit none
  type info_lig
     sequence!To force to store the information continuously
     character(len=3) :: symbol!Symbol of the ligand
     double precision :: rmin!minimal radius 
     double precision :: rmax!maximal radius
     integer :: number!Number of identical ligand (same symbol, same rmin, same rmax)
     integer :: start!Number used for the output filename if there are 2 ligands with the same symbol but with different bond length
     integer :: max!Number used if there are two ligands with the same symbol with different bond length
  end type info_lig

  type info_result
     sequence!To force to store the information continuously
     integer :: number!Number of the atoms considered as a ligand
     double precision :: x
     double precision :: y
     double precision :: z
     double precision :: r
  end type info_result

end module derived_type



program main

  use derived_type
  implicit none

  interface    
     subroutine atom_selection(nb_at,tab_at,tab_symb,symbol)
       integer , intent(in) ::nb_at
       logical , intent(inout) :: tab_at
       character(len=1) , intent(in) :: symbol
       character(len=3) , intent(in) :: tab_symb
     end subroutine atom_selection

     subroutine record(nb_at,nb_lig,max_lig_same,tab_center,tab_lig,&
          a,b,c,tabx,taby,tabz,tab_result,tab_info_lig)
       use derived_type
       integer, intent(in) :: nb_at,nb_lig,max_lig_same
       logical, dimension(nb_at), intent(in) :: tab_center
       logical, dimension(nb_at,nb_lig), intent(in) :: tab_lig
       double precision, intent(in) :: a,b,c
       double precision, dimension(nb_at), intent(in) :: tabx,taby,tabz
       type (info_lig), dimension(nb_lig),intent(in) :: tab_info_lig
       type (info_result), dimension(max_lig_same,nb_lig,nb_at), intent(inout):: tab_result
     end subroutine record

     recursive subroutine permutation(position,position_max,nb_tot_lig,valid,tab_r_ref,tab_permut,tab_order)
       use derived_type
       integer,intent(in) ::  position,position_max,nb_tot_lig
       integer, dimension(nb_tot_lig),intent(inout) :: tab_order
       logical,intent(inout) ::  valid
       double precision, dimension(2,nb_tot_lig),intent(in) :: tab_r_ref
       type (info_result), dimension(nb_tot_lig),intent(inout) :: tab_permut
     end subroutine permutation

     subroutine sharing(nb_at,nb_lig,max_lig_same,tab_result,tab_info_lig)
       use derived_type
       integer, intent(in) :: nb_at,nb_lig,max_lig_same
       type(info_lig),dimension(nb_lig),intent(in) :: tab_info_lig
       type (info_result), dimension(max_lig_same,nb_lig,nb_at), intent(inout):: tab_result
     end subroutine sharing
     
     subroutine lifetime(nb_at,nb_lig,nb_tot_lig,max_lig_same,tab_start_life,istep,tab_life,tab_info_lig,tab_result)
       use derived_type
       integer,intent(in) :: nb_at,nb_lig,nb_tot_lig,max_lig_same,istep
       integer, dimension(2,nb_at),intent(inout) :: tab_start_life
       integer,dimension(max_lig_same,nb_lig,nb_at),intent(inout) :: tab_life
       type(info_lig),dimension(nb_lig),intent(in) :: tab_info_lig
       type (info_result),dimension(max_lig_same,nb_lig,nb_at),intent(in) :: tab_result
     end subroutine lifetime

  end interface

  integer i,j,k,l,m,n,iligand,position,nb_at,nb_lig,nb_tot_lig,istep,step_i,step_f,end_file,ifile,max_file,statistics,max_lig_same
  integer(8), dimension(2) :: pos_file
  integer, dimension(:,:), allocatable :: tab_start_life
  integer, dimension(:), allocatable :: tab_order
  integer, dimension(:,:,:), allocatable :: tab_life
  logical valid,share
  logical, dimension(:), allocatable :: tab_center
  logical, dimension(:,:), allocatable :: tab_lig
  double precision a,b,c,pi
  double precision, dimension(:),allocatable :: tabx,taby,tabz
  double precision, dimension(:,:), allocatable :: tab_r_ref
  character char
  character(len=3)symb_center
  character(len=3),dimension(:),allocatable::tab_symb
  character(len=200)filein,fileout,tmp
  type(info_lig),dimension(:),allocatable :: tab_info_lig
  type(info_result),dimension(:),allocatable:: tab_permut
  type(info_result),dimension(:,:,:),allocatable:: tab_result

  pi=atan2(1.,1.)*4
  statistics=0!Number of central atoms which have not the good ligands


  !==============
  !Basic qestions
  !==============
  write(*,*)"Name of the input file"
  read(*,*)filein
  open(10,file=filein,status='old')
  read(10,*)nb_at
  write(*,*)"There are",nb_at,"atoms."
  allocate(tab_center(nb_at),tab_symb(nb_at),tabx(nb_at),taby(nb_at),tabz(nb_at))

  read(10,*)char,char,istep
  !  write(*,*)"The starting step is", istep
  !  write(*,*)

  write(*,*)
  write(*,*)"Symbol of the central atom:"
  read(*,*)symb_center
  write(*,*)"How many ligands?"
  read(*,*)nb_tot_lig
  if(nb_tot_lig>0)then
     allocate(tab_info_lig(nb_tot_lig),tab_permut(nb_tot_lig),tab_r_ref(2,nb_tot_lig),tab_order(nb_tot_lig))
  endif
  do iligand=1,nb_tot_lig
     write(*,*)"Symbol of the ligand number",iligand,":"
     read(*,*)tab_info_lig(iligand)%symbol
     tab_info_lig(iligand)%symbol=adjustl(tab_info_lig(iligand)%symbol)
     tab_info_lig(iligand)%number=1
     write(*,*)"Minimal distance defining this ",trim(symb_center),"-",trim(tab_info_lig(iligand)%symbol),":"
     read(*,*)tab_info_lig(iligand)%rmin
     write(*,*)"Maximal distance defining this ",trim(symb_center),"-",trim(tab_info_lig(iligand)%symbol),":"
     read(*,*)tab_info_lig(iligand)%rmax
  enddo

  !We arrange tab_info_lig to have the name of the ligands (without redundance) and their amount.
  !nb_tot_lig: Number of ligand surrounding the central molecule
  !nb_lig: Number of kind of ligand
  i=1
  nb_lig=nb_tot_lig
  do while (i<nb_lig)
     do j=i+1,nb_lig
        do while(  (tab_info_lig(i)%symbol==tab_info_lig(j)%symbol).and.&!If the symbol is the same
                   (tab_info_lig(i)%rmin  ==tab_info_lig(j)%rmin  ).and.&!... if rmin too and
                   (tab_info_lig(i)%rmax  ==tab_info_lig(j)%rmax  ))  !... if rmax too, therefore the two ligands are the same and we just increase the increment tab_info_lig(i)%number
        
           tab_info_lig(i)%number=tab_info_lig(i)%number+1!... we consider that the ligands are the same   

           do k=j+1,nb_lig
              tab_info_lig(k-1)=tab_info_lig(k)
           enddo
           tab_info_lig(nb_lig)%symbol=""
           tab_info_lig(nb_lig)%number=0
           nb_lig=nb_lig-1
        enddo
     end do
     i=i+1
  enddo

  !We complete tab_info_lig%max which is equals to tab_info_lig%number if:
  !   1) 2 different ligands have 2 different symbol
  !   2) 2 different ligands with the same symbol have rmax(1)<rmin(2)
  !If there are at least 2 ligands with the same symbol but with bonds which are almost the same:
  !   1) rmin(1)<rmin(2)<rmax(1)<rmax(2)
  !   2) rmin(1)<rmin(2)<rmax(2)<rmax(1)
  !therefore tab_info_lig%max is equals to the the sum of the tab_info_lig%number
  !The interest of this value is to speed up the code if there is no intersection between the 2 radius windows
  tab_info_lig(:)%max=tab_info_lig(:)%number
  do i=1,nb_lig
     do j=i+1,nb_lig
        if(tab_info_lig(i)%symbol==tab_info_lig(j)%symbol)then!If the symbol are the same
           if(  ((tab_info_lig(i)%rmin<=tab_info_lig(j)%rmin).and.&!If there is an intersection between two radius
                 (tab_info_lig(i)%rmax>=tab_info_lig(j)%rmin))&
                .or.&
                ((tab_info_lig(j)%rmin<=tab_info_lig(i)%rmin).and.&
                 (tab_info_lig(j)%rmax>=tab_info_lig(i)%rmin)))then
              tab_info_lig(i)%max=tab_info_lig(i)%max+tab_info_lig(j)%number
              tab_info_lig(j)%max=tab_info_lig(j)%max+tab_info_lig(i)%number
           endif!Radius
        endif!Symbol
     enddo
  enddo
  if(nb_lig>0)then
     max_lig_same=maxval(tab_info_lig(:)%max)
     allocate(tab_result(max_lig_same,nb_lig,nb_at),tab_life(max_lig_same,nb_lig,nb_at),tab_start_life(2,nb_at))
  endif

  !Now we have to take care about the filename of the output (e.g. numerotation of 2 O-H bonds of ~1Ang and 2 O-H bonds of ~1.5 Ang: the two first H will be noted 1, 2 and the two last will be noted 3, 4)
  do i=1,nb_lig
     tab_info_lig(i)%start=1!It is the role of the "start" parameter to avoid this problem
     do j=1,i-1
        if(tab_info_lig(i)%symbol==tab_info_lig(j)%symbol)then!If the two ligands have the same symbol (but different bonds because i/=j)         
           tab_info_lig(i)%start=tab_info_lig(i)%start+tab_info_lig(j)%number!... therefore we take them into account for the output name
        endif
     enddo
  enddo


  if(nb_lig>0)then
     allocate(tab_lig(nb_at,nb_lig))!Tab_lig will contain the info about the atoms which have to be considered as ligand (with the value nb_lig)
     
     write(*,*)"Do you accept to have one ligand shared by two central atom? (y/n - default n)"
     read(*,'(a)')tmp
     tmp=adjustl(tmp)
     if((trim(tmp)=="").or.(trim(tmp)=="n"))then
        share=.false.
     else
        share=.true.
     endif

     write(*,*)
     write(*,*)"!!!    Only MXn molecules with the 'good' bonds will be taken into account    !!!"   
     write(*,*)"!!!    MXn-1, MXn+1 or MXn with short bonds will not be taken into account    !!!"
     if(.not.share)then
        write(*,*)"!!!If a X is shared by 2 M, none of these molecules will be taken into account!!!"
     endif
     write(*,*)

     do i=1,nb_lig
        write(*,'(a a a f6.2 a f6.2 a)',advance='no')"Label of the ligand ", trim(tab_info_lig(i)%symbol),&
             " with rmin=", tab_info_lig(i)%rmin, "Ang and rmax=", tab_info_lig(i)%rmax,"Ang: "
        do j=tab_info_lig(i)%start,tab_info_lig(i)%start+tab_info_lig(i)%number-1
           write(char,'(i1)')j!Convert an integer into a character (j is lower than 10)
           write(*,'(a a a)',advance='no')trim(tab_info_lig(i)%symbol),char," "
        enddo
        write(*,*)
     enddo
  endif
  write(*,*)


  write(*,*)"Write the root name of the output files"
  read(*,*)fileout
  open(11,file=trim(fileout)//"_x"//trim(symb_center)//".dat")
  open(12,file=trim(fileout)//"_y"//trim(symb_center)//".dat")
  open(13,file=trim(fileout)//"_z"//trim(symb_center)//".dat")
  write(*,*)"The output files will be:"
  write(*,*)"    -",trim(fileout),"_x",trim(symb_center),".dat"
  write(*,*)"    -",trim(fileout),"_y",trim(symb_center),".dat"
  write(*,*)"    -",trim(fileout),"_z",trim(symb_center),".dat"
  max_file=13
  do i=1,nb_lig
     do j=tab_info_lig(i)%start,tab_info_lig(i)%start+tab_info_lig(i)%number-1
        write(char,'(i1)')j!Convert an integer into a character (j is lower than 10)
        open(max_file+1,file=trim(fileout)//"_vx"//trim(tab_info_lig(i)%symbol)//char//".dat")
        open(max_file+2,file=trim(fileout)//"_vy"//trim(tab_info_lig(i)%symbol)//char//".dat")
        open(max_file+3,file=trim(fileout)//"_vz"//trim(tab_info_lig(i)%symbol)//char//".dat")
        write(*,*)"    -",trim(fileout),"_vx",trim(tab_info_lig(i)%symbol),char,".dat"
        write(*,*)"    -",trim(fileout),"_vy",trim(tab_info_lig(i)%symbol),char,".dat"
        write(*,*)"    -",trim(fileout),"_vz",trim(tab_info_lig(i)%symbol),char,".dat"
        max_file=max_file+3
     enddo
  enddo
  !File about lifetime
  if(nb_lig>0)then
     open(max_file+1,file=trim(fileout)//"_life-over.dat")
     open(max_file+2,file=trim(fileout)//"_life-under.dat")
     max_file=max_file+2
     write(*,*)"    -",trim(fileout),"_life-over.dat"
     write(*,*)"    -",trim(fileout),"_life-under.dat"
  endif


  write(*,*)
  write(*,*)"What are the crystal parameters? (Ang)"
  read(*,*)a,b,c

  !Reading the first step to be able to associate a symbol to a line number
  do i=1,nb_at
     read(10,*)tab_symb(i),tabx(i),taby(i),tabz(i)
  enddo

  !Central atom selection
  write(*,*)
  write(*,*)trim(symb_center)," atoms:"
  tab_center=.false.
  call atom_selection(nb_at,tab_center(1),tab_symb(1),symb_center)

  !Ligand selection
  tab_lig=.false.
  do i=1,nb_lig
     write(*,*)
     write(*,*)trim(tab_info_lig(i)%symbol)," atoms:"
     call atom_selection(nb_at,tab_lig(1,i),tab_symb(1),tab_info_lig(i)%symbol)
  enddo

  !First lines of the output files
  write(11,*)"#Projection along x of the ",trim(symb_center)," atoms"
  write(12,*)"#Projection along y of the ",trim(symb_center)," atoms"
  write(13,*)"#Projection along z of the ",trim(symb_center)," atoms"
  do ifile=11,13
     write(ifile,'(a)',advance='no')"#Step / Associated central atom:"
     do k=1,nb_at
        if(tab_center(k))write(ifile,'(i6)',advance='no')k
     enddo
     write(ifile,*)
  enddo
  do i=1,nb_lig
     do j=tab_info_lig(i)%start,tab_info_lig(i)%start+tab_info_lig(i)%number-1
        do k=0,2
           select case(k) 
           case(0)
              write(ifile,*)"#Projection along x of ",trim(symb_center),"-",trim(tab_info_lig(i)%symbol),j
           case(1)
              write(ifile,*)"#Projection along y of ",trim(symb_center),"-",trim(tab_info_lig(i)%symbol),j
           case(2)
              write(ifile,*)"#Projection along z of ",trim(symb_center),"-",trim(tab_info_lig(i)%symbol),j
           end select
           write(ifile,'(a)',advance='no')"#Step / Associated central atom "
           do l=1,nb_at
              if(tab_center(l))write(ifile,'(i6)',advance='no')l
           enddo
           write(ifile,*)
           ifile=ifile+1
        enddo
     enddo
  enddo

  !Step boundaries
  write(*,*)
  write(*,*)"Starting step (default:",istep,")"
  read(*,'(i30)')step_i
  if(step_i<istep)step_i=istep
  write(*,*)"Final step (default: end of file)"
  read(*,'(a)')tmp!I want to make a difference between step_f=0 (the final step is 0) and step_f="" (the final step is the last step of the input file)
  if(tmp=="")then
     step_f=-1!The input file will be read untill the end if we do not precise step_f
  else
     read(tmp,'(i10)')step_f
  endif
  if(step_f<step_i)step_f=-1!The input file will be read untill the end if we have step_f < step_i

  !Writing the first lines of the lifetime files
  if(nb_lig>0)then
     write(14+nb_tot_lig*3,*)"#Overestimated lifetime of the molecule"
     write(14+nb_tot_lig*3,*)"#Starting / Final step of the dynamic"!These data can be different from the minimal/maximal steps when molecules are created/destroyed 
     write(14+nb_tot_lig*3,'(i10)',advance='no')step_i
     pos_file(1)=ftell(14+nb_tot_lig*3)!We record this position to put at the end the value of the final step of the dynamic
     write(14+nb_tot_lig*3,*)step_i!step_i is writen just to save some space in the file (the final step will be writen on this space)
     write(15+nb_tot_lig*3,*)"#Underestimated lifetime of the molecule"
     write(15+nb_tot_lig*3,*)"#Starting / Final step of the dynamic" 
     write(15+nb_tot_lig*3,'(i10)',advance='no')step_i
     pos_file(2)=ftell(15+nb_tot_lig*3)
     write(15+nb_tot_lig*3,*)step_i
   endif



  !=======================
  !Passing the first steps
  !=======================
  !This code is not optimized according to the number of writen lines but according to the execution time
  !"read(10,*)" is arround 3 times faster than "read(10,*)tab_symb(i),tabx(i),taby(i),tabz(i)"
  if(istep<step_i)then
     read(10,*,iostat=end_file)!Have we reached the end of the file before step_i?
     if(end_file/=0)then
        write(*,*)"Starting step too high! No treated data!"
        stop
     endif
     read(10,*)char,char,istep

     !I do not use a while loop because...
     !Either I do one extra loop (because I have not "else")
     !Or I will not take into account the very last step (because of if(end_file/=0))
     do  
        if(istep<step_i)then
           do i=1,nb_at
              read(10,*)
           enddo
           read(10,*,iostat=end_file)!Have we reached the end of the file before step_i?
           if(end_file/=0)then
              write(*,*)"Starting step too high! No treated data!"
              stop
           endif
           read(10,*)char,char,istep
        else
           do i=1,nb_at
              read(10,*)tab_symb(i),tabx(i),taby(i),tabz(i)!Slower than read(10,*)           
           enddo
           exit
        endif
     enddo
  endif

  !==========================
  !Recording of the good data
  !==========================
  do while((istep<=step_f).or.(step_f==-1))
     if(nb_tot_lig>0)then
        !This table has to be totally free because when tab_info_lig%number/=tab_info_lig%max the data can be stored anywhere in the table
        tab_result(:,:,:)%number=0

        !A first record of the data is done: only the molecules with the good amount of ligands and possible good bonds are recorded
        call record(nb_at,nb_lig,max_lig_same,tab_center(1),tab_lig(1,1),&
             a,b,c,tabx(1),taby(1),tabz(1),tab_result(1,1,1),tab_info_lig(1))


        !If there are two ligands with the same symbol and bonds distance which are more or less the same (r1min<r2min<r1max<r2max or r1min<r2min<r2max<r1max) some false positive errors can be done. The permutation subroutine will test all the combinations to know if an arrangement of data can satisfy tab_info_lig(1)%number, tab_info_lig(2)%number, ...
        do i=1,nb_at 
           if(tab_result(1,1,i)%number/=0)then
              do j=1,nb_lig
                 if(  (tab_info_lig(j)%number/=tab_info_lig(j)%max).and.&!These two values are different only if there are two ligands with close properties (same symbol and close bonds)
                      (.not. any(tab_info_lig(:j-1)%symbol==tab_info_lig(j)%symbol)))then!We do not need to take into account a symbol which was already treated
                    k=0!k stands for the number of problematic ligands with the symbol tab_info_lig(j)%symbol

                    do l=1,tab_info_lig(j)%max
                       if(tab_result(l,j,i)%number/=0)then
                          k=k+1
                          !Information to do the permutations and to know if we can find a good arrangement when there are problematic bonds
                          tab_permut(k)=tab_result(l,j,i)
                       endif
                    enddo
                    !tab_r_ref contains the information about the r_min and r_max expected and their amount
                    tab_r_ref(1,1:tab_info_lig(j)%number)=tab_info_lig(j)%rmin
                    tab_r_ref(2,1:tab_info_lig(j)%number)=tab_info_lig(j)%rmax
                    n=tab_info_lig(j)%number+1!n is just an intermediate variable to fill tab_r_ref

                    do m=j+1,nb_lig
                       if(  (tab_info_lig(m)%number/=tab_info_lig(m)%max).and.&!These two values are different only if there are two ligands with close properties (same symbol and close bonds)
                            (tab_info_lig(m)%symbol==tab_info_lig(j)%symbol))then!The symbol of i and j atoms are the same
                          do l=1,tab_info_lig(m)%max
                             if(tab_result(l,m,i)%number/=0)then
                                k=k+1
                                tab_permut(k)=tab_result(l,m,i)
                             endif
                          enddo
                          !tab_r_ref contains the information about the r_min and r_max expected and their amount
                          tab_r_ref(1,n:n+tab_info_lig(m)%number-1)=tab_info_lig(m)%rmin
                          tab_r_ref(2,n:n+tab_info_lig(m)%number-1)=tab_info_lig(m)%rmax
                          n=n+tab_info_lig(m)%number!n is just an intermediate variable to fill tab_r_ref
                       endif
                    enddo!Now in tab_permut there are the k ligands which can cause problems


                    position=1
                    valid=.false.
                    call permutation(position,k,nb_tot_lig,valid,tab_r_ref(1,1),tab_permut(1),tab_order(1))
                    
                    if(valid)then!If there is a good combination to satisfy the bond length we store change the order of the data
                       m=j
                       n=1
                       do l=1,k
                          if(n>tab_info_lig(m)%number)then!We have to convert the "k" value into the "n" atom of the "m" ligand
                             m=m+1
                             do while (tab_info_lig(m)%symbol/=tab_info_lig(j)%symbol)
                                m=m+1
                             enddo
                             n=1
                          endif
                          !tab_result contains the information (number, x, y, z, r) about the ligands of the atoms for a specific step
                          !tab_permut contains the information (number, x, y, z, r) about the problematic ligands for a specific ligand
                          !tab_order is used to reorganize tab_result (if necessary) if there are problematic ligands
                          tab_result(n,m,i)=tab_permut(tab_order(l))
                          n=n+1
                       enddo
                    else!In the other case we do not take into account the atom
                       tab_result(1,1,i)%number=0
                       exit
                    endif

                 end if!If there are ligands with close properties
              enddo!loop for the ligands j
           endif
        enddo!Loop for the central atoms i


        !The final test is done only if the user does not want two central atoms which share the same ligand
        if(.not. share)call sharing(nb_at,nb_lig,max_lig_same,tab_result(1,1,1),tab_info_lig(1))

     endif 


     !=====================
     !Writing the new files
     !=====================
     do i=11,13+nb_tot_lig*3
        write(i,'(i10)',advance='no')istep
     enddo
     do i=1,nb_at
        if(tab_center(i))then
           valid=.true.
           if(nb_tot_lig>0)then!It must not be done if tab_result is not allocated
              !Study of the molecule lifetime
              call lifetime(nb_at,nb_lig,nb_tot_lig,max_lig_same,tab_start_life(1,1)&
                   ,istep,tab_life(1,1,1),tab_info_lig(1),tab_result(1,1,1))

              !The central atom has not the good ligands
              if(tab_result(1,1,i)%number==0)valid=.false.
           endif

           if(valid)then!If the central atom has the good ligands (or if we just study an atom)
              
              !Position of the central atom with pbc
              tabx(i)=tabx(i)-floor(tabx(i)/a)*a
              taby(i)=taby(i)-floor(taby(i)/b)*b
              tabz(i)=tabz(i)-floor(tabz(i)/c)*c
              write(11,'(f10.5)',advance='no')tabx(i)
              write(12,'(f10.5)',advance='no')taby(i)
              write(13,'(f10.5)',advance='no')tabz(i)

              iligand=0
              do j=1,nb_lig
                 do k=1,tab_info_lig(j)%number
                    !Vector M-X
                    write(14+3*iligand,'(f10.5)',advance='no')tab_result(k,j,i)%x
                    write(15+3*iligand,'(f10.5)',advance='no')tab_result(k,j,i)%y
                    write(16+3*iligand,'(f10.5)',advance='no')tab_result(k,j,i)%z
                    iligand=iligand+1
                 enddo
              enddo
           else
              do ifile=11,13+3*nb_tot_lig
                 write(ifile,'(a)',advance='no')"         x"
              enddo
              statistics=statistics+1
           endif
        endif
     enddo
     do ifile=11,13+3*nb_tot_lig
        write(ifile,*)
     enddo

     !========================
     !Reading of the next step
     !========================
     read(10,*,iostat=end_file)nb_at
     if(end_file/=0) exit
     read(10,*)char,char,istep
     if(mod(istep,1000)==0)write(*,*)"Step",istep," treated."
     do i=1,nb_at
        read(10,*)tab_symb(i),tabx(i),taby(i),tabz(i)
     enddo

  enddo


  !Adjustment to have istep="the true final step"
  if((istep>step_f).and.(step_f/=-1))istep=step_f
  
  if(nb_tot_lig>0)then
     !Quantity of good central atoms
     write(*,*)(1.-dble(statistics)/(count(tab_center)*(istep-step_i+1)))*100,"% of ",trim(symb_center)," have the good ligands."

     !We write one more line in the files about lifetime when we reach the end of the record
     do i=1,nb_at
        if(tab_life(1,1,i)/=0)then

           write(14+nb_tot_lig*3,'(i10)',advance='no')i
           do j=1,nb_lig
              do k=1,tab_info_lig(j)%number
                 write(14+nb_tot_lig*3,'(i10)',advance='no')tab_life(k,j,i)
              enddo
           enddo
           write(14+nb_tot_lig*3,'(i10)',advance='no')tab_start_life(1,i)
           write(14+nb_tot_lig*3,'(i10)')istep
           write(15+nb_tot_lig*3,'(i10)',advance='no')i
           do j=1,nb_lig
              do k=1,tab_info_lig(j)%number
                 write(15+nb_tot_lig*3,'(i10)',advance='no')tab_life(k,j,i)
              enddo
           enddo
           write(15+nb_tot_lig*3,'(i10)',advance='no')tab_start_life(2,i)
           write(15+nb_tot_lig*3,'(i10)')istep

        endif
     enddo
     !AND we write the value of the final step
     call fseek(14+nb_tot_lig*3,pos_file(1),0)
     write(14+nb_tot_lig*3,*)istep
     call fseek(15+nb_tot_lig*3,pos_file(2),0)
     write(15+nb_tot_lig*3,*)istep
  endif

  do ifile=10,max_file
     close(ifile)
  enddo
  deallocate(tab_center,tab_symb,tabx,taby,tabz)
  if(nb_tot_lig>0)deallocate(tab_lig,tab_info_lig,tab_permut,tab_result,tab_start_life,tab_life,tab_r_ref,tab_order)


end program main





! 子程序
subroutine atom_selection(nb_at,tab_at,tab_symb,symbol)

  implicit none

  integer nb_at
  logical :: tab_at(nb_at)
  character(len=3)symbol
  character(len=3) :: tab_symb(nb_at)

  !local
  integer i,end_line
  character(len=200)type_ref

  write(*,*)"How do you want to select these atoms? (ALL/NUMBER)"
  read(*,*)type_ref

  if(type_ref=="ALL")then
     do i=1,nb_at
        if(symbol==trim(tab_symb(i)))then
           tab_at(i)=.true.
        endif
     enddo

  elseif(type_ref=="NUMBER")then
     write(*,*)"Write the list of the atoms separated by a coma."
     end_line=0
     do while(end_line==0)
        read(*,'(i30)',advance='no',iostat=end_line)i
        tab_at(i)=.true.
     enddo

  endif

end subroutine atom_selection





subroutine record(nb_at,nb_lig,max_lig_same,tab_center,tab_lig,a,b,c,tabx,taby,tabz,tab_result,tab_info_lig)

  use derived_type
  implicit none

  integer nb_at,nb_lig,max_lig_same
  logical, dimension(nb_at) :: tab_center
  logical, dimension(nb_at,nb_lig) :: tab_lig
  double precision a,b,c
  double precision, dimension(nb_at) :: tabx,taby,tabz
  type (info_lig), dimension(nb_lig) :: tab_info_lig
  type (info_result), dimension(max_lig_same,nb_lig,nb_at) :: tab_result

  !local
  integer i,j,k,bond,iligand
  double precision x,y,z,r

  do i=1,nb_at
     iligand=0!Total amount of ligands surrounding the ith atom
     if(tab_center(i))then!If the ith atom is a central atom
        ligand: do j=1,nb_lig
           bond=0!the values of bond are [1;tab_info_lig(j)%max]
           do k=1,nb_at
              if(tab_lig(k,j))then!if the kth atom is a ligand
                 if(.not. any(tab_result(:,:,i)%number==k))then!... and if this atoms has not been already recorded
                    
                    !Applying pbc
                    x=tabx(k)-tabx(i)
                    x=x-nint(x/a)*a
                    y=taby(k)-taby(i)
                    y=y-nint(y/b)*b
                    z=tabz(k)-tabz(i)
                    z=z-nint(z/c)*c

                    !Recording molecules
                    r=sqrt(x**2+y**2+z**2)
                    if(r<=tab_info_lig(j)%rmax)then
                       bond=bond+1
                       !If the central atom has exactly the good number of ligand with the good distance then we save this molecule
                       if((r>=tab_info_lig(j)%rmin).and.(bond<=tab_info_lig(j)%max))then
                          iligand=iligand+1
                          tab_result(bond,j,i)%number=k
                          tab_result(bond,j,i)%x=x
                          tab_result(bond,j,i)%y=y
                          tab_result(bond,j,i)%z=z
                          tab_result(bond,j,i)%r=r
                       else!If there are too many bonds or if one bond is too short, we do not take into account the molecule
                          tab_result(1,1,i)%number=0!Remark I do not write 0 tab_result(bond,j,i)%number=0 because in the code I just test the value of tab_result(1,1,i)%number
                          exit ligand ! I do not test the other ligands, I pass to the next central atom
                       endif
                    endif
                 endif
              endif

           enddo!End of loop testing all the atoms
        enddo ligand !End of loop passing the ligand
        if(iligand/=sum(tab_info_lig(:)%number))then
           tab_result(1,1,i)%number=0!If not enough M-X bond as been recorded, then we do not take the central atom into account
        endif
        !At this step, we have bond<=tab_info_lig(j)%max for the accepted molecules and the number of ligand is good.        
     endif
  enddo!End of loop pasing the O


end subroutine record



!Subroutine to do permutations
!This subroutine is used only if there are two ligands with the same symbol and bonds distance which are more or less the same: r1min<r2min<r1max<r2max or r1min<r2min<r2max<r1max
recursive subroutine permutation(position,position_max,nb_tot_lig,valid,tab_r_ref,tab_permut,tab_order)

  use derived_type
  implicit none

  integer position,position_max,nb_tot_lig
  integer, dimension(nb_tot_lig) :: tab_order
  logical valid
  double precision, dimension(2,nb_tot_lig) :: tab_r_ref
  type (info_result), dimension(nb_tot_lig) :: tab_permut
  !local
  integer value

  if (position>position_max)then!Once a permutation is done, we test if it is possible to arrange the values to fit the bond length conditions
     
     if(  all(tab_permut(tab_order(:position_max))%r>=tab_r_ref(1,:position_max)).and.&
          all(tab_permut(tab_order(:position_max))%r<=tab_r_ref(2,:position_max)))then!If this order is good
        valid=.true.
     endif

  else
     !A new permutation is doing
     do value = 1,position_max
        if (.not. any(tab_order(:position-1)==value)) then
           tab_order(position)= value
           call permutation (position+1,position_max,nb_tot_lig,valid,tab_r_ref(1,1),tab_permut(1),tab_order(1))
           if(valid)exit
        end if
     end do
  end if!if(pos>pos_max)
  
end subroutine permutation



subroutine sharing(nb_at,nb_lig,max_lig_same,tab_result,tab_info_lig)

  use derived_type
  implicit none

  integer nb_at,nb_lig,max_lig_same
  type(info_lig),dimension(nb_lig) :: tab_info_lig
  type (info_result), dimension(max_lig_same,nb_lig,nb_at) :: tab_result

  !local
  integer i,j,k,l,m,n
  logical valid
  
  !Removing the molecules which share the same ligand
  do i=1,nb_at-1
     if(tab_result(1,1,i)%number/=0)then
        valid=.true.
        do j=i+1,nb_at
           if(tab_result(1,1,j)%number/=0)then
              ligandk: do k=1,nb_lig
                 do l=1,tab_info_lig(k)%number
                    do m=1,nb_lig
                       do n=1,tab_info_lig(m)%number
                          if(tab_result(l,k,i)%number==tab_result(n,m,j)%number)then
                             tab_result(1,1,j)%number=0!Remark I do not write 0 tab_result(m,n,j)%number=0 because in the code I just test the value of tab_result(1,1,j)%number
                             valid=.false.
                             exit ligandk
                          end if
                       enddo!ligand m number n
                    enddo!ligand m
                 enddo!ligand k number l
              enddo ligandk !ligand k
           endif
        enddo!center j
        !We have to remove the ith atom at the end in the case where one ligand is shared by 3, 4, ... atoms
        if(.not.valid)then
           tab_result(1,1,i)%number=0
        endif
     endif
  enddo!center i

end subroutine sharing





subroutine lifetime(nb_at,nb_lig,nb_tot_lig,max_lig_same,tab_start_life,istep,tab_life,tab_info_lig,tab_result)

  use derived_type
  implicit none

  integer nb_at,nb_lig,nb_tot_lig,max_lig_same,istep
  integer, dimension(2,nb_at) :: tab_start_life
  integer,dimension(max_lig_same,nb_lig,nb_at) :: tab_life
  type(info_lig),dimension(nb_lig) :: tab_info_lig
  type (info_result),dimension(max_lig_same,nb_lig,nb_at) :: tab_result
  !local
  integer i,j,k
  logical valid_under,valid_over

  do i=1,nb_at
     if(tab_life(1,1,i)/=0)then!Condition to know if this central atom had good ligands during the previous step
        valid_over=.true.
        valid_under=.true.
        
        ligand:do j=1,nb_lig
           do k=1,tab_info_lig(j)%number
              !Two files about lifetime will be created:
              !-One which can be overestimated (we do not do any difference between all the ligands of the molecule)
              !-One which can be underestimated (we do a difference between all the kind of ligands. It can be problematic if two ligands with the same symbol can be associated with two different bonds)
              if(.not.any(tab_result(:,:,i)%number==tab_life(k,j,i)))then
                 valid_over=.false.
                 valid_under=.false.
                 exit ligand
              elseif(.not.any(tab_result(:,j,i)%number==tab_life(k,j,i)))then
                 valid_under=.false.
              endif
           enddo!k
        enddo ligand!j

        if(.not.valid_over)then!If the ligands change we write the information about the litfe time of the molecule
           write(14+nb_tot_lig*3,'(i10)',advance='no')i
           do j=1,nb_lig
              do k=1,tab_info_lig(j)%number
                 write(14+nb_tot_lig*3,'(i10)',advance='no')tab_life(k,j,i)
              enddo
           enddo
           write(14+nb_tot_lig*3,'(i10)',advance='no')tab_start_life(1,i)
           write(14+nb_tot_lig*3,'(i10)')istep-1
           write(15+nb_tot_lig*3,'(i10)',advance='no')i
           do j=1,nb_lig
              do k=1,tab_info_lig(j)%number
                 write(15+nb_tot_lig*3,'(i10)',advance='no')tab_life(k,j,i)
              enddo
           enddo
           write(15+nb_tot_lig*3,'(i10)',advance='no')tab_start_life(2,i)
           write(15+nb_tot_lig*3,'(i10)')istep-1
           tab_life(:,:,i)=tab_result(:,:,i)%number
           tab_start_life(:,i)=istep
        elseif(.not.valid_under)then
           write(15+nb_tot_lig*3,'(i10)',advance='no')i
           do j=1,nb_lig
              do k=1,tab_info_lig(j)%number
                 write(15+nb_tot_lig*3,'(i10)',advance='no')tab_life(k,j,i)
              enddo
           enddo
           write(15+nb_tot_lig*3,*)" ",istep-1
           tab_life(:,:,i)=tab_result(:,:,i)%number
           tab_start_life(2,i)=istep
        endif

     elseif(tab_result(1,1,i)%number/=0)then!True if the central atom had no ligand but if it has now
        tab_life(:,:,i)=tab_result(:,:,i)%number
        tab_start_life(:,i)=istep
     endif

  enddo!i

end subroutine lifetime
