program velocities
  !autocorrelation of velocities

  implicit none


  interface
     subroutine vdos(nmovie,natoms,nvdos,ivdos,tab_vdos,vx,vy,vz,corr,tab_filename,tab_proj)
       integer, intent(in) :: nmovie,natoms,nvdos,ivdos
       integer, intent(in) :: tab_vdos
       double precision,intent(in) :: vx,vy,vz
       double precision,intent(in) :: corr
       character(len=200),intent(in) :: tab_filename
       character(len=3), intent(in)::tab_proj
     end subroutine vdos
  end interface


  !Definition of a structure to store the informations about common elements (already used atoms) and new elements (to define when the program runs)
  type info_atoms
     character(len=3) name
     double precision mw
     integer amount
  end type info_atoms


  logical file_exist
  integer nmax,nmovie_max,nmovie,natoms,natoms_vdos,nvdos,kindmax,step_min,step_max,end_file
  parameter(nmax=500)        ! max number of atoms
  parameter(nmovie_max = 10000) ! max number of movie steps
  integer i,j,iatom,imovie,ivdos
  integer, dimension(:,:), allocatable:: tab_vdos
  double precision amass_tot
  double precision, dimension(:), allocatable :: corr,vcmx,vcmy,vcmz !Use a 2D table vcm(3,nmovie_max) instead of vcmx, vcmy, vcmz does not change the calculation time
  double precision, dimension(:), allocatable :: fmass
  double precision, dimension(:,:), allocatable :: vx,vy,vz !Use a 3D table v(3,nmax,nmovie_max) instead of vx,vy and vz is not efficient
  character(len=200) type,filetr
  character(len=200), dimension(:), allocatable :: tab_filename
  character(len=3), dimension(:), allocatable :: atom_type,tab_proj
  type (info_atoms), dimension(:), allocatable :: tab_info,tab_tmp 
  character char
  
  step_min=0
  step_max=0
  nmovie=0
  nvdos=0

  !==============================
  !Information about common atoms
  !==============================
  kindmax=11
  allocate (tab_info(kindmax))
  tab_info(1)%name='H'; tab_info(1)%mw=1.0d0 ; tab_info(1)%amount=0
  tab_info(2)%name='Li'; tab_info(2)%mw=7.0d0 ; tab_info(2)%amount=0
  tab_info(3)%name='C'; tab_info(3)%mw=12.0d0 ; tab_info(3)%amount=0
  tab_info(4)%name='N'; tab_info(4)%mw=14.0d0 ; tab_info(4)%amount=0
  tab_info(5)%name='O'; tab_info(5)%mw=16.0d0 ; tab_info(5)%amount=0
  tab_info(6)%name='F'; tab_info(6)%mw=19.0d0 ; tab_info(6)%amount=0
  tab_info(7)%name='Na'; tab_info(7)%mw=22.99d0; tab_info(7)%amount=0
  tab_info(8)%name='S'; tab_info(8)%mw=32.1d0 ; tab_info(8)%amount=0
  tab_info(9)%name='K'; tab_info(9)%mw=39.1d0 ; tab_info(9)%amount=0
  tab_info(10)%name='Ca'; tab_info(10)%mw=40.1d0 ; tab_info(10)%amount=0
  tab_info(11)%name='I' ; tab_info(11)%mw=126.9d0; tab_info(11)%amount=0
  

  !=============================================
  !Reading of the parameters from standard input
  !=============================================
  write(6,*)'Name of the velocity file (default: traj_vel.xyz):'
  read(5,'(a)')filetr
  filetr=adjustl(filetr)
  if(trim(filetr)=='')filetr='traj_vel.xyz'
  inquire(file=filetr, exist=file_exist)
  if(.not. file_exist)then 
     write(6,*)'!!!File ', trim(filetr),' does not exist!!!'
     stop
  endif
  open(10,file=filetr,action='read')
  read(10,*)natoms
  
  
  write(6,*)'First step (default: 0):'
  read(5,'(i30)')step_min
  write(6,*)'Final step (default: maximal size allowed by the amount of data):'
  read(5,'(i30)')step_max
  
  
  write(6,*)'Number of VDOS to do (default a VDOS of all the atoms will writen in vdos.dat):'
  read(5,'(i30)')nvdos
  if(nvdos==0)then!A normal full VDOS is done
     nvdos=1
     allocate(tab_filename(1),tab_vdos(natoms,1),tab_proj(1))
     tab_filename(1)='vdos.dat'
     tab_vdos=1
     tab_proj(1)='no'
  else!In this case, we need the filename, the atoms concerned (ALL or specific) and the projection
     allocate(tab_filename(nvdos),tab_vdos(natoms,nvdos),tab_proj(nvdos))
     tab_vdos=0
     do ivdos=1,nvdos
        write(6,*)'VDOS number: ',ivdos
        write(6,'(a)',advance='no')'  -Filename (default vdos.dat): '
        read(5,'(a)')tab_filename(ivdos)
        tab_filename(ivdos)=adjustl(tab_filename(ivdos))
        if(trim(tab_filename(ivdos))=='')tab_filename(ivdos)='vdos.dat'

        write(6,'(a)',advance='no')'  -Atoms selection type ALL / NUMBER (default ALL): '
        read(5,'(a)')type
        type=adjustl(type)
        if(trim(type)=='')type='ALL'
        if(type=='ALL')then!The VDOS deals with all atoms
           do iatom=1,natoms 
              tab_vdos(iatom,ivdos)=1
           enddo
        elseif(type=='NUMBER')then!The VDOS of specific atoms is done
           write(6,*)'How many atoms are concerned?'
           read(5,'(i30)')natoms_vdos
           write(6,*)'Write the number associated with each atom separated by coma.'
           do iatom=1,natoms_vdos-1 
              read(5,'(i30)',advance='no')i
              tab_vdos(i,ivdos)=1
           enddo
           !It is necessary to read the last line without advance no
           read(5,'(i30)')i
           tab_vdos(i,ivdos)=1
        endif

        write(6,'(a)',advance='no')'  -VDOS projection no / x / y / z (default no projection):'
        read(5,'(a)')tab_proj(ivdos)
        tab_proj(ivdos)=adjustl(tab_proj(ivdos))
        if(trim(tab_proj(ivdos))=='')tab_proj(ivdos)='no'
     enddo
  endif


  !===========================================
  !Tables allocation and velocity file reading
  !===========================================
  if(step_max>0)nmovie=step_max-step_min+1
  if(nmovie*natoms>48000000)then
     write(6,*)'Be carefull, a dynamic table with a size over 3GB will be created!'
  elseif(nmovie<=1) then
     nmovie=floor(16000000./natoms)
     !nmovie=floor(85000000./natoms) !A dynamic table of 5.7GB (approximatively 180,000 steps with 500 atoms) can be created
     write(6,*)'MAX value not precised.'
     write(6,*)'So, MAX value fixed to ',nmovie,'.'
     write(6,*)'A dynamic table of 1GB will be created.'
  endif

  allocate(fmass(natoms),atom_type(natoms))
  allocate(vx(natoms,nmovie),vy(natoms,nmovie),vz(natoms,nmovie))
  write(6,*)'End of INPUT reading'
  write(6,*)
  write(6,*)'Number of atoms in molecule:',natoms
  write(6,*)'Maximal number of movie steps:',nmovie!The actual value of nmovie is >= to the true value

  imovie=1
  end_file=0
  do while (end_file==0)
     if(mod(imovie,1000)==0)write(6,*)'imovie',imovie
     read(10,*,iostat=end_file)char,char,i

     if((step_max==0).and.(step_min<i))then!If step_max is not precised by user, nmovie steps will be done from the beggining
        step_max=i+nmovie-1
     elseif(step_max==0)then!If step_max is not precised by user, nmovie steps will be done from step_min
        step_max=step_min+nmovie-1
     elseif ((imovie==1).and.(step_max<i))then!If the step_max value precised by the user is too low the program stops
        write(6,*)
        write(6,*)'!!!MAX value (final step) too low!!!'
        write(6,*)'!!!The velocity file starts at ',i,'!!!'
        stop
     endif

     if((i>=step_min).and.(i<=step_max))then
        do iatom = 1,natoms
           read(10,*)atom_type(iatom),vx(iatom,imovie),vy(iatom,imovie),vz(iatom,imovie)
        enddo
        imovie=imovie+1
     elseif (i<step_min)then
        do iatom = 1,natoms
           read(10,*)
        enddo
     else
        exit
     endif

     read(10,*,iostat=end_file)
  enddo


  nmovie=imovie-1!Here is the good value of nmovie
  allocate(corr(nmovie),vcmx(nmovie),vcmy(nmovie),vcmz(nmovie))
  close(10)


  write(*,*)'Number of steps recorded: ',nmovie
  write(6,*)'End of TRAJECTORY reading'
  write(6,*)


  !=====================================================================
  !Remove ensemble translation from velocities of atoms 
  !Problems could come from the standard deviation which is pretty high:
  !- vx,vy,vz = 1E-4 ua
  !- drift average = 1E-6 ua
  !- drift std. dev = 1E-5 ua
  !=====================================================================

  amass_tot = 0.d0
  do iatom = 1,natoms
     do i=1,kindmax
        if(atom_type(iatom).eq.tab_info(i)%name)then 
           amass_tot = amass_tot + tab_info(i)%mw
           fmass(iatom) = tab_info(i)%mw
           tab_info(i)%amount=tab_info(i)%amount+1
           exit
        endif

        !Procedure if an non-referenced atom is present
        if(i==kindmax)then
           
           allocate (tab_tmp(kindmax))
           tab_tmp=tab_info
           deallocate (tab_info)
           kindmax=kindmax+1
           allocate (tab_info(kindmax))

           do j=1,kindmax-1
              tab_info(j)=tab_tmp(j)
           enddo

           deallocate(tab_tmp)

           tab_info(kindmax)%name=atom_type(iatom)
           tab_info(kindmax)%amount=1
           write(6,*)'What is the molar weight of ', atom_type(iatom),'?'
           read(5,*) tab_info(kindmax)%mw
           amass_tot = amass_tot + tab_info(kindmax)%mw
           fmass(iatom) = tab_info(kindmax)%mw

           exit
        endif
     enddo   
  enddo

  do iatom= 1,natoms
     fmass(iatom)=fmass(iatom)/amass_tot
  enddo
  
  do i=1,kindmax
     if(tab_info(i)%amount>0)write(6,*)'M(',tab_info(i)%name,')=',tab_info(i)%mw,' g.mol-1'
  enddo
  write(6,*)'There are:'
  do i=1,kindmax
     if(tab_info(i)%amount>0)write(6,*)tab_info(i)%amount,tab_info(i)%name,&
     ' %w=',100*tab_info(i)%amount*tab_info(i)%mw/amass_tot
  enddo
  write(6,*)'Total weight',amass_tot
  write(6,*)

  do imovie = 1,nmovie
     vcmx(imovie) = 0.d0
     vcmy(imovie) = 0.d0
     vcmz(imovie) = 0.d0
     do iatom = 1,natoms
        vcmx(imovie) = vcmx(imovie) + fmass(iatom)*vx(iatom,imovie)
        vcmy(imovie) = vcmy(imovie) + fmass(iatom)*vy(iatom,imovie)
        vcmz(imovie) = vcmz(imovie) + fmass(iatom)*vz(iatom,imovie)
     enddo

     do iatom=1,natoms
        vx(iatom,imovie) = vx(iatom,imovie) - vcmx(imovie)
        vy(iatom,imovie) = vy(iatom,imovie) - vcmy(imovie)
        vz(iatom,imovie) = vz(iatom,imovie) - vcmz(imovie)
     enddo

  enddo
  
  !==========================================
  !The ensemble translation has been removed,
  !the vdos can be calculate.
  !==========================================
  do ivdos=1,nvdos
     call vdos(nmovie,natoms,nvdos,ivdos,tab_vdos(1,1),vx(1,1),vy(1,1),vz(1,1),corr(1),tab_filename(1),tab_proj(1))
  enddo


  deallocate(tab_info,fmass,atom_type,tab_vdos,tab_filename)
  deallocate(corr,vcmx,vcmy,vcmz)
  deallocate(vx,vy,vz)

end program velocities






!VDOS calculation
subroutine vdos(nmovie,natoms,nvdos,ivdos,tab_vdos,vx,vy,vz,corr,tab_filename,tab_proj)
  implicit none

  !global
  integer nmovie,natoms,nvdos,ivdos
  integer,dimension(natoms,nvdos)::tab_vdos
  double precision,dimension(natoms,nmovie):: vx,vy,vz
  double precision,dimension(nmovie)::corr
  character(len=200), dimension(nvdos)::tab_filename
  character(len=3), dimension(nvdos)::tab_proj

  !local
  integer iatom,mt,nts
  double precision scalar,norm
  
  open(10,file=tab_filename(ivdos))

  select case (tab_proj(ivdos))
  !=============
  !No projection
  !=============
  case('no')
     do mt = 0,nmovie-1
        corr(mt+1) = 0.d0
        do nts = 1,nmovie-mt
           scalar = 0.d0

           do iatom=1,natoms
              if(tab_vdos(iatom,ivdos)==1)then!All the atoms have to be tested even if only one is studied so it is time consuming
                 scalar = scalar + vx(iatom,nts)*vx(iatom,nts+mt) + vy(iatom,nts)*vy(iatom,nts+mt)&
                      + vz(iatom,nts)*vz(iatom,nts+mt) ! scalar product of velocities
              endif
           enddo
           corr(mt+1) = corr(mt+1) + scalar
        enddo
        if(mt==0)norm=corr(1)/dble(nmovie)
        write(10,*)mt,corr(mt+1)/dble(nmovie-mt)/norm
        if(modulo(mt,1000)==0 .and. (mt/=0))write(6,*)mt,' steps already writen in ', trim(tab_filename(ivdos))
     enddo

  !============
  !x projection
  !============
  case('x')
     do mt = 0,nmovie-1
        corr(mt+1) = 0.d0
        do nts = 1,nmovie-mt
           scalar = 0.d0

           do iatom=1,natoms
              if(tab_vdos(iatom,ivdos)==1)then!All the atoms have to be tested even if only one is studied so it is time consuming
                 scalar = scalar + vx(iatom,nts)*vx(iatom,nts+mt)! scalar product of velocities
              endif
           enddo
           corr(mt+1) = corr(mt+1) + scalar
        enddo
        if(mt==0)norm=corr(1)/dble(nmovie)
        write(10,*)mt,corr(mt+1)/dble(nmovie-mt)/norm
        if(modulo(mt,1000)==0 .and. (mt/=0))write(6,*)mt,' steps already writen in ', trim(tab_filename(ivdos))
     enddo

  !============
  !y projection
  !============
  case('y')
     do mt = 0,nmovie-1
        corr(mt+1) = 0.d0
        do nts = 1,nmovie-mt
           scalar = 0.d0

           do iatom=1,natoms
              if(tab_vdos(iatom,ivdos)==1)then!All the atoms have to be tested even if only one is studied so it is time consuming
                 scalar = scalar + vy(iatom,nts)*vy(iatom,nts+mt)! scalar product of velocities
              endif
           enddo
           corr(mt+1) = corr(mt+1) + scalar
        enddo
        if(mt==0)norm=corr(1)/dble(nmovie)
        write(10,*)mt,corr(mt+1)/dble(nmovie-mt)/norm
        if(modulo(mt,1000)==0 .and. (mt/=0))write(6,*)mt,' steps already writen in ', trim(tab_filename(ivdos))
     enddo

  !============
  !z projection
  !============
  case('z')
     do mt = 0,nmovie-1
        corr(mt+1) = 0.d0
        do nts = 1,nmovie-mt
           scalar = 0.d0

           do iatom=1,natoms
              if(tab_vdos(iatom,ivdos)==1)then!All the atoms have to be tested even if only one is studied so it is time consuming
                 scalar = scalar + vz(iatom,nts)*vz(iatom,nts+mt)! scalar product of velocities
              endif
           enddo
           corr(mt+1) = corr(mt+1) + scalar
        enddo
        if(mt==0)norm=corr(1)/dble(nmovie)
        write(10,*)mt,corr(mt+1)/dble(nmovie-mt)/norm
        if(modulo(mt,1000)==0 .and. (mt/=0))write(6,*)mt,' steps already writen in ', trim(tab_filename(ivdos))
     enddo

  end select

  close(10)
  write(6,*)'Correlation results written in file ',trim(tab_filename(ivdos))

  return
end subroutine vdos
