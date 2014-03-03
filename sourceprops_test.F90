module sourceprops

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR, YEAR
  use cosmology_parameters, only: Omega_B, Omega0
  use sizes, only: Ndim, mesh ! WW added
  use nbody, only:  dir_src
  use material, only: xh
  use grid, only: x,y,z
  use c2ray_parameters, only: phot_per_atom, lifetime, &
       S_star_nominal, StillNeutral, Number_Sourcetypes,&
       control_rank

  integer :: NumSrc_Glob, NumSrc_Loc
  integer :: srcpos_xyz(3)
  integer :: target_rank,rank_counter
  integer,dimension(:,:),allocatable :: srcpos
  real(kind=dp),dimension(:,:),allocatable :: rsrcpos
  real(kind=dp),dimension(:),allocatable :: NormFlux
  real(kind=dp) :: NormFlux1,Sum_NormFlux
  integer,dimension(:),allocatable :: srcTarget ! WW: target rank for source
  integer,dimension(:),allocatable :: target_count_by_node ! WW: array tracking how many sources each node will have
  integer,dimension(:),allocatable :: tag_array

contains
  
  ! =======================================================================

  subroutine source_properties(zred_now,nz,lifetime2,restart)

    ! Input routine: establish the source properties
    ! Authors: Garrelt Mellema, Ilian Iliev
    ! Update: 30-Jan-2008 (20-Sep-2006 (3-jan-2005, 15-Apr-2004))


    !> WW: I have updated this file to only distribute the sources to 
    !! the node they belong to spatially

    ! For random permutation of sources
    use  m_ctrper

    real(kind=dp),intent(in) :: zred_now ! current redshift
    real(kind=dp),intent(in) :: lifetime2 ! time step
    integer,intent(in) :: nz
    integer,intent(in) :: restart

    character(len=512) :: sourcelistfile,sourcelistfilesuppress
    integer :: ns,ns0

#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPILOG     
    write(logf,*) "Check sourceprops: ",zred_now,nz,lifetime2,restart
#endif 

    ! Deallocate arrays
    if (allocated(srcpos)) deallocate(srcpos)
    if (allocated(rsrcpos)) deallocate(rsrcpos)
    if (allocated(NormFlux)) deallocate(NormFlux)
    if (allocated(srcTarget)) deallocate(srcTarget)
    if (allocated(target_count_by_node)) deallocate(target_count_by_node)
    if (allocated(tag_array)) deallocate(tag_array)


#ifdef MPI
    allocate(target_count_by_node(0:npr-1))
    allocate(tag_array(0:npr-1))
    tag_array = 1
#else
    allocate(target_count_by_node(0:0))
#endif
    target_count_by_node = 0
    Sum_NormFlux = 0

    ! control rank reads in sources
    if (rank == control_rank) then
       ! Construct the file names
       sourcelistfile=trim(adjustl(dir_src))//"test_sources.dat"

       open(unit=50,file=sourcelistfile,status="old")
       ! Number of sources
       read(50,*) NumSrc_Glob

    endif 
    
#ifdef MPI
    ! Distribute source number to all other nodes
    call MPI_BCAST(NumSrc_Glob,1,MPI_INTEGER,control_rank,MPI_COMM_NEW,mympierror)
#endif

    if (NumSrc_Glob > 0) then

      allocate(srcTarget(NumSrc_Glob))
       
      if (rank == control_rank) then
	do ns=1,NumSrc_Glob
	  read(50,*) srcpos_xyz,NormFlux1
#ifdef MPI
	  call find_target_rank(srcpos_xyz(1),srcpos_xyz(2),srcpos_xyz(3))
#else
	  target_rank = 0
#endif
	  srcTarget(ns) = target_rank

	  target_count_by_node(target_rank) = &
	      &target_count_by_node(target_rank) + 1
	enddo

	close(50)

      endif
    else ! if (NumSrc_Glob > 0) 
      NumSrc_Glob = 0
      NumSrc_Loc = 0
    endif

    ! Broadcast source info
#ifdef MPI


    call MPI_BCAST(target_count_by_node,npr,MPI_INTEGER,&
	&control_rank,MPI_COMM_NEW,mympierror)

    NumSrc_Loc = target_count_by_node(rank)

    do rank_counter = 0, npr

      if(rank .eq. 0 .and. rank_counter .eq. 0) then
	write(*,*) 'Source counts distributed across tasks...'
	write(*,*) 'Total sources = ',NumSrc_Glob
      endif

      if(rank .eq. rank_counter) then
	write(*,*) 'Rank = ', rank, 'Source count = ',NumSrc_Loc
      endif
      call MPI_BARRIER(MPI_COMM_NEW,mympierror)
    enddo

#else
    NumSrc_Loc = target_count_by_node(0)
#endif


    if (NumSrc_Glob > 0) then

      ! now allocate source array on each rank

      allocate(srcpos(3,NumSrc_Loc))
      allocate(rsrcpos(3,NumSrc_Loc))
      allocate(NormFlux(NumSrc_Loc))

      ! Control rank reads in sources again to distribute them
      if(rank == control_rank) then

	open(unit=50,file=sourcelistfile,status="old")
	! Number of sources
	read(50,*) NumSrc_Glob

	do ns=1,NumSrc_Glob
	  read(50,*) srcpos_xyz,NormFlux1
	  NormFlux1=NormFlux1/S_star_nominal
	  Sum_NormFlux = Sum_NormFlux + NormFlux1

#ifdef MPI
	  call find_target_rank(srcpos_xyz(1),srcpos_xyz(2),srcpos_xyz(3))

	  if(target_rank .ne. control_rank) then
	    
	    call MPI_SEND(srcpos_xyz, sizeof(srcpos_xyz), &
		&MPI_INTEGER, target_rank,tag_array(target_rank), &
		&MPI_COMM_NEW, mympierror)

	    call MPI_SEND(NormFlux1, 1, MPI_DOUBLE_PRECISION, target_rank,&
	        & tag_array(target_rank)+target_count_by_node(target_rank),&
	        &MPI_COMM_NEW, mympierror)

	  else ! target_rank = control_rank

  	    srcpos(1,tag_array(target_rank)) = srcpos_xyz(1)
	    srcpos(2,tag_array(target_rank)) = srcpos_xyz(2)
	    srcpos(3,tag_array(target_rank)) = srcpos_xyz(3) 
            NormFlux(tag_array(target_rank)) = NormFlux1

	  endif

	  tag_array(target_rank) = tag_array(target_rank) + 1


#else
	  target_rank = 0
	  srcpos(1,ns) = srcpos_xyz(1)
	  srcpos(2,ns) = srcpos_xyz(2)
	  srcpos(3,ns) = srcpos_xyz(3)
	  NormFlux(ns) = NormFlux1
#endif

	enddo

	close(50)	

      else ! (rank == control_rank) then
#ifdef MPI
	do ns = 1,NumSrc_Loc

	  call MPI_RECV(srcpos_xyz, sizeof(srcpos_xyz), MPI_INTEGER, control_rank,&
	     &ns, MPI_COMM_NEW, mympi_status, mympierror)

	  call MPI_RECV(NormFlux1, 1, MPI_DOUBLE_PRECISION, control_rank,&
	     &ns+NumSrc_Loc, MPI_COMM_NEW, mympi_status, mympierror)

          srcpos(1,ns) = srcpos_xyz(1) 
          srcpos(2,ns) = srcpos_xyz(2) 
          srcpos(3,ns) = srcpos_xyz(3) 
	  NormFlux(ns) = NormFlux1

	enddo
#endif
      endif ! of rank 0

      do rank_counter = 0, npr-1

	if(rank .eq. rank_counter) then

	  write(*,*) 'SOURCES ON RANK',rank

	  do ns = 1,NumSrc_Loc
	    write(*,*) srcpos(1:3,ns),NormFlux(ns)
	  enddo

	endif
#ifdef MPI
	call MPI_BARRIER(MPI_COMM_NEW,mympierror)
#endif
      enddo


	! Source is always at cell centre!!
	do ns=1,NumSrc_Loc
	  rsrcpos(1,ns)=x(srcpos(1,ns))
	  rsrcpos(2,ns)=y(srcpos(2,ns))
	  rsrcpos(3,ns)=z(srcpos(3,ns))
	enddo

      if (rank == control_rank) then
	write(logf,*) 'Total flux= ',Sum_NormFlux*S_star_nominal,' s^-1'
	flush(logf)
      endif
    endif
    
  end subroutine source_properties
     
  ! =======================================================================

  !> Initialization routine: dummy
  !! Author: Garrelt Mellema
  
  subroutine source_properties_ini()
    
  end subroutine source_properties_ini

  ! =======================================================================

  !> Determines the appropriate target node for a source
  !! N.B. Currently needs reorder = .TRUE. to be set in 
  !! MPI.F90.
  !!
  !! Author: William Watson
#ifdef MPI
  subroutine find_target_rank(pos_x,pos_y,pos_z)

      integer :: pos_x, pos_y, pos_z
      real(kind=dp) :: cell_test
      real(kind=dp) :: cells_per_node_per_dim1
      real(kind=dp) :: cells_per_node_per_dim2
      real(kind=dp) :: cells_per_node_per_dim3
      real(kind=dp) :: c_per_n_per_d
      integer :: nodes_dim

      nodes_dim = dims(1)

      cells_per_node_per_dim1 = mesh(1)
      cells_per_node_per_dim2 = mesh(2)
      cells_per_node_per_dim3 = mesh(3)

      if(cells_per_node_per_dim1 .ne. cells_per_node_per_dim2 .or.&
	  cells_per_node_per_dim2 .ne. cells_per_node_per_dim3) then
      
	write(*,*) 'IRREGULAR GRID!!! EXITING!'

	call MPI_ABORT(MPI_COMM_NEW,mympierror)

      else

	c_per_n_per_d =  cells_per_node_per_dim1

      endif

!      cell_test =  real(mesh(1))/real(dims(1)) - &
!	  &floor(real(mesh(1))/real(dims(1)))

!      if(abs(cell_test) .gt. 0.001) then

!	write(*,*) 'ERROR - CHOICE OF GRID SIZE NOT &
!	    & COMPATIBLE WITH NUMBER OF MPI NODES.'
!	write(*,*) 'GRID SIZE = ',mesh(1),"Nodes per dim = ",dims(1)

!	call MPI_ABORT(MPI_COMM_NEW,mympierror)

!      endif

      ! find target file for halo:

    if(int(pos_x/c_per_n_per_d).lt.nodes_dim .and.&
	& int(pos_y/c_per_n_per_d).lt.nodes_dim.and.&
	& int(pos_z/c_per_n_per_d).lt.nodes_dim) then

      target_rank = int(pos_x/c_per_n_per_d)+nodes_dim*&
	  &int(pos_y/c_per_n_per_d)+(nodes_dim**2.0)*&
	  int(pos_z/c_per_n_per_d)

    elseif(int(pos_x/c_per_n_per_d).ge.nodes_dim.and.&
	  & int(pos_y/c_per_n_per_d).lt.nodes_dim.and.&
	  & int(pos_z/c_per_n_per_d).lt.nodes_dim) then

      target_rank = (nodes_dim-1)+nodes_dim*&
	  &int(pos_y/c_per_n_per_d)+(nodes_dim**2.0)*&
	  &int(pos_z/c_per_n_per_d)

    elseif(int(pos_x/c_per_n_per_d).lt.nodes_dim.and.&
	  & int(pos_y/c_per_n_per_d).ge.nodes_dim.and.&
	  & int(pos_z/c_per_n_per_d).lt.nodes_dim) then

      target_rank = int(pos_x/c_per_n_per_d)+&
	  &nodes_dim*(nodes_dim-1)+(nodes_dim**2.0)*&
	  &int(pos_z/c_per_n_per_d)

    elseif(int(pos_x/c_per_n_per_d).lt.nodes_dim.and.&
	  & int(pos_y/c_per_n_per_d).lt.nodes_dim.and.&
	  & int(pos_z/c_per_n_per_d).ge.nodes_dim) then

      target_rank = int(pos_x/c_per_n_per_d)+&
	  &nodes_dim*int(pos_y/c_per_n_per_d)+&
	  &(nodes_dim**2.0)*(nodes_dim-1)

    elseif(int(pos_x/c_per_n_per_d).ge.nodes_dim.and.&
	  & int(pos_y/c_per_n_per_d).ge.nodes_dim.and.&
	  & int(pos_z/c_per_n_per_d).lt.nodes_dim) then

      target_rank = (nodes_dim-1)+nodes_dim*&
	  &(nodes_dim-1)+(nodes_dim**2.0)*int(pos_z/c_per_n_per_d)

    elseif(int(pos_x/c_per_n_per_d).ge.nodes_dim.and.&
	  & int(pos_y/c_per_n_per_d).lt.nodes_dim.and.&
	  & int(pos_z/c_per_n_per_d).ge.nodes_dim) then

      target_rank = (nodes_dim-1)+nodes_dim*&
	  &int(pos_y/c_per_n_per_d)+(nodes_dim**2.0)*(nodes_dim-1)

    elseif(int(pos_x/c_per_n_per_d).lt.nodes_dim.and.&
	  & int(pos_y/c_per_n_per_d).ge.nodes_dim.and.&
	  & int(pos_z/c_per_n_per_d).ge.nodes_dim) then

      target_rank = int(pos_x/c_per_n_per_d)+&
	  &nodes_dim*(nodes_dim-1)+(nodes_dim**2.0)*(nodes_dim-1)

    elseif(int(pos_x/c_per_n_per_d).ge.nodes_dim.and.&
	  & int(pos_y/c_per_n_per_d).ge.nodes_dim.and.&
	  & int(pos_z/c_per_n_per_d).ge.nodes_dim) then

      target_rank = (nodes_dim-1)+nodes_dim*&
	  &(nodes_dim-1)+(nodes_dim**2.0)*(nodes_dim-1)


    endif



  end subroutine find_target_rank
#endif

end module sourceprops


