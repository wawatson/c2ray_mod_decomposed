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
       S_star_nominal, StillNeutral, Number_Sourcetypes

  integer :: NumSrc
  integer :: srcpos_x, srcpos_y, srcpos_z
  integer :: target_rank 
  integer,dimension(:,:),allocatable :: srcpos
  real(kind=dp),dimension(:,:),allocatable :: rsrcpos
  real(kind=dp),dimension(:),allocatable :: NormFlux
  integer,dimension(:),allocatable :: srcSeries
  integer,dimension(:),allocatable :: srcTarget ! WW: target rank for source

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
    
    ! Deallocate source arrays
    if (allocated(srcpos)) deallocate(srcpos)
    if (allocated(rsrcpos)) deallocate(rsrcpos)
    if (allocated(NormFlux)) deallocate(NormFlux)
    if (allocated(srcSeries)) deallocate(srcSeries)
    if (allocated(srcTarget)) deallocate(srcTarget)


    ! Rank 0 reads in sources
    if (rank == 0) then
       ! Construct the file names
       sourcelistfile=trim(adjustl(dir_src))//"test_sources.dat"

       open(unit=50,file=sourcelistfile,status="old")
       ! Number of sources
       read(50,*) NumSrc

    endif ! end of rank 0 test
    
#ifdef MPI
    ! Distribute source number to all other nodes
    call MPI_BCAST(NumSrc,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

    ! Allocate arrays for this NumSrc
    if (NumSrc > 0) then
       allocate(srcpos(3,NumSrc))
       allocate(rsrcpos(3,NumSrc))
       allocate(NormFlux(NumSrc))
       allocate(SrcSeries(NumSrc))
       allocate(srcTarget(NumSrc))
       
       ! Fill in the source arrays
       if (rank == 0) then
          do ns=1,NumSrc
             read(50,*) srcpos(1,ns),srcpos(2,ns),srcpos(3,ns),NormFlux(ns)
             NormFlux(ns)=NormFlux(ns)/S_star_nominal
#ifdef MPI
	     call find_target_node(srcpos(1,ns),srcpos(2,ns),srcpos(3,ns))
#else
	     target_rank = 0
#endif
	     srcTarget(ns) = target_rank

          enddo
          close(50)
          
          ! STANDARD TEST
          ! Source positions in file start at 1!
          !           srcpos(1:3,1)=(/ 50, 50, 50 /)
          !           srcpos(1:3,2)=(/ 51, 50, 50 /)
          !           srcpos(1:3,3)=(/ 52, 50, 50 /)
          !           srcpos(1:3,4)=(/ 53, 50, 50 /)
          !           NormFlux(1:4)=1e55_dp/S_star_nominal
          
          !           srcpos(1:3,5)=(/ 20, 10, 10 /)
          !           NormFlux(5)=1e57_dp/S_star_nominal
          
          
          !           srcpos(1:3,6)=(/ 70, 70, 50 /)
          !           srcpos(1:3,7)=(/ 72, 70, 50 /)
          !           srcpos(1:3,8)=(/ 70, 72, 50 /)
          !           srcpos(1:3,9)=(/ 72, 72, 50 /)
          !           NormFlux(6:8)=1e55_dp/S_star_nominal
          !           NormFlux(9)=1e56_dp/S_star_nominal
          
          !           srcpos(1:3,10)=(/ 20, 10, 90 /)
          !           NormFlux(10)=1e54_dp/S_star_nominal
          
          ! Source is always at cell centre!!
          do ns=1,NumSrc
             rsrcpos(1,ns)=x(srcpos(1,ns))
             rsrcpos(2,ns)=y(srcpos(2,ns))
             rsrcpos(3,ns)=z(srcpos(3,ns))
          enddo
       endif ! of rank 0 test
       
       
#ifdef MPI
       ! Distribute the source parameters to the other nodes
       call MPI_BCAST(srcpos,3*NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(rsrcpos,3*NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(NormFlux,NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(srcTarget,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)       
#endif

       if (rank == 0) then
          write(logf,*) 'Total flux= ',sum(NormFlux)*S_star_nominal,' s^-1'

          ! Create array of source numbers for generating random order
          do ns=1,NumSrc
             SrcSeries(ns)=ns
          enddo
          ! Make a random order
          call ctrper(SrcSeries(1:NumSrc),1.0)
       endif
       
#ifdef MPI
       ! Distribute the source series to the other nodes
       call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
     
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

  subroutine find_target_node(pos_x,pos_y,pos_z)

      integer :: pos_x, pos_y, pos_z
      real(kind=dp) :: cell_test
      real(kind=dp) :: cells_per_node_per_dim1
      real(kind=dp) :: cells_per_node_per_dim2
      real(kind=dp) :: cells_per_node_per_dim3
      real(kind=dp) :: c_per_n_per_d
      integer :: nodes_dim

      nodes_dim = dims(1)

      cells_per_node_per_dim1 = mesh(1)/dims(1)
      cells_per_node_per_dim2 = mesh(2)/dims(2)
      cells_per_node_per_dim3 = mesh(3)/dims(3)

      if(cells_per_node_per_dim1 .ne. cells_per_node_per_dim2 .or.&
	  cells_per_node_per_dim2 .ne. cells_per_node_per_dim3) then
      
	write(*,*) 'IRREGULAR GRID!!! EXITING!'

	call MPI_ABORT(MPI_COMM_NEW,mymierror)

      else

	c_per_n_per_d =  cells_per_node_per_dim1

      endif

      cell_test =  real(mesh(1))/real(dims(1)) - &
	  &floor(real(mesh(1))/real(dims(1)))

      if(abs(cell_test) .gt. 0.001) then

	write(*,*) 'ERROR - CHOICE OF GRID SIZE NOT &
	    & COMPATIBLE WITH NUMBER OF MPI NODES.'
	write(*,*) 'GRID SIZE = ',mesh(1),"Nodes per dim = ",dims(1)

	call MPI_ABORT(MPI_COMM_NEW,mymierror)

      endif


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



  end subroutine find_target_node


end module sourceprops


