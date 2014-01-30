program Make_Source_List

  ! reads in halo catalogues and bins them in coarser grid cells
  ! Authors: Garrelt Mellema and Ilian Iliev
  ! Modified from code by Ilian Iliev
  ! Date: 11-Oct-2005
  ! Modified further by Ilian Iliev (2 types of sources, large and
  ! small)
  ! Modified to correct the source positions for the random shifts
  ! done by CubeP3M (II: 21/05/2011)
  ! Random shifts correction not needed anymore for new runs (corrected
  ! inside code now) (II:01/04/2013)

  ! This program constructs lists of collapsed halos and their mass,
  ! using output (halo list) from a CubeP3M simulation. The data is
  ! regridded to 2-3 different lower resolutions.

  implicit none

  ! Where the CubeP3M data lives
  character(len=180),parameter :: data_dir="/ccc/scratch/cont005/pa0442&
       /ilievi/cubepm_130315_6_1728_47Mpc_ext2/results/full_catalogues/"
  !Where the results will be put
  character(len=180),parameter :: results_dir="/ccc/scratch/cont005/pa0442&                                                                                   /ilievi/cubepm_130315_6_1728_47Mpc_ext2/RT_postprocessing"

  ! Cosmological parameters
  real*8,parameter :: omegabh2=0.02156!WMAP5
  real*8,parameter :: omega0=0.27!DM+baryons
  real*8,parameter :: lambda0=0.73
  real*8,parameter :: h=0.70
  real*8,parameter :: an=0.96d0
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  real*8,parameter :: tcmb=2.732
  real*8,parameter :: M_sun = 1.989d33
  real*8,parameter :: cm_to_Mpc=3.086d24
  !  real*8,parameter :: pi=3.1415927
  !  real*8,parameter :: grav_const=6.672d-8 !cm^3/g/s^2
  real*8,parameter :: h0 = h*100.
  real*8,parameter :: omegab = omegabh2/h**2
  real*8,parameter :: rho_crit_0=1.88d-29*h**2 !g/cm^3
  real*8,parameter :: rho_crit_0_MpcM_sun=rho_crit_0*cm_to_Mpc**3&
       &/M_sun
  real*8,parameter :: rho_bar = rho_crit_0_MpcM_sun*omega0

  ! Simulation parameters
  real*8,parameter  :: boxsize=47. ! Mpc h^-1
  integer,parameter :: n_box=3456 ! cells/side in original data
  real(kind=8),parameter :: M_box = rho_bar*(boxsize/h)**3
  real(kind=8),parameter :: M_grid = M_box/(real(n_box))**3

  ! Smaller grid sizes
  !  integer,parameter :: n_coarse=750
  !  integer,parameter :: n_coarser=375
  !integer,parameter :: n_coarse= 
  integer,parameter :: n_coarser=612
  integer,parameter :: n_coarsest=306

  ! Parameters which distinguish between different types of halos
  !  real*8,parameter :: large_halo=8000.
  !  real*8,parameter :: small_halo=800.
  real*8,parameter :: large_halo=1.e9
  real*8,parameter :: small_halo=1.e8 !M_grid

  ! counters
  integer(4) :: numhalo, num_z

  ! Data as found in the PMFAST/CubeP3M halo catalogue
  real(4), dimension(3) :: halo_pos ! halo position (cells)
  real(4), dimension(3) :: x_mean ! centre of mass position 
  real(4), dimension(3) :: v_mean ! velocity
  real(4), dimension(3) :: l      ! angular momentum
  real(4), dimension(3) :: v_disp      ! velocity dispersion                   
  real(4) :: radius_calc ! radius                                              
  real(4) :: halo_mass,halo_mass_uncorrected ! mass calculated on the grid (in grid masses)                                                                  
  real(4) :: imass       ! mass calculated from the particles (in grid masses)
  real(4), dimension(3) :: var      ! some kind of xyz variance of the halo shape                                                                            
  !real(4), dimension(6) :: I_ij     ! some kind of mass matrix of the halo shape                                              

  real(8) :: coll_frac(2) ! collapsed fractions in halos: low-mass
  ! (1) and high-mass(2)

  ! positions in the new grids
  !integer,dimension(3) :: p4 ! coarse
  integer,dimension(3) :: p8 ! coarser
  integer,dimension(3) :: p16 ! coarsest

  ! loop counters and other useful integers
  integer :: i,j,k,lob,ll,halo_type
  integer :: igrid
  integer :: ngrid

  integer, parameter :: max_input=150 !! maximum number of checkpoints
  integer, parameter :: max_halos=176000000 !! maximum number of
  !! halos to read
                                            ! might need increasing!!
  real(4),dimension(max_input) :: z
  real(4) :: z_write

  integer :: num_sources,nz,num_sources_arr(3)

  real(4) :: z_offset, offset(3)

  ! Gridded data: numbers and masses of halos
  !integer,dimension(:,:,:,:),target,allocatable :: n1halo ! coarse 
  integer,dimension(:,:,:,:),target,allocatable :: n2halo ! coarser
  integer,dimension(:,:,:,:),target,allocatable :: n3halo ! coarsest
  integer,dimension(:,:,:,:),pointer :: nhalo
  !real,dimension(:,:,:,:),target,allocatable :: m1halo ! coarse 
  real,dimension(:,:,:,:),target,allocatable :: m2halo ! coarser
  real,dimension(:,:,:,:),target,allocatable :: m3halo ! coarsest
  real,dimension(:,:,:,:),pointer :: mhalo

  ! Various character variables
  integer(kind=4),parameter :: MSL = 180
  character(len=MSL) :: ifile
  character(len=MSL) :: mfile
  character(len=7) :: z_str ! (string(z))
  character(len=10) :: id 

  !
  ! Simulation parameters

  print*,M_grid, M_box, small_halo, large_halo
  ! Allocate grids

  !allocate(m1halo(n_coarse,n_coarse,n_coarse,2))
  allocate(m2halo(n_coarser,n_coarser,n_coarser,2))
  allocate(m3halo(n_coarsest,n_coarsest,n_coarsest,2))
  !allocate(n1halo(n_coarse,n_coarse,n_coarse,2))
  allocate(n2halo(n_coarser,n_coarser,n_coarser,2))
  allocate(n3halo(n_coarsest,n_coarsest,n_coarsest,2))

  ! Open input file with list of redshifts
  open(8,file='../code/input/halofinds')

  ! Open output file for source number statistics
  open(104,file='z_sourcenum_fcoll.dat')

 ! Open input file for random shifts
  !open(105,file='random_offsets_final.dat')

  do num_z=1,max_input
     read(unit=8,err=51,end=41,fmt='(f20.10)') z(num_z)
  enddo
41 num_z=num_z-1
51 close(8)

  print*,'num_z= ',num_z


  ! Loop over redshifts
  redshift:  do nz=1,num_z

     coll_frac=0.0d0

     z_write = z(nz)
     write(*,*) 'processing z=',z_write

     ! Construct file names of halo lists

     write(z_str,'(f7.3)') z_write
     z_str=adjustl(z_str)
     ifile=trim(adjustl(z_str))//"halo.dat"

!     print*,ifile,trim(adjustl(data_dir))//ifile

     ! Open the halo catalogues
     open(unit=31,status='old',file=trim(adjustl(data_dir))//ifile)!
     !,form='formatted')

!     print*,'opened ',trim(adjustl(data_dir))//ifile

     ! Initialize to zero
     !m1halo(:,:,:,:)=0.0d0
     m2halo(:,:,:,:)=0.0d0
     m3halo(:,:,:,:)=0.0d0
     !n1halo(:,:,:,:)=0
     n2halo(:,:,:,:)=0
     n3halo(:,:,:,:)=0

     ! Inform about progress
!     print*,'doing z=',z_write

     !read(105,*) z_offset, offset

     !if(abs(z_write-z_offset)>1e-3)then
     !   print*,' random offsets are not at matching redshifts'
     !   stop
     !end if

     ! Read in catalogues and count sources
     !     read(31) numhalo ! number of halos
     !     print*,'number of halos=', numhalo
     do i=1,max_halos
        read(unit=31,end=133,fmt=*) &!,fmt='(17f20.10)') 
             halo_pos(:), x_mean(:), v_mean(:), &
             l(:), v_disp(:), radius_calc, halo_mass, imass, halo_mass_uncorrected, var(:)!, I_ij(:)
!        print*, halo_pos(:), x_mean(:), v_mean(:), &
!             l(:), v_disp, radius_calc, halo_mass, imass, halo_mass_no_corr
!        pause

        if(halo_mass > 400. .and. imass >0)then
           
           ! Treat different halo masses differently
!           if ( halo_mass*M_grid > small_halo ) then 
              
              if ( halo_mass*M_grid > large_halo ) then 
                 halo_type=2
                 ! these are the larger, never suppressed halos
              else
                 halo_type=1
                 ! these are the smaller, Pop. II/III, suppressed halos
                 ! include only well-resolved halos, with >20 particles
              endif

              coll_frac(halo_type)=coll_frac(halo_type)+halo_mass*M_grid
              !              print*, halo_type
              
              !These checks are for the rare occasions when a halo is in last cell, in which case
              !it goes out of bounds of the coarsened arrays.
!!$              if(halo_pos(1).ge.n_box)then
!!$                 print*,'out of bounds',halo_pos(:), x_mean(:), v_mean(:), &
!!$                      l(:), v_disp, radius_calc, halo_mass, imass
!!$                 halo_pos(1)=n_box-1
!!$              end if
!!$              if(halo_pos(2).ge.n_box)then
!!$                 print*,'out of bounds',halo_pos(:), x_mean(:), v_mean(:), &
!!$                      l(:), v_disp, radius_calc, halo_mass, imass
!!$                 halo_pos(2)=n_box-1
!!$              end if
!!$              if(halo_pos(3).ge.n_box)then
!!$                 print*,'out of bounds',halo_pos(:), x_mean(:), v_mean(:), &
!!$                      l(:), v_disp, radius_calc, halo_mass, imass
!!$                 halo_pos(3)=n_box-1
!!$              end if

              !correct halo positions for the random offsets
              !halo_pos = halo_pos-offset

              !re-map negative numbers periodically
              halo_pos=modulo(halo_pos,real(n_box))
              
              ! Positions on different resolution grids
              ! need real division in case the two grids do not divide neatly 
!              p4(:)=modulo(int(halo_pos(:)/(n_box/real(n_coarse))),n_coarse)+1
              p8(:)=modulo(int(halo_pos(:)/(n_box/real(n_coarser))),n_coarser)+1
              p16(:)=modulo(int(halo_pos(:)/(n_box/real(n_coarsest))),n_coarsest)+1
!              p4(:)=modulo(int(halo_pos(:)/(n_box/n_coarse)),n_coarse)+1
!              p8(:)=modulo(int(halo_pos(:)/(n_box/n_coarser)),n_coarser)+1
!              p16(:)=modulo(int(halo_pos(:)/(n_box/n_coarsest)),n_coarsest)+1 
              !              print*, p4,p8,p16
              ! fill array with halo numbers
!              n1halo(p4(1),p4(2),p4(3),halo_type) = &
!                   n1halo(p4(1),p4(2),p4(3),halo_type)+1
              n2halo(p8(1),p8(2),p8(3),halo_type) = &
                   n2halo(p8(1),p8(2),p8(3),halo_type)+1
              n3halo(p16(1),p16(2),p16(3),halo_type) = &
                   n3halo(p16(1),p16(2),p16(3),halo_type)+1
              
              ! add masses
!              m1halo(p4(1),p4(2),p4(3),halo_type)= &
!                   m1halo(p4(1),p4(2),p4(3),halo_type)+ &
!                   halo_mass
              m2halo(p8(1),p8(2),p8(3),halo_type)= &
                   m2halo(p8(1),p8(2),p8(3),halo_type)+ &
                   halo_mass
              m3halo(p16(1),p16(2),p16(3),halo_type)= &
                   m3halo(p16(1),p16(2),p16(3),halo_type)+ &
                   halo_mass
              
 !          end if
        end if

     end do
133  continue

     coll_frac=coll_frac/M_box
     print*,'coll. fracs =',coll_frac

     close(31)

     ! Construct the output files

     !do igrid=1,3 ! loop over three different resolutions
     do igrid=2,3 ! loop over two different resolutions
        num_sources=0
        select case (igrid)
!        case (1)
!           id="coarse"
!           ngrid = n_coarse
!           nhalo => n1halo
!           mhalo => m1halo
        case (2) 
           id="coarser"
           ngrid = n_coarser
           nhalo => n2halo
           mhalo => m2halo
        case (3)
           id="coarsest"
           ngrid = n_coarsest
           nhalo => n3halo
           mhalo => m3halo
        end select

        ! count sources
        do k=1,ngrid
           do j=1,ngrid
              do i=1,ngrid
                 if (nhalo(i,j,k,1).ne.0.or.nhalo(i,j,k,2).ne.0) &
                      num_sources=num_sources+1
              enddo
           enddo
        enddo

        ! Write source list
        mfile=trim(adjustl(z_str))//"-"//trim(adjustl(id))//"_sources.dat"
        open(unit=80,file=mfile,status='unknown')
        write(80,*) num_sources
        do k=1,ngrid
           do j=1,ngrid
              do i=1,ngrid
                 if (nhalo(i,j,k,1).ne.0.or.nhalo(i,j,k,2).ne.0) & 
                      write(80,*) i,j,k,mhalo(i,j,k,2),mhalo(i,j,k,1)
              enddo
           enddo
        enddo
        close(80)
        num_sources_arr(igrid)=num_sources
     enddo

     ! write how many sources there are vs. z
     !write(104,*) z(nz),(num_sources_arr(i),i=1,3),(real(coll_frac(i)),i=1,2)  
     !write(*,*) z(nz),(num_sources_arr(i),i=1,3),(real(coll_frac(i)),i=1,2)  
     write(104,*) z(nz),(num_sources_arr(i),i=2,3),(real(coll_frac(i)),i=1,2)  
     write(*,*) z(nz),(num_sources_arr(i),i=2,3),(real(coll_frac(i)),i=1,2)  
     !write(104,*) z,(num_sources_arr(i),i=2,3) 

  end do redshift

  ! 
  close(unit=8) ! redshift file
  close(unit=104)

end program Make_Source_List
