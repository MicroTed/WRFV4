
!IDEAL:MODEL_LAYER:INITIALIZATION
!

!  This MODULE holds the routines which are used to perform various initializations
!  for the individual domains.  

!  This MODULE CONTAINS the following routines:

!  initialize_field_test - 1. Set different fields to different constant
!                             values.  This is only a test.  If the correct
!                             domain is not found (based upon the "id")
!                             then a fatal error is issued.               

!-----------------------------------------------------------------------

MODULE module_initialize_ideal

   USE module_domain
   USE module_io_domain
   USE module_state_description
   USE module_model_constants
   USE module_bc
   USE module_timing
   USE module_configure
   USE module_init_utilities
   USE module_soil_pre
#ifdef DM_PARALLEL
   USE module_dm
#endif

  private rndnum

  real, private :: ugrid = 0.0        ! east-west grid motion (subtracted from U wind component)
  real, private :: vgrid = 0.0        ! north-south grid motion (subtracted from V wind component)
  logical :: use_snd_plevels = .false. ! whether to read in pressure from the sounding file and use the same levels
  logical :: use_snd_zlevels = .false. ! whether to use Z levels the sounding file and create pressure

 CONTAINS


!-------------------------------------------------------------------
! this is a wrapper for the solver-specific init_domain routines.
! Also dereferences the grid variables and passes them down as arguments.
! This is crucial, since the lower level routines may do message passing
! and this will get fouled up on machines that insist on passing down
! copies of assumed-shape arrays (by passing down as arguments, the
! data are treated as assumed-size -- ie. f77 -- arrays and the copying
! business is avoided).  Fie on the F90 designers.  Fie and a pox.
! NOTE:  Modified to remove all but arrays of rank 4 or more from the
!        argument list.  Arrays with rank>3 are still problematic due to the
!        above-noted fie- and pox-ities.  TBH 20061129.  

   SUBROUTINE init_domain ( grid )

   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid
   !  Local data.
   INTEGER :: idum1, idum2

   
   CALL set_scalar_indices_from_config ( head_grid%id , idum1, idum2 )

     CALL init_domain_rk( grid &
!
#include "actual_new_args.inc"
!
                        )

   END SUBROUTINE init_domain

!-------------------------------------------------------------------

   SUBROUTINE init_domain_rk ( grid &
!
# include "dummy_new_args.inc"
!
)
   USE module_optional_input
   IMPLICIT NONE

   !  Input data.
   TYPE (domain), POINTER :: grid

# include "dummy_new_decl.inc"

   TYPE (grid_config_rec_type)              :: config_flags

   !  Local data
   INTEGER                             ::                       &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte, &
                                  i, j, k, kk, n

   ! Local data

   INTEGER, PARAMETER :: nl_max = 1000
   REAL, DIMENSION(nl_max) :: zk, p_in, theta, rho, u, v, qv, pd_in, etal
   REAL, DIMENSION(nl_max) :: alttde
   REAL    :: tmp_v_1, tmp_v_2, tmp_1, tmp_2
   INTEGER :: nl_in, icount


   INTEGER :: icm,jcm, ii, im1, jj, jm1, loop, error, fid
   REAL    :: u_mean,v_mean, f0, p_surf, p_level, qvf, z_at_v, z_at_u
   REAL    :: z_scale, xrad, yrad, zrad, rad, cof1, cof2
!   REAL, EXTERNAL :: interp_0
   REAL    :: hm, xa
   REAL    :: pi, rnd

!  stuff from original initialization that has been dropped from the Registry
   REAL    :: vnu, xnu, xnus, dinit0, cbh, p0_temp, t0_temp, zd, zt
   REAL    :: qvf1, qvf2, pd_surf, theta_surf
   INTEGER :: it
   real :: thtmp, ptmp, tem, temp(3)
   real :: t_min, t_max, xpos, xposml, xpospl ! used by grav2d_x

   LOGICAL :: moisture_init
   LOGICAL :: stretch_grid, dry_sounding
   character (len=256) :: mminlu2

   REAL    :: xa1, xal1,pii,hm1  !  data for intercomparison setup from dale
!  space for initial jet in b_wave

   INTEGER, parameter :: nz_jet=64, ny_jet=80
   REAL, DIMENSION(nz_jet, ny_jet) :: u_jet, rho_jet, th_jet, z_jet

!  perturbation parameters for b_wave

   REAL, PARAMETER :: htbub=8000., radbub=2000000., radz=8000., tpbub=1.0
   REAL :: piov2, tp
   INTEGER :: icen, jcen

   REAL    :: B1, B2, B3, B4, B5, sin_arg

   REAL    :: Nsq, z, z1, z2
   INTEGER :: iter_loop
   INTEGER :: xs , xe , ys , ye
   REAL :: mtn_ht = 500.
#if defined ( MTN )
   integer :: ido_terrain = 1 ! 1 = mountain, 2 = EW Ridge, 3 = NS Ridge
#elseif defined ( EW_RIDGE )
   integer :: ido_terrain = 2 ! 1 = mountain, 2 = EW Ridge, 3 = NS Ridge
#elseif defined ( NS_RIDGE )
   integer :: ido_terrain = 3 ! 1 = mountain, 2 = EW Ridge, 3 = NS Ridge
#else
   integer :: ido_terrain = 0 ! 1 = mountain, 2 = EW Ridge, 3 = NS Ridge
#endif
   REAL :: randx       ! for LES
   INTEGER :: ks, ke, id

   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   
! soil input
   INTEGER :: ns_input
   REAL    :: tmn_input, tsk_input
   REAL    :: zs_input(100),tslb_input(100),smois_input(100)
   
!    CALL nl_set_cen_lat(1,40.)
!    CALL nl_set_cen_lon(1,-105.)

    real :: centlat = 40.0, centlon = -105.0
    real :: delt = 0.0
   integer :: lu_index = 2
   real :: vegfra = 0.5, canwat = 0.0
   integer :: isltyp = 4
    
   INTEGER :: nxc=-1, nyc=-1
   real    :: pres, qvs, rh, qvs1
   real    :: zradbub = 1500.
   real    :: xradbub = 10000.
   real    :: yradbub = 10000.
   logical :: use_eta_levels = .false.
   integer, parameter :: maxbub = 10
   integer :: numbub = 1
   integer :: nxb(maxbub) = -1 , nyb(maxbub) = -1

   logical :: setccn = .false.
   
   integer :: istat
    
   logical :: moisten = .false.
   logical :: do_rebalance = .true.
   
   integer :: lm = 0
   
   integer :: lu_input = 0
   
   real :: sst_input = 300.
   LOGICAL :: real_soil = .false.

   integer            :: ibbleseed = -1
   integer :: bubbletype = 1 ! =2 to add random perts
   real, allocatable :: ranarray(:,:),ranarraysmth(:,:)
!   real :: rndnum
  integer           :: iranseed, i1, j1
  double precision  :: dthtot, dthrantot, dthrantotv, dthran
  double precision  :: dthrantotn, dthrantotp
  double precision  :: dthrantotvn, dthrantotvp

   NAMELIST /init_ideal/ nxc,nyc,delt,use_eta_levels, &
        ugrid,vgrid,setccn,use_snd_plevels,use_snd_zlevels,numbub,nxb,nyb, &
        xradbub,yradbub,zradbub,moisten,mtn_ht,ido_terrain,centlat, centlon, lm, &
        lu_input, vegfra, canwat, isltyp, sst_input, bubbletype, ibbleseed, &
        real_soil,do_rebalance

   SELECT CASE ( model_data_order )
         CASE ( DATA_ORDER_ZXY )
   kds = grid%sd31 ; kde = grid%ed31 ;
   ids = grid%sd32 ; ide = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   kms = grid%sm31 ; kme = grid%em31 ;
   ims = grid%sm32 ; ime = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   kts = grid%sp31 ; kte = grid%ep31 ;   ! note that tile is entire patch
   its = grid%sp32 ; ite = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XYZ )
   ids = grid%sd31 ; ide = grid%ed31 ;
   jds = grid%sd32 ; jde = grid%ed32 ;
   kds = grid%sd33 ; kde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   jms = grid%sm32 ; jme = grid%em32 ;
   kms = grid%sm33 ; kme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   jts = grid%sp32 ; jte = grid%ep32 ;   ! note that tile is entire patch
   kts = grid%sp33 ; kte = grid%ep33 ;   ! note that tile is entire patch
         CASE ( DATA_ORDER_XZY )
   ids = grid%sd31 ; ide = grid%ed31 ;
   kds = grid%sd32 ; kde = grid%ed32 ;
   jds = grid%sd33 ; jde = grid%ed33 ;

   ims = grid%sm31 ; ime = grid%em31 ;
   kms = grid%sm32 ; kme = grid%em32 ;
   jms = grid%sm33 ; jme = grid%em33 ;

   its = grid%sp31 ; ite = grid%ep31 ;   ! note that tile is entire patch
   kts = grid%sp32 ; kte = grid%ep32 ;   ! note that tile is entire patch
   jts = grid%sp33 ; jte = grid%ep33 ;   ! note that tile is entire patch

   END SELECT

  open(15,file="namelist.input",status='old',form='formatted')
  rewind(15)
  read(15,NML=init_ideal,iostat=istat) ! returning iostat prevents bombing if the namelist does not exist
  close(15)

   CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )

!-----------------------------------------------------------------------------
! TEMPERATURE perturbations

    dthtot      = 0.0d0
    dthrantot   = 0.0d0
    dthrantotv  = 0.0d0
    dthrantotn  = 0.0d0
    dthrantotp  = 0.0d0
    dthrantotvn = 0.0d0
    dthrantotvp = 0.0d0
   IF ( bubbletype > 1 ) THEN
   !  write(0,*) 'allocate ranarray: ids:ide,jds:jde = ',ids,ide,jds,jde
     allocate ( ranarray(ids:ide,jds:jde) )
     allocate ( ranarraysmth(ids:ide,jds:jde) )

    iranseed = -Abs(ibbleseed)
!    dfac = 0.0
    ranarray(:,:) = 0.0
    ranarraysmth(:,:) = 0.0

    DO j = jds,jde
     DO i = ids,ide
       ranarray(i,j) = 2.0*(rndnum(iranseed) - 0.5)  ! values of -1 to 1
     ENDDO
    ENDDO

    DO j = jds,jde-4,4
     DO i = ids,ide-4,4
       DO  j1 = 0,4
       DO  i1 = 0,4
         ranarraysmth(i+i1,j+j1) = ranarray(i,j)
       ENDDO
       ENDDO
     ENDDO
    ENDDO
     
   ENDIF
   
  ideal_constants: SELECT CASE ( model_config_rec%ideal_case )
  CASE ( hill2d_x )
   hm = 100.
   xa = 5.0

   icm = ide/2


   xa1  = 5000./500.
   xal1 = 4000./500.
   pii  = 2.*asin(1.0)
   hm1  = 250.
!   hm1  = 1000.


   stretch_grid = .true.
!   stretch_grid = .not. ( use_eta_levels .or. use_snd_plevels .or. use_snd_zlevels ) ! .true.
   delt = 0.
!   z_scale = .40
   z_scale = 8000./config_flags%ztop
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2

  CASE ( quarter_ss, squall2d_x, squall2d_y )

!   stretch_grid = .true.
   stretch_grid = .not. ( use_eta_levels .or. use_snd_plevels .or. use_snd_zlevels ) ! .true.
   IF (delt == 0.0 ) THEN
     delt = 3.
   ENDIF
!   z_scale = .40
   z_scale = 8000./config_flags%ztop
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   IF ( nxc <= 0 .and. nxb(1) <= 0) THEN

      IF ( nxc == -2 ) THEN ! put bubble in center of patch
      nxc = (ime+ims)/2
      write(6,*) 'nxc = ',nxc
      nxb(1) = nxc
      ELSE
      nxc = (ide-ids)/3 + 1
      nxb(1) = nxc
      ENDIF
      
   ELSEIF  ( nxc > 0 .and. nxb(1) <= 0) THEN
     nxb(1) = nxc
   ENDIF
   IF ( nyc <= 0 .and. nyb(1) <= 0 ) THEN
      IF ( nyc == -2 ) THEN ! put bubble in center of patch
      nyc = (jme+jms)/2
      write(6,*) 'nyc = ',nyc
      nyb(1) = nyc
      ELSE
      nyc = (jde-jds)/2 + 1
      nyb(1) = nyc
      ENDIF
   ELSEIF ( nyc > 0 .and. nyb(1) <= 0 ) THEN
      nyb(1) = nyc
   ENDIF


  CASE (grav2d_x)

   stretch_grid = .true.
!    z_scale = 1.675
   z_scale = 8000./config_flags%ztop
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2

  CASE (convrad)
   stretch_grid = .true.
   delt = 1.
!   z_scale = .32
   z_scale = 8000./config_flags%ztop
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = jde/2
   icm = ide/2
! lm is the half width of the land in terms of grid points
   !IF ( lm == 0 ) lm = 25
   !write(6,*) 'lm,icm-lm,icm+lm = ', lm,icm-lm,icm+lm

  CASE (b_wave)

! Note: htbub, radbub, radz, tpbub set above as parameters

   piov2 = 2.*atan(1.0)
   icen = ide/4
   jcen = jde/2

   stretch_grid = .true.
   delt = 0.
!   z_scale = .50
   z_scale = 8000./config_flags%ztop
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/4
   nyc = (jde-jds)/2

  CASE (seabreeze2d_x)

   stretch_grid = .true.
   delt = 6.
!   z_scale = .40
   z_scale = 8000./config_flags%ztop
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = jde/2
   icm = ide/2
! lm is the half width of the land in terms of grid points
   IF ( lm == 0 ) lm = 25
   write(6,*) 'lm,icm-lm,icm+lm = ', lm,icm-lm,icm+lm

  CASE (les)
  !  FOR LES, set stretch to false
   stretch_grid = .false.
   delt = 3.
!   z_scale = .40
   z_scale = 8000./config_flags%ztop  ! not used if stretch_grid = .false.
   pi = 2.*asin(1.0)
   write(6,*) ' pi is ',pi
   nxc = (ide-ids)/2
   nyc = (jde-jds)/2


  CASE DEFAULT

      WRITE( wrf_err_message , * ) 'Need to choose valid non-zero ideal_case:  ideal_case = ', model_config_rec%ideal_case
      CALL wrf_error_fatal ( wrf_err_message )

  END SELECT ideal_constants




! here we check to see if the boundary conditions are set properly

   CALL boundary_condition_check( config_flags, bdyzone, error, grid%id )

   moisture_init = .true.

    grid%itimestep=0

#ifdef DM_PARALLEL
   CALL wrf_dm_bcast_bytes( icm , IWORDSIZE )
   CALL wrf_dm_bcast_bytes( jcm , IWORDSIZE )
#endif

  ideal_landmap: SELECT CASE ( model_config_rec%ideal_case )

   CASE(quarter_ss)

    IF ( model_config_rec%sf_surface_physics(grid%id) == 0 ) THEN
    
    CALL nl_set_mminlu(1,'    ')
    CALL nl_set_iswater(1,0)
    CALL nl_set_cen_lat(1,centlat)
    CALL nl_set_cen_lon(1,centlon)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)
    
    ELSE
    
    mminlu2 = ' '
    mminlu2(1:4) = 'USGS'
    CALL nl_set_mminlu(1, mminlu2)
!   CALL nl_set_mminlu(1, 'USGS')
    CALL nl_set_iswater(1,16)
    CALL nl_set_isice(1,3)
!    CALL nl_set_cen_lat(1,20.)
!    CALL nl_set_cen_lon(1,-105.)
    CALL nl_set_cen_lat(1,centlat)
    CALL nl_set_cen_lon(1,centlon)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)
!   CALL model_to_grid_config_rec(1,model_config_rec,config_flags)
    CALL nl_get_iswater(1,grid%iswater)
    
    ENDIF

!    CALL nl_set_mminlu(1,'    ')
!    CALL nl_set_iswater(1,0)
!    CALL nl_set_truelat1(1,0.)
!    CALL nl_set_truelat2(1,0.)
!    CALL nl_set_moad_cen_lat (1,0.)
!    CALL nl_set_stand_lon (1,0.)
!    CALL nl_set_pole_lon (1,0.)
!    CALL nl_set_pole_lat (1,90.)
!    CALL nl_set_map_proj(1,0)


!  here we initialize data that currently is not initialized
!  in the input data

    IF ( lu_input == 0 ) lu_input = 7

    DO j = jts, jte
      DO i = its, ite
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%msfvy(i,j)    = 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.

         grid%e(i,j)        = 0.
         grid%f(i,j)        = 0.
         grid%xlat(i,j)     = centlat
         grid%xlong(i,j)     = centlon

       IF ( model_config_rec%sf_surface_physics(grid%id) /= 0 ) THEN
         grid%xland(i,j)     = 1
         grid%lu_index(i,j)  = lu_input
         grid%tsk(i,j) = 280.0
         grid%tmn(i,j) = 280.0
       ENDIF

      END DO
   END DO

! for Noah LSM, additional variables need to be initialized

   IF ( model_config_rec%sf_surface_physics(grid%id) /= 0 ) THEN
   quarter_masked_fields : SELECT CASE ( model_config_rec%sf_surface_physics(grid%id) )

      CASE (SLABSCHEME)

      CASE (LSMSCHEME)

        DO j = jts , MIN(jde-1,jte)
           DO i = its , MIN(ide-1,ite)
              IF (grid%xland(i,j) .lt. 1.5) THEN
                 grid%vegfra(i,j) =  vegfra*100. ! 50.0
                 grid%canwat(i,j) = canwat
                 grid%ivgtyp(i,j) = 18
                 grid%isltyp(i,j) = isltyp ! 8
                 grid%xice(i,j) = 0.
                 grid%snow(i,j) = 0.
              ELSE
                 grid%vegfra(i,j) = 0.
                 grid%canwat(i,j) = 0.
                 grid%ivgtyp(i,j) = 16
                 grid%isltyp(i,j) = 14
                 grid%xice(i,j) = 0.
                 grid%snow(i,j) = 0.
              ENDIF
           END DO
        END DO

      CASE (RUCLSMSCHEME)

   END SELECT quarter_masked_fields

! Process the soil; note that there are some things hard-wired into share/module_soil_pre.F
      CALL process_soil_ideal(grid%xland,grid%xice,grid%vegfra,grid%snow,grid%canwat, &
                     grid%ivgtyp,grid%isltyp,grid%tslb,grid%smois, &
                     grid%tsk,grid%tmn,grid%zs,grid%dzs,model_config_rec%num_soil_layers, &
                     model_config_rec%sf_surface_physics(grid%id), &
                                   ids,ide, jds,jde, kds,kde,&
                                   ims,ime, jms,jme, kms,kme,&
                                   its,ite, jts,jte, kts,kte )

    ENDIF
    
   CASE(hill2d_x, squall2d_x, squall2d_y, grav2d_x)


    CALL nl_set_mminlu(1,'    ')
    CALL nl_set_iswater(1,0)
    CALL nl_set_cen_lat(1,centlat)
    CALL nl_set_cen_lon(1,centlon)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)


!  here we initialize data that currently is not initialized
!  in the input data

    DO j = jts, jte
      DO i = its, ite
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%msfvy(i,j)    = 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.

         grid%e(i,j)        = 0.
         grid%f(i,j)        = 0.

      END DO
   END DO

   CASE(b_wave)

    CALL nl_set_mminlu(1,'    ')
    CALL nl_set_iswater(1,0)
    CALL nl_set_cen_lat(1,centlat)
    CALL nl_set_cen_lon(1,centlon)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)


!  here we initialize data that currently is not initialized
!  in the input data

    DO j = jts, jte
      DO i = its, ite
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%msfvy(i,j)    = 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
         grid%f(i,j)        = 1.e-04

      END DO
   END DO

   CASE(convrad)
    mminlu2 = ' '
    mminlu2(1:4) = 'USGS'
    CALL nl_set_mminlu(1, mminlu2)
!   CALL nl_set_mminlu(1, 'USGS')
    CALL nl_set_iswater(1,16)
    CALL nl_set_isice(1,3)
    CALL nl_set_cen_lat(1,20.)
    CALL nl_set_cen_lon(1,-105.)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)
!   CALL model_to_grid_config_rec(1,model_config_rec,config_flags)
    CALL nl_get_iswater(1,grid%iswater)

!  here we initialize data that currently is not initialized
!  in the input data

    IF ( lu_input == 0 ) lu_input = 18
    
    DO j = jts, jte
      DO i = its, ite
         grid%msft(i,j)     = 1.
         grid%msfu(i,j)     = 1.
         grid%msfv(i,j)     = 1.
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvy(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
         grid%xlat(i,j)     = 10.
         grid%f(i,j)        = 2.5e-5
         grid%xlong(i,j)     = 0.
! Hard-wire the ocean-land configuration
! all ocean
       if ( lm > 0 .and. i .ge. (icm-lm) .and. i .lt. (icm+lm)) then ! land
         grid%xland(i,j)     = 1.
         grid%lu_index(i,j)  = lu_input
         grid%landmask(i,j)  = 1.
       else ! sea
         grid%xland(i,j)     = 2.
         grid%lu_index(i,j)  = 16
         grid%sst(i,j) = sst_input
       end if
      END DO
   END DO

! for Noah LSM, additional variables need to be initialized

   convrad_masked_fields : SELECT CASE ( model_config_rec%sf_surface_physics(grid%id) )

      CASE (SLABSCHEME)

      CASE (LSMSCHEME)

        
        DO j = jts , MIN(jde-1,jte)
           DO i = its , MIN(ide-1,ite)
              IF (grid%xland(i,j) .lt. 1.5) THEN ! land
                 grid%vegfra(i,j) = vegfra*100. ! 50.
                 grid%canwat(i,j) = 0.
                 grid%ivgtyp(i,j) = 18
                 grid%isltyp(i,j) = 8
                 grid%xice(i,j) = 0.
                 grid%snow(i,j) = 0.
              ELSE ! sea
                 grid%vegfra(i,j) = 0.
                 grid%canwat(i,j) = 0.
                 grid%ivgtyp(i,j) = 16
                 grid%isltyp(i,j) = 14
                 grid%xice(i,j) = 0.
                 grid%snow(i,j) = 0.
                 grid%smois(i,:,j) = 1.
                 grid%tslb(i,:,j) = 273.16
              ENDIF
           END DO
        END DO

      CASE (RUCLSMSCHEME)

   END SELECT convrad_masked_fields
! Process the soil; note that there are some things hard-wired into share/module_soil_pre.F
      
    IF ( real_soil ) THEN ! from input file

      write(0,*) 'read input_soil'
      
      IF (config_flags%sf_surface_physics .NE. 2) WRITE (6, *)   &
         'If using LSM option other than Noah, must edit input_soil file in test/em_scm_xy/ directory'
    
      CALL read_soil(100,ns_input,tmn_input,tsk_input,zs_input,tslb_input,smois_input)

      CALL init_module_optional_input(grid,config_flags)
      num_st_levels_input = ns_input
      num_sm_levels_input = ns_input
      num_sw_levels_input = ns_input
      DO k = 1,ns_input
         st_levels_input(k) = zs_input(k)*100.0 ! to cm
         sm_levels_input(k) = zs_input(k)*100.0 ! to cm
         sw_levels_input(k) = zs_input(k)*100.0 ! to cm
         st_input(:,k+1,:) = tslb_input(k)
         sm_input(:,k+1,:) = smois_input(k)
         sw_input(:,k+1,:) = smois_input(k)
      ENDDO
  
      grid%tsk = tsk_input
      grid%sst = tsk_input
      grid%tmn = tmn_input

      flag_soil_layers  = 0 ! go ahead and put skin temp in
      flag_soil_levels  = 0 ! go ahead and put skin moisture in
      flag_sst          = 0 ! don't modify for ocean
      flag_tavgsfc      = 0 
      flag_soilhgt      = 0 

      CALL process_soil_real ( grid%tsk , grid%tmn , grid%tavgsfc, &
                   grid%landmask , grid%sst , grid%ht, grid%toposoil, &
                   st_input , sm_input , sw_input , &
                   st_levels_input , sm_levels_input , sw_levels_input , &
                   grid%zs , grid%dzs , grid%flag_sm_adj , grid%tslb , grid%smois , grid%sh2o , &
                   flag_sst , flag_tavgsfc, flag_soilhgt, flag_soil_layers, flag_soil_levels,  &
                   ids , ide , jds , jde , kds , kde , &
                   ims , ime , jms , jme , kms , kme , &
                   its , ite , jts , jte , kts , kte , &
                   model_config_rec%sf_surface_physics(grid%id) , &
                   model_config_rec%num_soil_layers , &
                   model_config_rec%real_data_init_type , &
                   num_st_levels_input , num_sm_levels_input , num_sw_levels_input , &
                   num_st_levels_alloc , num_sm_levels_alloc , num_sw_levels_alloc )
      
    ELSE
    
      CALL process_soil_ideal(grid%xland,grid%xice,grid%vegfra,grid%snow,grid%canwat, &
                     grid%ivgtyp,grid%isltyp,grid%tslb,grid%smois, &
                     grid%tsk,grid%tmn,grid%zs,grid%dzs,model_config_rec%num_soil_layers, &
                     model_config_rec%sf_surface_physics(grid%id), &
                                   ids,ide, jds,jde, kds,kde,&
                                   ims,ime, jms,jme, kms,kme,&
                                   its,ite, jts,jte, kts,kte )
     ENDIF
   CASE(seabreeze2d_x)

    mminlu2 = ' '
    mminlu2(1:4) = 'USGS'
    CALL nl_set_mminlu(1, mminlu2)
!   CALL nl_set_mminlu(1, 'USGS')
    CALL nl_set_iswater(1,16)
    CALL nl_set_isice(1,3)
    CALL nl_set_cen_lat(1,20.)
    CALL nl_set_cen_lon(1,-105.)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)
!   CALL model_to_grid_config_rec(1,model_config_rec,config_flags)
    CALL nl_get_iswater(1,grid%iswater)

!  here we initialize data that currently is not initialized
!  in the input data

    DO j = jts, jte
      DO i = its, ite
         grid%msft(i,j)     = 1.
         grid%msfu(i,j)     = 1.
         grid%msfv(i,j)     = 1.
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvy(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
         grid%f(i,j)        = 0.
         grid%xlat(i,j)     = centlat
         grid%xlong(i,j)     = centlon
! Hard-wire the ocean-land configuration
        if (i .ge. (icm-lm) .and. i .lt. (icm+lm)) then
         grid%xland(i,j)     = 1.
         grid%lu_index(i,j)  = 18
         grid%tsk(i,j) = 280.0
         grid%tmn(i,j) = 280.0
        else
         grid%xland(i,j)     = 2.
         grid%lu_index(i,j)  = 16
         grid%tsk(i,j) = 287.0
         grid%tmn(i,j) = 280.0
        end if
      END DO
   END DO

! for Noah LSM, additional variables need to be initialized

   seabreeze_masked_fields : SELECT CASE ( model_config_rec%sf_surface_physics(grid%id) )

      CASE (SLABSCHEME)

      CASE (LSMSCHEME)

        DO j = jts , MIN(jde-1,jte)
           DO i = its , MIN(ide-1,ite)
              IF (grid%xland(i,j) .lt. 1.5) THEN
                 grid%vegfra(i,j) = vegfra*100. ! 50.
                 grid%canwat(i,j) = 0.
                 grid%ivgtyp(i,j) = 18
                 grid%isltyp(i,j) = 8
                 grid%xice(i,j) = 0.
                 grid%snow(i,j) = 0.
              ELSE
                 grid%vegfra(i,j) = 0.
                 grid%canwat(i,j) = 0.
                 grid%ivgtyp(i,j) = 16
                 grid%isltyp(i,j) = 14
                 grid%xice(i,j) = 0.
                 grid%snow(i,j) = 0.
              ENDIF
           END DO
        END DO

      CASE (RUCLSMSCHEME)

   END SELECT seabreeze_masked_fields

! Process the soil; note that there are some things hard-wired into share/module_soil_pre.F
      CALL process_soil_ideal(grid%xland,grid%xice,grid%vegfra,grid%snow,grid%canwat, &
                     grid%ivgtyp,grid%isltyp,grid%tslb,grid%smois, &
                     grid%tsk,grid%tmn,grid%zs,grid%dzs,model_config_rec%num_soil_layers, &
                     model_config_rec%sf_surface_physics(grid%id), &
                                   ids,ide, jds,jde, kds,kde,&
                                   ims,ime, jms,jme, kms,kme,&
                                   its,ite, jts,jte, kts,kte )

   CASE (les)

    CALL nl_set_mminlu(1, '    ')
    CALL nl_set_iswater(1,0)
    CALL nl_set_cen_lat(1,40.)
    CALL nl_set_cen_lon(1,-105.)
    CALL nl_set_truelat1(1,0.)
    CALL nl_set_truelat2(1,0.)
    CALL nl_set_moad_cen_lat (1,0.)
    CALL nl_set_stand_lon (1,0.)
    CALL nl_set_pole_lon (1,0.)
    CALL nl_set_pole_lat (1,90.)
    CALL nl_set_map_proj(1,0)


!  here we initialize data that currently is not initialized
!  in the input data

    DO j = jts, jte
      DO i = its, ite
         grid%msftx(i,j)    = 1.
         grid%msfty(i,j)    = 1.
         grid%msfux(i,j)    = 1.
         grid%msfuy(i,j)    = 1.
         grid%msfvx(i,j)    = 1.
         grid%msfvx_inv(i,j)= 1.
         grid%msfvy(i,j)    = 1.
         grid%sina(i,j)     = 0.
         grid%cosa(i,j)     = 1.
         grid%e(i,j)        = 0.
!  for LES, include Coriolis force
         grid%f(i,j)        = 1.e-4

      END DO
   END DO

   END SELECT ideal_landmap

    DO j = jts, jte
    DO k = kts, kte
      DO i = its, ite
         grid%ww(i,k,j)     = 0.
      END DO
   END DO
   END DO

   grid%step_number = 0

! set up the grid

   ideal_levels: SELECT CASE ( model_config_rec%ideal_case )
   CASE(hill2d_x, quarter_ss, squall2d_x, squall2d_y, grav2d_x, b_wave)
   IF (stretch_grid) THEN ! exponential stretch for eta (nearly constant dz)
     DO k=1, kde
      grid%znw(k) = (exp(-(k-1)/float(kde-1)/z_scale) - exp(-1./z_scale))/ &
                                (1.-exp(-1./z_scale))
     ENDDO
   ELSE

    IF ( use_snd_plevels .or. use_snd_zlevels ) THEN
    CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, etal, dry_sounding, nl_max, nl_in, theta_surf)
     grid%znw(1) = 1.0
     write(6,*) 'eta levels:'
     write(6,*) 'k, eta'
     write(6,*) 1,grid%znw(1)
     DO k=2, kde-1
      grid%znw(k) = 0.5*( etal(k) + etal(k+1) )
     write(6,*) k,grid%znw(k)
     ENDDO
     k = kde
     grid%znw(k) = etal(k+1) ! use the added top w-point for model top
     write(6,*) k,grid%znw(k)
     
      write(6,"( 5(f10.7,',',1x))") (grid%znw(k), k=1,kde)
    
    ELSEIF ( use_eta_levels ) THEN

     DO k=1, kde
!      grid%znw(k) = 1. - float(k-1)/float(kde-1)
!      grid%znw(k) = eta_levels(k)
      grid%znw(k) = model_config_rec%eta_levels(k)
!      write(6,*) 'k,eta = ',k,grid%znw(k) 
     ENDDO
    ELSE
     DO k=1, kde
!      grid%znw(k) = 1. - float(k-1)/float(kde-1)
      grid%znw(k) = model_config_rec%eta_levels(k)
     ENDDO
    ENDIF
   ENDIF
   CASE(convrad, seabreeze2d_x)
   IF (stretch_grid) THEN ! read eta_levels
     DO k=1, kde
      grid%znw(k) = model_config_rec%eta_levels(k)
     ENDDO
   ELSE
     DO k=1, kde
      grid%znw(k) = 1. - float(k-1)/float(kde-1)
     ENDDO
   ENDIF
   CASE(les)
       IF (model_config_rec%eta_levels(1) .EQ. -1) THEN !we do not have eta_levels from namelist
   IF (stretch_grid) THEN ! exponential stretch for eta (nearly constant dz)
     DO k=1, kde
      grid%znw(k) = (exp(-(k-1)/float(kde-1)/z_scale) - exp(-1./z_scale))/ &
                                (1.-exp(-1./z_scale))
     ENDDO
   ELSE
     DO k=1, kde
      grid%znw(k) = 1. - float(k-1)/float(kde-1)
     ENDDO
   ENDIF
      ELSE
          CALL wrf_debug(0,"module_initialize_les: vertical nesting is enabled, using eta_levels specified in namelist.input")
          ks = 0
          DO id=1,grid%id
             ks = ks+model_config_rec%e_vert(id)
          ENDDO
          IF (ks .GT. max_eta) THEN
             CALL wrf_error_fatal("too many vertical levels, increase max_eta in frame/module_driver_constants.F")
          ENDIF
          !Now set the eta_levels to what we specified in the namelist. We've
          !packed all the domains' eta_levels into a 'vector' and now we need
          !to pull only the section of the vector associated with our domain
          !of interest, which is between indicies ks and ke.
          IF (grid%id .EQ. 1) THEN
            ks = 1
            ke = model_config_rec%e_vert(1)
          ELSE
            id = 1
            ks = 1
            ke = 0
            DO WHILE (grid%id .GT. id)
              id = id+1
              ks = ks+model_config_rec%e_vert(id-1)
              ke = ks+model_config_rec%e_vert(id)
            ENDDO
          ENDIF
          DO k=1,kde
            grid%znw(k) = model_config_rec%eta_levels(ks+k-1)
          ENDDO
          !Check the value of the first and last eta level for our domain,
          !then check that the vector of eta levels is only decreasing
          IF (grid%znw(1) .NE. 1.0) THEN
            CALL wrf_error_fatal("error with specified eta_levels, first level is not 1.0")
          ENDIF
          IF (grid%znw(kde) .NE. 0.0) THEN
            CALL wrf_error_fatal("error with specified eta_levels, last level is not 0.0")
          ENDIF
          DO k=2,kde
            IF (grid%znw(k) .GT. grid%znw(k-1)) THEN
              CALL wrf_error_fatal("eta_levels are not uniformly decreasing from 1.0 to 0.0")
            ENDIF
          ENDDO
      ENDIF
   END SELECT ideal_levels

   DO k=1, kde-1
    grid%dnw(k) = grid%znw(k+1) - grid%znw(k)
    grid%rdnw(k) = 1./grid%dnw(k)
    grid%znu(k) = 0.5*(grid%znw(k+1)+grid%znw(k))
   ENDDO
   DO k=2, kde-1
    grid%dn(k) = 0.5*(grid%dnw(k)+grid%dnw(k-1))
    grid%rdn(k) = 1./grid%dn(k)
    grid%fnp(k) = .5* grid%dnw(k  )/grid%dn(k)
    grid%fnm(k) = .5* grid%dnw(k-1)/grid%dn(k)
   ENDDO

   cof1 = (2.*grid%dn(2)+grid%dn(3))/(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(2)
   cof2 =     grid%dn(2)        /(grid%dn(2)+grid%dn(3))*grid%dnw(1)/grid%dn(3)
   grid%cf1  = grid%fnp(2) + cof1
   grid%cf2  = grid%fnm(2) - cof1 - cof2
   grid%cf3  = cof2       

   grid%cfn  = (.5*grid%dnw(kde-1)+grid%dn(kde-1))/grid%dn(kde-1)
   grid%cfn1 = -.5*grid%dnw(kde-1)/grid%dn(kde-1)
   grid%rdx = 1./config_flags%dx
   grid%rdy = 1./config_flags%dy

!  get the sounding from the ascii sounding file, first get dry sounding and
!  calculate base state
   ideal_sounding: SELECT CASE ( model_config_rec%ideal_case )

   CASE (b_wave)

  write(6,*) ' reading input jet sounding '
  call read_input_jet( u_jet, rho_jet, th_jet, z_jet, nz_jet, ny_jet )

  write(6,*) ' getting dry sounding for base state '
  write(6,*) ' using middle column in jet sounding, j = ',ny_jet/2

  dry_sounding   = .true.
!  debug = .true.  !  this will produce print of the sounding (last argument)
  CALL get_sounding_b_wave( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, &
                      nl_max, nl_in, u_jet, rho_jet, th_jet, z_jet,      &
                      nz_jet, ny_jet, ny_jet/2, .true.                   )

  write(6,*) ' returned from reading sounding, nl_in is ',nl_in

!  find ptop for the desired ztop (ztop is input from the namelist),
!  and find surface pressure

!  For the jet, using the middle column for the base state means that
!  we will be extrapolating above the highest height data to the south
!  of the centerline.


   CASE DEFAULT

  IF ( wrf_dm_on_monitor() ) THEN
  write(6,*) ' getting dry sounding for base state '
  dry_sounding = .true.
  CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, etal, dry_sounding, &
                     nl_max, nl_in, theta_surf)
  ENDIF
  CALL wrf_dm_bcast_real( zk , nl_max )
  CALL wrf_dm_bcast_real( p_in , nl_max )
  CALL wrf_dm_bcast_real( pd_in , nl_max )
  CALL wrf_dm_bcast_real( theta , nl_max )
  CALL wrf_dm_bcast_real( rho , nl_max )
  CALL wrf_dm_bcast_real( u , nl_max )
  CALL wrf_dm_bcast_real( v , nl_max )
  CALL wrf_dm_bcast_real( qv , nl_max )
  CALL wrf_dm_bcast_integer ( nl_in , 1 )

  write(6,*) ' returned from reading sounding, nl_in is ',nl_in


!  find ptop for the desired ztop (ztop is input from the namelist),
!  and find surface pressure
   END SELECT ideal_sounding

  grid%p_top = interp_0( p_in, zk, config_flags%ztop, nl_in )

  IF ( use_eta_levels ) THEN
    grid%p_top = Max( grid%p_top_requested, p_in(nl_in) )
    write(6,*) 'ptop,p_top_requested,p_in(nl_in) = ',grid%p_top,grid%p_top_requested,p_in(nl_in)
  ENDIF

!  Fill in the hybrid coordinate coefficients

   DO k=1, kde
      IF      ( config_flags%hybrid_opt .EQ. 0 ) THEN
         grid%c3f(k) = grid%znw(k)
      ELSE IF ( config_flags%hybrid_opt .EQ. 1 ) THEN
         grid%c3f(k) = grid%znw(k)
      ELSE IF ( config_flags%hybrid_opt .EQ. 2 ) THEN
         B1 = 2. * grid%etac**2 * ( 1. - grid%etac )
         B2 = -grid%etac * ( 4. - 3. * grid%etac - grid%etac**3 )
         B3 = 2. * ( 1. - grid%etac**3 )
         B4 = - ( 1. - grid%etac**2 )
         B5 = (1.-grid%etac)**4
         grid%c3f(k) = ( B1 + B2*grid%znw(k) + B3*grid%znw(k)**2 + B4*grid%znw(k)**3 ) / B5
         IF ( grid%znw(k) .LT. grid%etac ) THEN
            grid%c3f(k) = 0.
         END IF
         IF ( k .EQ. kds ) THEN
            grid%c3f(k) = 1.
         ELSE IF ( k .EQ. kde ) THEN
            grid%c3f(k) = 0.
         END IF
      ELSE IF ( config_flags%hybrid_opt .EQ. 3 ) THEN
         IF ( grid%znw(k) .GE. grid%etac ) THEN
            sin_arg = (1./(1.-grid%etac))*(grid%znw(k)-1.)+1
            grid%c3f(k) = (sin(sin_arg*3.14159265358/2.))**2
         ELSE
            grid%c3f(k) = 0.
         END IF
         IF ( k .EQ. kds ) THEN
            grid%c3f(k) = 1.
         ELSE IF ( k .EQ. kds ) THEN
            grid%c3f(kde) = 0.
         END IF
      ELSE
         CALL wrf_error_fatal ( 'ERROR: --- hybrid_opt=0 ===> Standard WRF Coordinate; hybrid_opt>=1 ===> Hybrid Vertical Coordinate' )
      END IF
   END DO

   DO k=1, kde
      grid%c4f(k) = ( grid%znw(k) - grid%c3f(k) ) * ( p1000mb - grid%p_top )
   ENDDO

   !  Now on half levels, just add up and divide by 2 (for c3h).  Use (eta-c3)*(p00-pt) for c4 on half levels.

   DO k=1, kde-1
      grid%c3h(k) = ( grid%c3f(k+1) + grid%c3f(k) ) * 0.5
      grid%c4h(k) = ( grid%znu(k) - grid%c3h(k) ) * ( p1000mb - grid%p_top )
   ENDDO

   !  c1 = d(B)/d(eta).  We define c1f as c1 on FULL levels.  For a vertical difference,
   !  we need to use B and eta on half levels.  The k-loop ends up referring to the
   !  full levels, neglecting the top and bottom.

   DO k=kds+1, kde-1
      grid%c1f(k) = ( grid%c3h(k) - grid%c3h(k-1) ) / ( grid%znu(k) - grid%znu(k-1) )
   ENDDO

   !  The boundary conditions to get the coefficients:
   !  1) At k=kts: define d(B)/d(eta) = 1.  This gives us the same value of B and d(B)/d(eta)
   !     when doing the sigma-only B=eta.
   !  2) At k=kte: with the new vertical coordinate, define d(B)/d(eta) = 0.  The curve B SMOOTHLY
   !     goes to zero, and at the very top, B continues to SMOOTHLY go to zero.  Note that for
   !     almost all cases of non B=eta, B is ALREADY=ZERO at the top, so this is a reasonable BC to
   !     assume.
   !  3) At k=kte: when trying to mimic the original vertical coordinate, since B = eta, then
   !     d(B)/d(eta) = 1.

   grid%c1f(kds) = 1.
   IF      ( ( config_flags%hybrid_opt .EQ. 0 ) .OR. ( config_flags%hybrid_opt .EQ. 1 ) ) THEN
      grid%c1f(kde) = 1.
   ELSE
      grid%c1f(kde) = 0.
   END IF

   !  c2 = ( 1. - c1(k) ) * (p00 - pt).  There is no vertical differencing, so we can do the
   !  full kds to kde looping.

   DO k=kds, kde
      grid%c2f(k) = ( 1. - grid%c1f(k) ) * ( p1000mb - grid%p_top )
   END DO

   !  Now on half levels for c1 and c2.  The c1h will result from the full level c3 and full
   !  level eta differences.  The c2 value use the half level c1(k).

   DO k=1, kde-1
      grid%c1h(k) = ( grid%c3f(k+1) - grid%c3f(k) ) / ( grid%znw(k+1) - grid%znw(k) )
      grid%c2h(k) = ( 1. - grid%c1h(k) ) * ( p1000mb - grid%p_top )
   END DO

  IF ( use_snd_plevels .or. use_snd_zlevels ) grid%p_top = p_in(nl_in)


#if 0
   DO k=1, kde
      grid%c3f(k) = grid%znw(k)
      grid%c4f(k) = 0.
      grid%c3h(k) = grid%znu(k)
      grid%c4h(k) = 0.
      grid%c1f(k) = 1.
      grid%c2f(k) = 0.
      grid%c1h(k) = 1.
      grid%c2h(k) = 0.
   END DO
#endif

  ideal_terrain: SELECT CASE ( model_config_rec%ideal_case )
  CASE (hill2d_x)
  DO j=jts,jte
  DO i=its,ite  
    grid%ht(i,j) = hm/(1.+(float(i-icm)/xa)**2)
!!    grid%ht(i,j) = hm1*exp(-(( float(i-icm)/xa1)**2))   &
!!               *( (cos(pii*float(i-icm)/xal1))**2 )
    grid%phb(i,1,j) = g*grid%ht(i,j)
    grid%php(i,1,j) = 0.
    grid%ph0(i,1,j) = grid%phb(i,1,j)
  ENDDO
  ENDDO

  CASE (quarter_ss, convrad, squall2d_x, squall2d_y, grav2d_x, b_wave, seabreeze2d_x)
  DO j=jts,jte
  DO i=its,ite  
    grid%ht(i,j) = 0.
    grid%phb(i,1,j) = g*grid%ht(i,j)
    grid%php(i,1,j) = 0.
    grid%ph0(i,1,j) = grid%phb(i,1,j)
  ENDDO
  ENDDO

  CASE (les)
  DO j=jts,jte
  DO i=its,ite
    grid%ht(i,j) = 0.
  ENDDO
  ENDDO

  xs=ide/2 -3
  xs=ids   -3
  xe=xs + 6
  ys=jde/2 -3
  ye=ys + 6
 ! mtn_ht = 500
  IF ( ido_terrain == 1 ) THEN ! MTN
  DO j=max(ys,jds),min(ye,jde-1)
  DO i=max(xs,ids),min(xe,ide-1)
     grid%ht(i,j) = mtn_ht * 0.25 * &
               ( 1. + COS ( 2*pi/(xe-xs) * ( i-xs ) + pi ) ) * &
               ( 1. + COS ( 2*pi/(ye-ys) * ( j-ys ) + pi ) )
  ENDDO
  ENDDO
  ENDIF
 
  IF ( ido_terrain == 2 ) THEN !  EW_RIDGE
  DO j=max(ys,jds),min(ye,jde-1)
  DO i=ids,ide
     grid%ht(i,j) = mtn_ht * 0.50 * &
               ( 1. + COS ( 2*pi/(ye-ys) * ( j-ys ) + pi ) )
  ENDDO
  ENDDO
  endif
  
  IF ( ido_terrain ==3 ) THEN ! NS_RIDGE
  DO j=jds,jde
  DO i=max(xs,ids),min(xe,ide-1)
     grid%ht(i,j) = mtn_ht * 0.50 * &
               ( 1. + COS ( 2*pi/(xe-xs) * ( i-xs ) + pi ) )
  ENDDO
  ENDDO
  endif
  
  DO j=jts,jte
  DO i=its,ite
    grid%phb(i,1,j) = g * grid%ht(i,j)
    grid%ph0(i,1,j) = g * grid%ht(i,j)
  ENDDO
  ENDDO

  END SELECT ideal_terrain

  DO J = jts, jte
  DO I = its, ite

    p_surf = interp_0( p_in, zk, grid%phb(i,1,j)/g, nl_in )
    grid%MUB(i,j) = p_surf-grid%p_top

!  this is dry hydrostatic sounding (base state), so given grid%p (coordinate),
!  interp theta (from interp) and compute 1/rho from eqn. of state

    DO K = 1, kte-1
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top
      grid%pb(i,k,j) = p_level
      grid%t_init(i,k,j) = interp_0( theta, p_in, p_level, nl_in ) - t0
      grid%alb(i,k,j) = (r_d/p1000mb)*(grid%t_init(i,k,j)+t0)*(grid%pb(i,k,j)/p1000mb)**cvpm
    ENDDO
!  calc hydrostatic balance (alternatively we could interp the geopotential from the
!  sounding, but this assures that the base state is in exact hydrostatic balance with
!  respect to the model eqns.
    DO kk  = 2,kte
      k=kk - 1
      grid%phb(i,kk,j) = grid%phb(i,kk-1,j) - grid%dnw(kk-1)*(grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))*grid%alb(i,kk-1,j)
    ENDDO
  ENDDO
  ENDDO
  IF ( wrf_dm_on_monitor() ) THEN
  write(6,*) ' ptop is ',grid%p_top
  write(6,*) ' base state grid%MUB(1,1), p_surf is ',grid%MUB(1,1),grid%c3f(kts)*grid%MUB(1,1)+grid%c4f(kts)+grid%p_top
  ENDIF
!  calculate full state for each column - this includes moisture.
  write(6,*) ' getting moist sounding for full state '
  IF ( model_config_rec%ideal_case .EQ. b_wave )THEN
    dry_sounding = .true.
    IF (config_flags%mp_physics /= 0)  dry_sounding = .false.
  ELSE
    dry_sounding = .false.
    CALL get_sounding( zk, p_in, pd_in, theta, rho, u, v, qv, etal, dry_sounding, &
                     nl_max, nl_in, theta_surf )
  ENDIF
  DO J = jts, min(jde-1,jte)
  
  IF ( model_config_rec%ideal_case .EQ. b_wave )THEN
!  get sounding for this point
!   debug = .false.  !  this will turn off print of the sounding
    CALL get_sounding_b_wave( zk, p_in, pd_in, theta, rho, u, v, qv, dry_sounding, &
                      nl_max, nl_in, u_jet, rho_jet, th_jet, z_jet,      &
                      nz_jet, ny_jet, j, .false.                          )
  ENDIF
  DO i = its, min(ide-1,ite)
!  At this point grid%p_top is already set. find the DRY mass in the column
!  by interpolating the DRY pressure.  
   pd_surf = interp_0( pd_in, zk, grid%phb(i,1,j)/g, nl_in )
!  compute the perturbation mass and the full mass
    grid%MU_1(i,j) = pd_surf-grid%p_top - grid%MUB(i,j)
    grid%MU_2(i,j) = grid%MU_1(i,j)
    grid%MU0(i,j) = grid%MU_1(i,j) + grid%MUB(i,j)
! given the dry pressure and coordinate system, interp the potential
! temperature and qv
    do k=1,kde-1
      p_level = grid%c3h(k)*(pd_surf - grid%p_top) + grid%c4h(k) + grid%p_top
      moist(i,k,j,P_QV) = interp_0( qv, pd_in, p_level, nl_in )
      grid%t_1(i,k,j)          = interp_0( theta, pd_in, p_level, nl_in ) - t0
      grid%t_2(i,k,j)          = grid%t_1(i,k,j)
      
      IF ( f_qnn .and. ( setccn .or. config_flags%nssl_ccn_on >= 2 ) ) THEN
       IF ( config_flags%mp_physics == wdm5scheme .or. config_flags%mp_physics == wdm6scheme ) THEN
         scalar(i,k,j,P_QNN) = grid%ccn_conc
       ELSEIF ( config_flags%mp_physics == nssl_2mom .and. config_flags%nssl_ccn_on > 0 ) THEN
!       ELSEIF ( config_flags%mp_physics == nssl_2momccn .or. config_flags%mp_physics == nssl_3mom ) THEN
         IF ( config_flags%nssl_ccn_is_ccna == 0 .or. config_flags%nssl_ccn_on >= 2 ) THEN ! qnn field is for UNactivated CCN
          ! The value of 1.225 converts number per m^3 to number per kg at sea level pressure in a perfectly mixed atmosphere
           scalar(i,k,j,P_QNN) = grid%nssl_cccn/1.225
           IF ( P_QNN_NU > 1 ) THEN
              scalar(i,k,j,P_QNN_NU) = grid%nssl_cccnnu/1.225
           ENDIF
           IF ( P_QNN_CO > 1 ) THEN
              scalar(i,k,j,P_QNN_CO) = grid%nssl_cccnco/1.225
           ENDIF
!  Alternatively one can set CCN values as a function of k
!   Values should be set as number per kilogram
!          IF ( k < 7 ) THEN
!          scalar(i,k,j,P_QNN) = grid%nssl_cccn/1.225
!          ELSE
!          scalar(i,k,j,P_QNN) = 0.1*grid%nssl_cccn/1.225
!          ENDIF
         ELSE ! qnn field is for Activated CCN (zero initial state)
           scalar(i,k,j,P_QNN) = 0.0
         ENDIF
       ELSE
         scalar(:,:,:,P_QNN) = grid%ccn_conc
       ENDIF
      ENDIF

    enddo

!  integrate the hydrostatic equation (from the RHS of the bigstep
!  vertical momentum equation) down from the top to get grid%p.
!  first from the top of the model to the top pressure
    kk = kte-1  ! top level
    k=kk+1
    qvf1 = 0.5*(moist(i,kk,j,P_QV)+moist(i,kk,j,P_QV))
    qvf2 = 1./(1.+qvf1)
    qvf1 = qvf1*qvf2
    grid%p(i,kk,j) = - 0.5*((grid%c1f(k)*grid%Mu_1(i,j))+qvf1*(grid%c1f(k)*grid%Mub(i,j)+grid%c2f(k)))/grid%rdnw(kk)/qvf2
    qvf = 1. + rvovrd*moist(i,kk,j,P_QV)
    grid%alt(i,kk,j) = (r_d/p1000mb)*(grid%t_1(i,kk,j)+t0)*qvf* &
                (((grid%p(i,kk,j)+grid%pb(i,kk,j))/p1000mb)**cvpm)
    grid%al(i,kk,j) = grid%alt(i,kk,j) - grid%alb(i,kk,j)
!  down the column
    do kk=kte-2,1,-1
      k = kk + 1
      qvf1 = 0.5*(moist(i,kk,j,P_QV)+moist(i,kk+1,j,P_QV))
      qvf2 = 1./(1.+qvf1)
      qvf1 = qvf1*qvf2
      grid%p(i,kk,j) = grid%p(i,kk+1,j) - ((grid%c1f(k)*grid%Mu_1(i,j)) + qvf1*(grid%c1f(k)*grid%Mub(i,j)+grid%c2f(k)))/qvf2/grid%rdn(kk+1)
      qvf = 1. + rvovrd*moist(i,kk,j,P_QV)
      grid%alt(i,kk,j) = (r_d/p1000mb)*(grid%t_1(i,kk,j)+t0)*qvf* &
                  (((grid%p(i,kk,j)+grid%pb(i,kk,j))/p1000mb)**cvpm)
      grid%al(i,kk,j) = grid%alt(i,kk,j) - grid%alb(i,kk,j)
    enddo
!  this is the hydrostatic equation used in the model after the
!  small timesteps.  In the model, grid%al (inverse density)
!  is computed from the geopotential.
    grid%ph_1(i,1,j) = 0.
    DO kk  = 2,kte
      k = kk-1
      grid%ph_1(i,kk,j) = grid%ph_1(i,kk-1,j) - (grid%dnw(kk-1))*(       &
                   ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,kk-1,j)+ &
                    (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,kk-1,j)  )
                                                   
      grid%ph_2(i,kk,j) = grid%ph_1(i,kk,j)
      grid%ph0(i,kk,j) = grid%ph_1(i,kk,j) + grid%phb(i,kk,j)
    ENDDO
    IF ( wrf_dm_on_monitor() ) THEN
    if((i==2) .and. (j==2)) then
     k=1
     write(6,*) ' grid%ph_1 k=1 calc ',grid%ph_1(2,k,2),&
                              (grid%c1h(k)*grid%mu_1(2,2))+(grid%c1h(k)*grid%mub(2,2)+grid%c2h(k)),(grid%c1h(k)*grid%mu_1(2,2)), &
                              grid%alb(2,k,2),grid%rdnw(k)
     k=2
     write(6,*) ' grid%ph_1 k=2 calc ',grid%ph_1(2,k,2),&
                              (grid%c1h(k)*grid%mu_1(2,2))+(grid%c1h(k)*grid%mub(2,2)+grid%c2h(k)),(grid%c1h(k)*grid%mu_1(2,2)), &
                              grid%alb(2,k,2)
    endif
    ENDIF
  IF ( model_config_rec%ideal_case .EQ. b_wave )THEN
    DO K = 1, kte
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top
      grid%u_1(i,k,j) = interp_0( u, p_in, p_level, nl_in )
      grid%u_2(i,k,j) = grid%u_1(i,k,j)
    ENDDO
  ENDIF

  ENDDO
  ENDDO

!===============
  ideal_pert: SELECT CASE ( model_config_rec%ideal_case )
  CASE (quarter_ss)
!  QSS thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

  DO J = jts, min(jde-1,jte)
    yrad = config_flags%dy*float(j-nyc)/yradbub
!   yrad = 0.
    DO I = its, min(ide-1,ite)
      xrad = config_flags%dx*float(i-nxc)/xradbub
!     xrad = 0.
      DO K = 1, kte-1

!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

               dthran = 0.0d0
           IF( bubbletype .eq. 2 ) dthran = dthran + 0.1*ranarray(i,j)

        DO n = 1,numbub
         yrad = config_flags%dy*float(j-nyb(n))/yradbub
         xrad = config_flags%dx*float(i-nxb(n))/xradbub
         
        zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
                   +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
        zrad = (zrad-zradbub)/zradbub
        RAD=SQRT(xrad*xrad+yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
           thtmp   = grid%t_1(i,k,j) + t0
           grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*COS(.5*PI*RAD)**2 + dthran
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
           IF ( moisten ) THEN

            pres  = grid%pb(i,k,j) ! psfc*pz(k)**(cp/rd)
         !   thtmp   = grid%t_1(i,k,j) - delt*COS(.5*PI*RAD)**2
            ptmp    = grid%p(i,k,j)+grid%pb(i,k,j)
            tem = thtmp * (ptmp/p1000mb)**rcp
            qvs   = 380.*exp(17.27*(tem-273.16) / (tem - 36.)) / ptmp
!            qvs   = 380.*exp(17.27*(pz(k)*tz(k)-273.16) / (pz(k)*tz(k)- 36.)) / pres
            rh    =  moist(i,k,j,P_QV) / qvs

            thtmp   = grid%t_1(i,k,j) + t0 ! this is the warmed theta
            tem = thtmp * (ptmp/p1000mb)**rcp
            
            qvs1  =  380.*exp(17.27*(tem-273.16) / (tem - 36.)) / ptmp
            
          !  write(0,*) 'i,j,k,pres,thtmp,tem,qvs,rh,qvs1 = ',i,j,k,pres,thtmp,tem,qvs,rh,qvs1
            
            moist(i,k,j,P_QV) = rh*qvs1
!            qv3(i,j,k) = rh*qvs1
             
           ENDIF
        ENDIF
        ENDDO
       ENDDO

!  rebalance hydrostatically

      IF ( do_rebalance ) THEN
      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )

        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO
      ENDIF

    ENDDO
  ENDDO

   IF ( allocated( ranarray ) ) deallocate( ranarray, ranarraysmth )

! ub>>
! Use the hypsometric equation to calculate altitudes (m) at the chosen
! eta (sigma-ish) coordinates.

   alttde(1) = 0.0
   do k = 2, kde-1
     write(6,*) ' '
     write(6,*) 'k: ', k
     tmp_1 = (grid%t_1(its,k-1,jts)+t0) * ( (grid%p(its,k-1,jts)+grid%pb(its,k-1,jts)) / 100000.0 )**0.286
     tmp_2 = (grid%t_1(its,k  ,jts)+t0) * ( (grid%p(its,k  ,jts)+grid%pb(its,k  ,jts)) / 100000.0 )**0.286
     tmp_v_1 = tmp_1 * ( 1.0 + rvovrd*moist(its,k-1,jts,P_QV) )
     tmp_v_2 = tmp_2 * ( 1.0 + rvovrd*moist(its,k  ,jts,P_QV) )
     alttde(k) = alttde(k-1) &
               +  ( ( ( 287.0*(tmp_v_1+tmp_v_2)/2.0 ) / 9.81 ) &
                    * log( (grid%p(its,k-1,jts)+grid%pb(its,k-1,jts))/(grid%p(its,k,jts)+grid%pb(its,k,jts)) ) )
     WRITE(6,*) 'alttde = ', alttde(k), alttde(k)-alttde(k-1), grid%p(its,k-1,jts), grid%pb(its,k-1,jts), grid%t_1(its,k-1,jts), t0
   end do
! ub<<

  CASE (squall2d_x)
!  QSS thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

  DO J = jts, min(jde-1,jte)
!   yrad = config_flags%dy*float(j-nyc)/4000.
    yrad = 0.
    DO I = its, min(ide-1,ite)
      xrad = config_flags%dx*float(i-nxc)/4000.
!     xrad = 0.
      DO K = 1, kte-1

!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

        zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
                   +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
        zrad = (zrad-1500.)/1500.
        RAD=SQRT(xrad*xrad+yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
           grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*COS(.5*PI*RAD)**2
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically

      DO k = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k-1)*grid%mub(i,j)+grid%c2h(k-1))+(grid%c1h(k-1)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k-1)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )

        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO
  CASE (squall2d_y)
!  QSS thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

  DO J = jts, min(jde-1,jte)
    yrad = config_flags%dy*float(j-nyc)/4000.
!   yrad = 0.
    DO I = its, min(ide-1,ite)
!     xrad = config_flags%dx*float(i-nxc)/4000.
      xrad = 0.
      DO K = 1, kte-1

!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

        zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
                   +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
        zrad = (zrad-1500.)/1500.
        RAD=SQRT(xrad*xrad+yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
           grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*COS(.5*PI*RAD)**2
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically

      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k-1)*grid%mub(i,j)+grid%c2h(k-1))+(grid%c1h(k-1)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k-1)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )

        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO
  CASE (convrad)
!  random low-level thermal perturbation to kick off convection
        call random_seed
  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

  DO J = jts, min(jde-1,jte)
    DO I = its, min(ide-1,ite)
      DO K = 1, 10

        call RANDOM_NUMBER(rnd)
          grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*(rnd-0.5)
         !  grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*COS(.5*PI*RAD)**2
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
      ENDDO

!  rebalance hydrostatically

      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k-1)*grid%mub(i,j)+grid%c2h(k-1))+(grid%c1h(k-1)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k-1)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )
        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO
  CASE (grav2d_x)
!  cold bubble input  (from straka et al, IJNMF, vol 17, 1993 pp 1-22)

  t_min = grid%t_1(its,kts,jts)
  t_max = t_min
  u_mean = 00.

  xpos = config_flags%dx*nxc - u_mean*900.
  xposml = xpos - config_flags%dx*(ide-1)
  xpospl = xpos + config_flags%dx*(ide-1)

  DO J = jts, min(jde-1,jte)
    DO I = its, min(ide-1,ite)
!      xrad = config_flags%dx*float(i-nxc)/4000.  !  4000 meter horizontal radius
!                                    !  centered in the domain

       xrad = min( abs(config_flags%dx*float(i)-xpos),   &
                   abs(config_flags%dx*float(i)-xposml), &
                   abs(config_flags%dx*float(i)-xpospl))/4000.

      DO K = 1, kte-1

!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

        zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
                   +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
        zrad = (zrad-3000.)/2000. !  2000 meter vertical radius,
                                  !  centered at z=3000,
        RAD=SQRT(xrad*xrad+zrad*zrad)
        IF(RAD <= 1.) THEN

           !  perturbation temperature is 15 C, convert to potential temperature

           delt = -15.0 / ((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**rcp

           grid%T_1(i,k,j)=grid%T_1(i,k,j)+delt*(COS(PI*RAD)+1.0)/2.
           grid%T_2(i,k,j)=grid%T_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF

        t_min = min(t_min, grid%t_1(i,k,j))
        t_max = max(t_max, grid%t_1(i,k,j))
      ENDDO

!  rebalance hydrostatically

      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k-1)*grid%mub(i,j)+grid%c2h(k-1))+(grid%c1h(k-1)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k-1)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )

        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO

  write(6,*) ' min and max theta perturbation ',t_min,t_max

! -- end bubble insert
  CASE (b_wave)
!  thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',tpbub

  DO J = jts, min(jde-1,jte)
    yrad = config_flags%dy*float(j-jde/2-1)/radbub
    DO I = its, min(ide-1,ite)
      xrad = float(i-1)/float(ide-ids)

      DO K = 1, kte-1

!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

        zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
                   +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
        zrad = (zrad-htbub)/radz
        RAD=SQRT(yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
           tp = tpbub*cos(rad*piov2)*cos(rad*piov2)*cos(xrad*2*pi+pi)
           grid%t_1(i,k,j)=grid%t_1(i,k,j)+tp
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*grid%moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically

      DO k  = 2,kte
        grid%ph_1(i,k,j) = grid%ph_1(i,k-1,j) - (grid%dnw(k-1))*(       &
                     ((grid%c1h(k-1)*grid%mub(i,j)+grid%c2h(k-1))+(grid%c1h(k-1)*grid%mu_1(i,j)))*grid%al(i,k-1,j)+ &
                      (grid%c1h(k-1)*grid%mu_1(i,j))*grid%alb(i,k-1,j)  )

        grid%ph_2(i,k,j) = grid%ph_1(i,k,j)
        grid%ph0(i,k,j) = grid%ph_1(i,k,j) + grid%phb(i,k,j)
      ENDDO

    ENDDO
  ENDDO

  CASE (les)
!  thermal perturbation to kick off convection

  write(6,*) ' nxc, nyc for perturbation ',nxc,nyc
  write(6,*) ' delt for perturbation ',delt

! For LES, change the initial random perturbations
! For 2D test, call randx outside I-loop
! For 3D runs, call randx inside both I-J loops

  DO J = jts, min(jde-1,jte)
!   yrad = config_flags%dy*float(j-nyc)/10000.
    yrad = 0.
    DO I = its, min(ide-1,ite)
!     xrad = config_flags%dx*float(i-nxc)/10000.
      xrad = 0.
      call random_number (randx)
      randx = randx - 0.5
!     DO K = 1, kte-1
      DO K = 1, 4

!  No bubbles for LES!
!  put in preturbation theta (bubble) and recalc density.  note,
!  the mass in the column is not changing, so when theta changes,
!  we recompute density and geopotential

!       zrad = 0.5*(grid%ph_1(i,k,j)+grid%ph_1(i,k+1,j)  &
!                  +grid%phb(i,k,j)+grid%phb(i,k+1,j))/g
!       zrad = (zrad-1500.)/1500.
        zrad = 0.
        RAD=SQRT(xrad*xrad+yrad*yrad+zrad*zrad)
        IF(RAD <= 1.) THEN
!          grid%t_1(i,k,j)=grid%t_1(i,k,j)+delt*COS(.5*PI*RAD)**2
           grid%t_1(i,k,j)=grid%t_1(i,k,j)+ 0.1 *randx
           grid%t_2(i,k,j)=grid%t_1(i,k,j)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                        (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
           grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
        ENDIF
      ENDDO

!  rebalance hydrostatically


      DO  kk  = 2,kte
        k = kk - 1
        grid%ph_1(i,kk,j) = grid%ph_1(i,kk-1,j) - (grid%dnw(kk-1))*(       &
                     ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,kk-1,j)+ &
                      (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,kk-1,j)  )

        grid%ph_2(i,kk,j) = grid%ph_1(i,kk,j)
        grid%ph0(i,kk,j) = grid%ph_1(i,kk,j) + grid%phb(i,kk,j)
      ENDDO

    ENDDO
  ENDDO

  END SELECT ideal_pert

#if 0
!===============

!  Test for resting atmosphere

  DO iter_loop = 1, 100
     DO J = jts, min(jde-1,jte)
        DO I = its, min(ide-1,ite)
            DO K = 2, kte-1
               z1 = (grid%phb(i,k+1,j)+grid%ph_1(i,k+1,j))/g
               z2 = (grid%phb(i,k  ,j)+grid%ph_1(i,k  ,j))/g
               z=(z1+z2)/2.
               IF ( z .LT. 10000 ) THEN
                  Nsq = 1.E-4
               ELSE
                  Nsq = 4.E-4
               END IF
               grid%t_1(i,k,j)=(grid%t_1(i,k-1,j)+t0) * EXP(Nsq / g * (z1-z2)) - t0
               grid%t_2(i,k,j)=grid%t_1(i,k,j)
               qvf = 1. + rvovrd*moist(i,k,j,P_QV)
               grid%alt(i,k,j) = (r_d/p1000mb)*(grid%t_1(i,k,j)+t0)*qvf* &
                                 (((grid%p(i,k,j)+grid%pb(i,k,j))/p1000mb)**cvpm)
               grid%al(i,k,j) = grid%alt(i,k,j) - grid%alb(i,k,j)
            ENDDO
      
            DO kk  = 2,kte
              k = kk - 1
              grid%ph_1(i,kk,j) = grid%ph_1(i,kk-1,j) - (grid%dnw(kk-1))*(       &
                           ((grid%c1h(k)*grid%mub(i,j)+grid%c2h(k))+(grid%c1h(k)*grid%mu_1(i,j)))*grid%al(i,kk-1,j)+ &
                            (grid%c1h(k)*grid%mu_1(i,j))*grid%alb(i,kk-1,j)  )
      
              grid%ph_2(i,kk,j) = grid%ph_1(i,kk,j)
              grid%ph0(i,kk,j) = grid%ph_1(i,kk,j) + grid%phb(i,kk,j)
            ENDDO
        ENDDO
     ENDDO
  END DO

!===============
#endif
   IF ( wrf_dm_on_monitor() ) THEN
   k=1
   write(6,*) ' grid%mu_1 from comp ', (grid%c1h(k)*grid%mu_1(1,1))
   write(6,*) ' full state sounding from comp, ph, grid%p, grid%al, grid%t_1, qv '
   do k=1,kde-1
     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1)+grid%phb(1,k,1), &
                                      grid%p(1,k,1)+grid%pb(1,k,1), grid%alt(1,k,1), &
                                      grid%t_1(1,k,1)+t0, moist(1,k,1,P_QV)
   enddo

   write(6,*) ' pert state sounding from comp, grid%ph_1, pp, alp, grid%t_1, qv '
   do k=1,kde-1
     write(6,'(i3,1x,5(1x,1pe10.3))') k, grid%ph_1(1,k,1), &
                                      grid%p(1,k,1), grid%al(1,k,1), &
                                      grid%t_1(1,k,1), moist(1,k,1,P_QV)
   enddo
   ENDIF


  IF ( model_config_rec%ideal_case .EQ. b_wave )THEN
! set v (u already interpolated)

  DO J = jts, jte
  DO I = its, min(ide-1,ite)

    DO K = 1, kte
      grid%v_1(i,k,j) = 0.
      grid%v_2(i,k,j) = grid%v_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO

!  fill out last i row for u

  DO J = jts, min(jde-1,jte)
  DO I = ite, ite

    DO K = 1, kte
      grid%u_1(i,k,j) = grid%u_1(its,k,j)
      grid%u_2(i,k,j) = grid%u_2(its,k,j)
    ENDDO

  ENDDO
  ENDDO

  ELSE
! interp v
  DO J = jts, jte
  DO I = its, min(ide-1,ite)

    IF (j == jds) THEN
      z_at_v = grid%phb(i,1,j)/g
    ELSE IF (j == jde) THEN
      z_at_v = grid%phb(i,1,j-1)/g
    ELSE
      z_at_v = 0.5*(grid%phb(i,1,j)+grid%phb(i,1,j-1))/g
    END IF

    p_surf = interp_0( p_in, zk, z_at_v, nl_in )

    DO K = 1, kte
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top
      grid%v_1(i,k,j) = interp_0( v, p_in, p_level, nl_in )
      grid%v_2(i,k,j) = grid%v_1(i,k,j)
    ENDDO
  ENDDO
  ENDDO
! interp u
  DO J = jts, min(jde-1,jte)
  DO I = its, ite
    IF (i == ids) THEN
      z_at_u = grid%phb(i,1,j)/g
    ELSE IF (i == ide) THEN
      z_at_u = grid%phb(i-1,1,j)/g
    ELSE
      z_at_u = 0.5*(grid%phb(i,1,j)+grid%phb(i-1,1,j))/g
    END IF
    p_surf = interp_0( p_in, zk, z_at_u, nl_in )
    DO K = 1, kte
      p_level = grid%c3h(k)*(p_surf - grid%p_top) + grid%c4h(k) + grid%p_top
      grid%u_1(i,k,j) = interp_0( u, p_in, p_level, nl_in )
      grid%u_2(i,k,j) = grid%u_1(i,k,j)
    ENDDO

  ENDDO
  ENDDO
  ENDIF

!  set w

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte
  DO I = its, min(ide-1,ite)
    grid%w_1(i,k,j) = 0.
    grid%w_2(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO

! ub>>
  write (6,*) ' kts,kte=',kts,kte
  write (6,'(A3,x,7A10)') 'k','z of eta','dz of eta','theta','qv','u'
  do k = kts, kte-1
    write(6,'(I3,x,7F10.3)') k, alttde(k), alttde(k)-alttde(Max(1,k-1)),grid%t_1(its,k,jts)+t0,moist(its,k,jts,P_QV)/(1.+moist(its,k,jts,P_QV))*1000.0,grid%u_1(its,k,jts)
  end do
! ub<<

!  set a few more things

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte-1
  DO I = its, min(ide-1,ite)
    grid%h_diabatic(i,k,j) = 0.
  ENDDO
  ENDDO
  ENDDO

  IF ( wrf_dm_on_monitor() ) THEN
  DO k=1,kte-1
    grid%t_base(k) = grid%t_1(1,k,1)
    grid%qv_base(k) = moist(1,k,1,P_QV)
    grid%u_base(k) = grid%u_1(1,k,1)
    grid%v_base(k) = grid%v_1(1,k,1)
    grid%z_base(k) = 0.5*(grid%phb(1,k,1)+grid%phb(1,k+1,1)+grid%ph_1(1,k,1)+grid%ph_1(1,k+1,1))/g
  ENDDO
  ENDIF
  CALL wrf_dm_bcast_real( grid%t_base , kte )
  CALL wrf_dm_bcast_real( grid%qv_base , kte )
  CALL wrf_dm_bcast_real( grid%u_base , kte )
  CALL wrf_dm_bcast_real( grid%v_base , kte )
  CALL wrf_dm_bcast_real( grid%z_base , kte )

  ideal_surfacet: SELECT CASE ( model_config_rec%ideal_case )
  CASE(hill2d_x, quarter_ss, squall2d_x, squall2d_y, grav2d_x, b_wave)
  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)
     thtmp   = grid%t_2(i,1,j)+t0
     ptmp    = grid%p(i,1,j)+grid%pb(i,1,j)
     temp(1) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,2,j)+t0
     ptmp    = grid%p(i,2,j)+grid%pb(i,2,j)
     temp(2) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,3,j)+t0
     ptmp    = grid%p(i,3,j)+grid%pb(i,3,j)
     temp(3) = thtmp * (ptmp/p1000mb)**rcp

     grid%tsk(I,J)=grid%cf1*temp(1)+grid%cf2*temp(2)+grid%cf3*temp(3)
     grid%tmn(I,J)=grid%tsk(I,J)-0.5
  ENDDO
  ENDDO

  CASE(seabreeze2d_x)
  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)
     thtmp   = grid%t_2(i,1,j)+t0
     ptmp    = grid%p(i,1,j)+grid%pb(i,1,j)
     temp(1) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,2,j)+t0
     ptmp    = grid%p(i,2,j)+grid%pb(i,2,j)
     temp(2) = thtmp * (ptmp/p1000mb)**rcp
     thtmp   = grid%t_2(i,3,j)+t0
     ptmp    = grid%p(i,3,j)+grid%pb(i,3,j)
     temp(3) = thtmp * (ptmp/p1000mb)**rcp

!    grid%tsk(I,J)=grid%cf1*temp(1)+grid%cf2*temp(2)+grid%cf3*temp(3)
     grid%tmn(I,J)=grid%tsk(I,J)-0.5
  ENDDO
  ENDDO

  CASE(convrad)
  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)
!         grid%tsk(i,j) = theta_surf * (p_surf/p1000mb)**rcp
         grid%tsk(i,j) = sst_input * (p_surf/p1000mb)**rcp
         grid%tmn(i,j) = grid%tsk(i,j)
  ENDDO
  ENDDO

  CASE(les)
  DO J = jts, min(jde-1,jte)
  DO I = its, min(ide-1,ite)
         grid%tsk(i,j) = theta_surf * (p_surf/p1000mb)**rcp
         grid%tmn(i,j) = grid%tsk(i,j)-0.5
  ENDDO
  ENDDO

  END SELECT ideal_surfacet


  trajectories: SELECT CASE ( model_config_rec%ideal_case )
  CASE (quarter_ss)
  !  Template for initializing trajectories.  The i, j, and k starting locations
  !  are specified.  Right now, a small plane in the middle of the domain is
  !  selected.

  grid%traj_i    = -9999
  grid%traj_j    = -9999
  grid%traj_k    = -9999
  grid%traj_lat  = -9999
  grid%traj_long = -9999

  IF (config_flags%num_traj .gt. 0 .and. config_flags%traj_opt .gt. 0) THEN
     icount = 1
     DO j = (jde + jds)/2 - 2, (jde + jds)/2 + 2, 1
        DO i = (ide + ids)/2 - 2, (ide + ids)/2 + 2, 1
           IF ( its .LE. i    .and. ite .GE. i   .and.  jts .LE. j    .and. jte .GE. j ) THEN
              grid%traj_i   (icount) = i
              grid%traj_j   (icount) = j
              grid%traj_k   (icount) = 10
              grid%traj_lat (icount) = grid%xlat(i,j)
              grid%traj_long(icount) = grid%xlong(i,j)
           END IF

#ifdef DM_PARALLEL
           grid%traj_i   (icount) = wrf_dm_max_real ( grid%traj_i   (icount) )
           grid%traj_j   (icount) = wrf_dm_max_real ( grid%traj_j   (icount) )
           grid%traj_k   (icount) = wrf_dm_max_real ( grid%traj_k   (icount) )
           grid%traj_lat (icount) = wrf_dm_max_real ( grid%traj_lat (icount) )
           grid%traj_long(icount) = wrf_dm_max_real ( grid%traj_long(icount) )
#endif

           icount = icount + 1
           IF (icount .GT. config_flags%num_traj) THEN
              EXIT
           END IF
        END DO
     END DO
  END IF
  END SELECT trajectories

  tracers: SELECT CASE ( config_flags%tracer_opt )
  CASE (tracer_test1)

  DO J = jts, min(jde-1,jte)
  DO K = kts, kte-1
  DO I = its, min(ide-1,ite)
    grid%h_diabatic(i,k,j) = 0.
    if(k.eq.kts)tracer(i,k,j,p_tr17_1)=1.
    if(k.eq.kts.and.grid%xland(i,j).lt.1.5)tracer(i,k,j,p_tr17_2)=1.
    if(k.eq.kts.and.grid%xland(i,j).gt.1.5)tracer(i,k,j,p_tr17_3)=1.
    if(k.le.5)tracer(i,k,j,p_tr17_4)=1.
    if(k.le.5.and.grid%xland(i,j).lt.1.5)tracer(i,k,j,p_tr17_5)=1.
    if(k.le.5.and.grid%xland(i,j).gt.1.5)tracer(i,k,j,p_tr17_6)=1.
    if(k.le.10)tracer(i,k,j,p_tr17_7)=1.
    if(k.le.10.and.k.gt.5)tracer(i,k,j,p_tr17_8)=1.
  ENDDO
  ENDDO
  ENDDO

  END SELECT tracers

      !  Save the dry perturbation potential temperature.

      DO j = jts, min(jde-1,jte)
         DO k = kts, kte
            DO i = its, min(ide-1,ite)
               grid%th_phy_m_t0(i,k,j) = grid%t_2(i,k,j)
            END DO
         END DO
      END DO

  !  Turn dry potential temperature into moist potential temperature
      !  at the very end of this routine
      !  This field will be in the model IC and and used to construct the 
      !  BC file.

      IF ( ( config_flags%use_theta_m .EQ. 1 ) .AND. (P_Qv .GE. PARAM_FIRST_SCALAR) ) THEN
      DO J  = jts, min(jde-1,jte)
         DO K = kts, kte-1
            DO I = its, min(ide-1,ite)
               grid%t_2(i,k,j) = ( grid%t_2(i,k,j) + T0 ) * (1. + (R_v/R_d) * moist(i,k,j,p_qv)) - T0
            END DO
         END DO
      END DO
      ENDIF

  RETURN

 END SUBROUTINE init_domain_rk

   SUBROUTINE init_module_initialize
   END SUBROUTINE init_module_initialize

!---------------------------------------------------------------------

!  test driver for get_sounding
!
!      implicit none
!      integer n
!      parameter(n = 1000)
!      real zk(n),p(n),theta(n),rho(n),u(n),v(n),qv(n),pd(n)
!      logical dry
!      integer nl,k
!
!      dry = .false.
!      dry = .true.
!      call get_sounding( zk, p, pd, theta, rho, u, v, qv, dry, n, nl )
!      write(6,*) ' input levels ',nl
!      write(6,*) ' sounding '
!      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
!      do k=1,nl
!        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), pd(k), theta(k), rho(k), u(k), v(k), qv(k)
!      enddo
!      end
!
!---------------------------------------------------------------------------

      subroutine get_sounding( zk, p, p_dry, theta, rho, &
                               u, v, qv, etal, dry, nl_max, nl_in, th_surf )
      implicit none

      integer nl_max, nl_in
      real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
           u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max), etal(nl_max)
      logical dry

      integer n, iz
      parameter(n=1000)
      logical debug
      parameter( debug = .true.)
      character*256 message

! input sounding data

      real p_surf, th_surf, qv_surf
      real pi_surf, pi(n)
      real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)

! diagnostics

      real rho_surf, p_input(n), rho_input(n)
      real pm_input(n)  !  this are for full moist sounding

! local data

      real r
      parameter (r = r_d)
      integer k, it, nl
      real qvf, qvf1, dz

!  first, read the sounding

      call read_sounding( p_surf, th_surf, qv_surf, &
                          h_input, th_input, qv_input, u_input, v_input, pm_input,n, nl, debug )

!        iz = 1
!        do k=2,nl
!          if(h_input(k) .lt. 12000.) iz = k
!        enddo
!        write(6,*) " tropopause ",iz,h_input(iz)
!        if(dry) then
!        write(6,*) ' nl is ',nl
!        do k=1,nl
!          th_input(k) = th_input(k)+10.+10*float(k)/nl
!        enddo
!        write(6,*) ' finished adjusting theta '
!        endif

!        do k=1,nl
!          u_input(k) = 2*u_input(k)
!        enddo
!
!      end if

      if(dry) then
       do k=1,nl
         qv_input(k) = 0.
       enddo
      endif

      if(debug) write(6,*) ' number of input levels = ',nl

        nl_in = nl
        if(nl_in .gt. nl_max ) then
          write(6,*) ' too many levels for input arrays ',nl_in,nl_max
          call wrf_error_fatal ( ' too many levels for input arrays ' )
        end if

!  compute diagnostics,
!  first, convert qv(g/kg) to qv(g/g)

      do k=1,nl
        qv_input(k) = 0.001*qv_input(k)
      enddo

      p_surf = 100.*p_surf  ! convert to pascals
      qvf = 1. + rvovrd*qv_input(1)
      rho_surf = 1./((r/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
      pi_surf = (p_surf/p1000mb)**(r/cp)

      if(debug) then
        write(6,*) ' surface density is ',rho_surf
        write(6,*) ' surface pi is      ',pi_surf
      end if


!  integrate moist sounding hydrostatically, starting from the
!  specified surface pressure
!  -> first, integrate from surface to lowest level

          qvf = 1. + rvovrd*qv_input(1)
          qvf1 = 1. + qv_input(1)
          rho_input(1) = rho_surf
          dz = h_input(1)
          do it=1,10
            IF ( .not. use_snd_plevels ) THEN
            pm_input(1) = p_surf &
                    - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1
            ENDIF
            rho_input(1) = 1./((r/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
          enddo

     IF ( use_snd_plevels ) etal(1) = 1.

! integrate up the column

          do k=2,nl
            rho_input(k) = rho_input(k-1)
            dz = h_input(k)-h_input(k-1)
            qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
            qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here
 
            IF ( use_snd_plevels ) THEN ! using pressure levels from input sounding
              etal(k) = (pm_input(k) - pm_input(nl))/(pm_input(1) - pm_input(nl))
            ENDIF
 
            do it=1,10
             IF ( .not. use_snd_plevels ) THEN ! calculate pressure by iteration (not using input pressure)
              pm_input(k) = pm_input(k-1) &
                      - 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
              IF(pm_input(k) .LE. 0. )THEN
                CALL wrf_message("Integrated pressure has gone negative - too cold for chosen height")
                WRITE(message,*)'k,pm_input(k),h_input(k),th_input(k) = ',k,pm_input(k),h_input(k),th_input(k)
                CALL wrf_error_fatal ( message )
              ENDIF
              ENDIF
              rho_input(k) = 1./((r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
            enddo
          enddo

          IF ( use_snd_zlevels ) THEN
            etal(1) = 1.
            do k=2,nl
              etal(k) = (pm_input(k) - pm_input(nl))/(pm_input(1) - pm_input(nl))
            enddo
          ENDIF

!  we have the moist sounding

!  next, compute the dry sounding using p at the highest level from the
!  moist sounding and integrating down.

        p_input(nl) = pm_input(nl)

          do k=nl-1,1,-1
            dz = h_input(k+1)-h_input(k)
            p_input(k) = p_input(k+1) + 0.5*dz*(rho_input(k)+rho_input(k+1))*g
          enddo


        do k=1,nl

          zk(k) = h_input(k)
          p(k) = pm_input(k)
          p_dry(k) = p_input(k)
          theta(k) = th_input(k)
          rho(k) = rho_input(k)
          u(k) = u_input(k) - ugrid
          v(k) = v_input(k) - vgrid
          qv(k) = qv_input(k)

        enddo

     if(debug) then
      write(6,*) ' sounding '

      IF ( use_snd_plevels ) THEN
      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) eta '
      do k=1,nl
        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k),etal(k)
      enddo

      ELSE

      write(6,*) '  k  height(m)  press (Pa) pd(Pa) theta (K) den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
      do k=1,nl
        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
      enddo

      ENDIF

     end if

      end subroutine get_sounding

!-------------------------------------------------------

      subroutine read_sounding( ps,ts,qvs,h,th,qv,u,v,p,n,nl,debug )
      implicit none
      integer n,nl
      real ps,ts,qvs,h(n),th(n),qv(n),u(n),v(n),p(n)
      real h0,pfac
      logical end_of_file
      logical debug
!      real :: hw(300),thw(300),qvw(300),uw(300),vw(300)
!      integer :: nltmp,j

      integer k

      pfac = 1.0
      open(unit=10,file='input_sounding',form='formatted',status='old')
      rewind(10)
      read(10,*) ps, ts, qvs
      IF ( ps > 10000. ) pfac = 0.01
      ps = ps*pfac
      if(debug) then
        write(6,*) ' input sounding surface parameters '
        write(6,*) ' surface pressure (mb) ',ps
        write(6,*) ' surface pot. temp (K) ',ts
        write(6,*) ' surface mixing ratio (g/kg) ',qvs
      end if

      end_of_file = .false.
      k = 0

      IF ( use_snd_plevels ) THEN
      
      do while (.not. end_of_file)

        read(10,*,end=200) h(k+1), th(k+1), qv(k+1), u(k+1), v(k+1), p(k+1)
        IF ( k == 0 ) h0 = h(k+1)
!        p(k+1) = p(k+1)*pfac
        k = k+1
        h(k) = h(k) - h0
        if(debug) write(6,'(1x,i3,5(1x,e10.3))') k, h(k), th(k), qv(k), u(k), v(k), p(k)
        go to 210
 200    end_of_file = .true.
 210    continue
      enddo
      
      ELSE
      
      do while (.not. end_of_file)

        read(10,*,end=100) h(k+1), th(k+1), qv(k+1), u(k+1), v(k+1)
        IF ( use_snd_zlevels .and. k == 0 .and. h(k+1) > 0.0001 ) THEN
          ! need to add a zero level.
           h0 = 0.0
           k = k+1
           h(k+1) = h(k)
           qv(k+1) = qv(k)
           u(k+1) = u(k)
           v(k+1) = v(k)
        ENDIF

        k = k+1


        if(debug) write(6,'(1x,i3,5(1x,e10.3))') k, h(k), th(k), qv(k), u(k), v(k)
        
        go to 110
 100    end_of_file = .true.
 110    continue
      enddo
              
        IF ( use_snd_zlevels ) THEN ! reading in a commas-produced sounding
                                    ! Sounding has surface (w-point) and then 
                                    ! scalar points, so we need to fake a top
                                    ! w-point to use as eta level 0
        h(k+1) = h(k) + 0.5*(h(k)-h(k-1))
        th(k+1) = th(k) + 0.5*(th(k)-th(k-1))
        
        k = k+1
        u(k) = u(k-1)
        v(k) = v(k-1)
        qv(k) = qv(k-1)

        if(debug) write(6,'(1x,i3,5(1x,e10.3))') k, h(k), th(k), qv(k), u(k), v(k)
        
       ENDIF

      ENDIF

      nl = k

      close(unit=10,status = 'keep')

      end subroutine read_sounding

!------------------------------------------------------------------

    subroutine get_sounding_b_wave( zk, p, p_dry, theta, rho,       &
                             u, v, qv, dry, nl_max, nl_in,  &
                             u_jet, rho_jet, th_jet, z_jet, &
                             nz_jet, ny_jet, j_point, debug )
    implicit none

    integer nl_max, nl_in
    real zk(nl_max), p(nl_max), theta(nl_max), rho(nl_max), &
         u(nl_max), v(nl_max), qv(nl_max), p_dry(nl_max)
    logical dry

    integer nz_jet, ny_jet, j_point
    real, dimension(nz_jet, ny_jet) :: u_jet, rho_jet, th_jet, z_jet

    integer n
    parameter(n=1000)
    logical debug

! input sounding data

    real p_surf, th_surf, qv_surf
    real pi_surf, pi(n)
    real h_input(n), th_input(n), qv_input(n), u_input(n), v_input(n)

! diagnostics

    real rho_surf, p_input(n), rho_input(n)
    real pm_input(n)  !  this are for full moist sounding

! local data

    real r
    parameter (r = r_d)
    integer k, it, nl
    real qvf, qvf1, dz

!  first, read the sounding

!    call read_sounding( p_surf, th_surf, qv_surf, &
!                          h_input, th_input, qv_input, u_input, v_input,n, nl, debug )

   call calc_jet_sounding( p_surf, th_surf, qv_surf,                             &
                           h_input, th_input, qv_input, u_input, v_input,        &
                           n, nl, debug, u_jet, rho_jet, th_jet, z_jet, j_point, &
                           nz_jet, ny_jet, dry                                  )

   nl = nz_jet

    if(dry) then
     do k=1,nl
       qv_input(k) = 0.
     enddo
    endif

    if(debug) write(6,*) ' number of input levels = ',nl

      nl_in = nl
      if(nl_in .gt. nl_max ) then
        write(6,*) ' too many levels for input arrays ',nl_in,nl_max
        call wrf_error_fatal ( ' too many levels for input arrays ' )
      end if

!  compute diagnostics,
!  first, convert qv(g/kg) to qv(g/g)
!
!      do k=1,nl
!        qv_input(k) = 0.001*qv_input(k)
!      enddo
!      p_surf = 100.*p_surf  ! convert to pascals

    qvf = 1. + rvovrd*qv_input(1)
    rho_surf = 1./((r/p1000mb)*th_surf*qvf*((p_surf/p1000mb)**cvpm))
    pi_surf = (p_surf/p1000mb)**(r/cp)

    if(debug) then
      write(6,*) ' surface density is ',rho_surf
      write(6,*) ' surface pi is    ',pi_surf
    end if


!  integrate moist sounding hydrostatically, starting from the
!  specified surface pressure
!  -> first, integrate from surface to lowest level

        qvf = 1. + rvovrd*qv_input(1)
        qvf1 = 1. + qv_input(1)
        rho_input(1) = rho_surf
        dz = h_input(1)
        do it=1,10
          pm_input(1) = p_surf &
                  - 0.5*dz*(rho_surf+rho_input(1))*g*qvf1
          rho_input(1) = 1./((r/p1000mb)*th_input(1)*qvf*((pm_input(1)/p1000mb)**cvpm))
        enddo

! integrate up the column

        do k=2,nl
          rho_input(k) = rho_input(k-1)
          dz = h_input(k)-h_input(k-1)
          qvf1 = 0.5*(2.+(qv_input(k-1)+qv_input(k)))
          qvf = 1. + rvovrd*qv_input(k)   ! qv is in g/kg here

          do it=1,10
            pm_input(k) = pm_input(k-1) &
                    - 0.5*dz*(rho_input(k)+rho_input(k-1))*g*qvf1
            rho_input(k) = 1./((r/p1000mb)*th_input(k)*qvf*((pm_input(k)/p1000mb)**cvpm))
          enddo
        enddo

!  we have the moist sounding

!  next, compute the dry sounding using p at the highest level from the
!  moist sounding and integrating down.

        p_input(nl) = pm_input(nl)

          do k=nl-1,1,-1
            dz = h_input(k+1)-h_input(k)
            p_input(k) = p_input(k+1) + 0.5*dz*(rho_input(k)+rho_input(k+1))*g
          enddo


        do k=1,nl

          zk(k) = h_input(k)
          p(k) = pm_input(k)
          p_dry(k) = p_input(k)
          theta(k) = th_input(k)
          rho(k) = rho_input(k)
          u(k) = u_input(k)
          v(k) = v_input(k)
          qv(k) = qv_input(k)

        enddo

     if(debug) then
      write(6,*) ' sounding '
      write(6,*) '  k  height(m)  press (Pa)   pd(Pa)   theta (K)  den(kg/m^3)  u(m/s)     v(m/s)    qv(g/g) '
      do k=1,nl
        write(6,'(1x,i3,8(1x,1pe10.3))') k, zk(k), p(k), p_dry(k), theta(k), rho(k), u(k), v(k), qv(k)
      enddo

     end if

     end subroutine get_sounding_b_wave

!------------------------------------------------------------------

  subroutine calc_jet_sounding( p_surf, th_surf, qv_surf,      &
                                h, th, qv, u, v, n, nl, debug, &
                                u_jet, rho_jet, th_jet, z_jet, &
                                jp, nz_jet, ny_jet, dry       )
  implicit none
  integer :: n, nl, jp, nz_jet, ny_jet

  real, dimension(nz_jet, ny_jet) :: u_jet, rho_jet, th_jet, z_jet
  real, dimension(n) :: h,th,qv,u,v
  real :: p_surf, th_surf, qv_surf
  logical :: debug, dry

  real, dimension(1:nz_jet) :: rho, rel_hum, p
  integer :: k

!  some local stuff

  real :: tmppi, es, qvs, temperature

!  get sounding from column jp

   do k=1,nz_jet
     h(k)  = z_jet(k,jp)
     th(k) = th_jet(k,jp)
     qv(k) = 0.
     rho(k) = rho_jet(k,jp)
     u(k) = u_jet(k,jp)
     v(k) = 0.
   enddo

   if (.not.dry) then
     DO k=1,nz_jet
       if(h(k) .gt. 8000.) then
         rel_hum(k)=0.1
       else
         rel_hum(k)=(1.-0.90*(h(k)/8000.)**1.25)
       end if
       rel_hum(k) = min(0.7,rel_hum(k))
     ENDDO
   else
     do k=1,nz_jet
       rel_hum(k) = 0.
     enddo
   endif

!  next, compute pressure

   do k=1,nz_jet
     p(k) = p1000mb*(R_d*rho(k)*th(k)/p1000mb)**cpovcv
   enddo

!  here we adjust for fixed moisture profile

     IF (.not.dry)  THEN

!  here we assume the input theta is th_v, so we reset theta accordingly

       DO k=1,nz_jet
         tmppi=(p(k)/p1000mb)**rcp
         temperature = tmppi*th(k)
         if (temperature .gt. svpt0) then
            es  = 1000.*svp1*exp(svp2*(temperature-svpt0)/(temperature-svp3))
            qvs = ep_2*es/(p(k)-es)
         else
            es  = 1000.*svp1*exp( 21.8745584*(temperature-273.16)/(temperature-7.66) )
            qvs = ep_2*es/(p(k)-es)
         endif
         qv(k) = rel_hum(k)*qvs
         th(k) = th(k)/(1.+.61*qv(k))
       ENDDO

     ENDIF

!  finally, set the surface data. We'll just do a simple extrapolation

   p_surf = 1.5*p(1) - 0.5*p(2)
   th_surf = 1.5*th(1) - 0.5*th(2)
   qv_surf = 1.5*qv(1) - 0.5*qv(2)

   end subroutine calc_jet_sounding

!---------------------------------------------------------------------

 SUBROUTINE read_input_jet( u, r, t, zk, nz, ny )
 implicit none

 integer, intent(in) :: nz,ny
 real, dimension(nz,ny), intent(out) :: u,r,t,zk
 integer :: ny_in, nz_in, j,k
 real, dimension(ny,nz) :: field_in
 character*256 message

! this code assumes it is called on processor 0 only

   OPEN(unit=10, file='input_jet', form='unformatted', status='old' )
   REWIND(10)
   read(10) ny_in,nz_in
   if((ny_in /= ny ) .or. (nz_in /= nz)) then
     write(message,*) ' error in input jet dimensions '
     CALL wrf_message (message)
     write(message,*) ' ny, ny_input, nz, nz_input ', ny, ny_in, nz,nz_in
     CALL wrf_message (message)
     write(message,*) ' error exit '
     CALL wrf_message (message)
     call wrf_error_fatal ( ' error in input jet dimensions ' )
   end if
   read(10) field_in
   do j=1,ny
   do k=1,nz
     u(k,j) = field_in(j,k)
   enddo
   enddo
   read(10) field_in
   do j=1,ny
   do k=1,nz
     t(k,j) = field_in(j,k)
   enddo
   enddo

   read(10) field_in
   do j=1,ny
   do k=1,nz
     r(k,j) = field_in(j,k)
   enddo
   enddo

   do j=1,ny
   do k=1,nz
     zk(k,j) = 125. + 250.*float(k-1)
   enddo
   enddo

 end subroutine read_input_jet

!-------------------------------------------------------

      subroutine read_soil( n,nl,tmn,tsk,zs,tslb,smois )
      implicit none
      integer n,nl
      real tmn,tsk
      real zs(n),tslb(n),smois(n)
      logical end_of_file
      logical debug 

      integer k
    
      debug = .true.

      open(unit=11,file='input_soil',form='formatted',status='old')
      rewind(11)

      read(11,*) zs(1),tmn,tsk

      if(debug) then
        write(6,*) ' input deep soil temperature (K) ',tmn
        write(6,*) ' input skin temperature (K) ',tsk
      end if

      end_of_file = .false.
      k = 0

      do while (.not. end_of_file)

        read(11,*,end=100) zs(k+1), tslb(k+1), smois(k+1)
        k = k+1
        if(debug) write(6,'(1x,i3,3(1x,f16.7))') k, zs(k), tslb(k), smois(k)
        go to 110
 100    end_of_file = .true.
 110    continue
      enddo

      nl = k

      close(unit=11,status = 'keep')

      end subroutine read_soil

! #####################################################################
!
!
! #####################################################################
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  STRAKAS ATMOSPHERIC MODEL  (SAM)
!    Designed by Jerry M. Straka
!
!  RNDNUM(ISEED) (must define iseed=-1 on intialization 
!                 in user routine)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!2345678901234567890123456789012345678901234567890123456789012345678912
!
!  start function for random numbers....
!
      real function rndnum(iseed)
!
      implicit none
!
!      real rndnum
!
      integer :: iseed
      real    ambig
      parameter (ambig=4000000)
      real    amseed
      parameter (amseed=1618033)
      real    amz,amj,amk
      parameter (amz=0)
!
      real, save ::    ama(55)
!
      real       fac
      parameter (fac=1./ambig)
!
      integer, save :: iff = 0, inext, inextp
      integer :: n, nn, m
!      save iff, inext, inextp, ama
!
!      data iff /0/
!
      if ( iseed .lt. 0 .or. iff .eq. 0 ) then
      iff = 1
      amj = amseed - abs(iseed)
      amj = mod(amj,ambig)
      ama(55) = amj
      amk = 1
      do n = 1,54
      nn = mod(21*n,55)
      ama(nn) = amk
      amk = amj-amk
      if (amk.lt.amz) amk=amk+ambig
      amj=ama(nn)
      end do
!
      do m = 1,4
      do n = 1,55
      ama(n)=ama(n)-ama(1+mod(n+30,55))
      if (ama(n).lt.amz) ama(n)=ama(n)+ambig
      end do
      end do
!
      inext = 0
      inextp = 31
      iseed = 1
!
      end if
!
! find random number
!
      inext=inext+1
      if (inext.eq.56) inext=1
      inextp=inextp+1
      if (inextp.eq.56) inextp=1
      amj=ama(inext)-ama(inextp)
      if (amj.lt.amz) amj=amj+ambig
      ama(inext)=amj
      rndnum=amj*fac
!
      return
      end function rndnum
!
!
!
!
!  end of random number function
!
! 
!
! #####################################################################


END MODULE module_initialize_ideal
