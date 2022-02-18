program cplbypass_metload
    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the input data from the coupler to the land model 
    !
  !
  !===============================================================================

    ! !USES:
    use netcdf
    implicit none

    include 'mpif.h'

    !------------------------------------------------------------------------

    integer myid, np, ierr
    character(len=150) forcdir
    !
    ! domain
    character(len=150) mydomain, fname
    integer :: i, j, g, ni, nj
    double precision, pointer :: xc(:,:), yc(:,:)
    integer, pointer:: mask(:,:)
    real, pointer :: latc(:)      => null()                         ! latitude of grid cell (deg)
    real, pointer :: lonc(:)      => null()                         ! longitude of grid cell (deg)
    integer, pointer :: gindx(:)  => null()                         ! 1-D grid cell index from 1
    integer :: ng_mpi, ng_res

    ! metdata from cpl_bypass_full
    character(len=20) :: met_type
    integer :: startyear_met                                        !staring driver met year
    integer :: endyear_met                                          !end driver met year
    integer*2, pointer :: atm_input            (:,:,:,:) => null()  !Single-site meteorological input
    real, pointer :: add_offsets                     (:) => null()  !offsets for compressed met drivers
    real, pointer :: scale_factors                   (:) => null()  !scale factors for compressed met drivers
    real, pointer :: timeres                         (:) => null()  !time resolution of driver met (hours)
    real, pointer :: var_offset                  (:,:,:) => null()  !correction offset for grid->site driver met (monthly)
    real, pointer :: var_mult                    (:,:,:) => null()  !correction factor for grid->site driver met (monthly)
    integer,  pointer :: timelen                     (:) => null()  !length of input meteorology
    integer,  pointer :: timelen_spinup              (:) => null()  !length of spinup meteorology
    integer,  pointer :: tindex                  (:,:,:) => null()  !current index for meteorolgoical data
    real, pointer :: npf                            (:)  => null()  !number of model timesteps per forcing timestep
    !
    ! locals
    character(len=3) zst
    character(len=300)  :: metdata_fname
    integer  :: topo,m,thism,nstep
    integer  :: thisng, num, nu_nml, nml_error
    integer  :: ng_all(100000)
    real     :: timetemp(2)
    real     :: latixy(500000), longxy(500000)
    integer ::  ncid, met_ncids(14), mask_ncid, thisncid, ng, tm, varid, dimid
    integer ::  aindex(2), starti(3), counti(3)
    integer ::  grid_map(500000), zone_map(500000)
    integer ::  met_nvars, nyears, starti_site, endi_site
    real    :: smap05_lat(360), smap05_lon(720)
    real    :: smapt62_lat(94), smapt62_lon(192)
    real    :: smap2_lat(96), smap2_lon(144)
    real    :: thisdist, mindist, thislon
    real    :: tempndep(1,1,158), thiscalday, wt1(14), wt2(14), thisdoy
    real    :: site_metdata(14,12)
    real    :: var_month_mean(12)
    integer*2 :: temp(1,500000)
    integer :: xtoget, ytoget, thisx, thisy
    character(len=200) metsource_str, thisline
    integer :: av, v, n, nummetdims, g3, gtoget, ztoget, line, mystart, tod_start, thistimelen  
    character(len=20) aerovars(14), metvars(14)

    real      :: ival  = 0.0
    real      :: ival_float = 0.0
    real      :: step_size
    integer   :: ival_int = 0
    integer*2 :: ival_short = 0

    integer   :: begg, endg
    logical   :: has_zonefile
    double precision      :: t0, t1, t2

    !------------------------------------------------------------------------

    call MPI_init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

    t0  = MPI_Wtime()

    ! -------------------------------------------------------------
    ! USER-SPECIFIED forcing data input directory, file header, domain, and period

    !Set the input data path
    forcdir = '.' // & !'/nfs/data/ccsi/f9y/GSWP3_daymet' // &
                '/TILE138/cpl_bypass_full'

    ! domain/mask
    mydomain = 'domain.nc'

    if (trim(mydomain) /= '') then
        fname = trim(forcdir) // '/' // trim(mydomain)
        ierr = nf90_open(trim(fname), NF90_NOWRITE, ncid)
        if (myid==0) print*, 'open domain file: ', trim(fname), ' - ', trim(nf90_strerror(ierr))
        ierr = nf90_inq_dimid(ncid, "ni", dimid)
        ierr = nf90_inquire_dimension(ncid, dimid, len = ni)
        if (ierr /= 0) then
           print *, 'error in dimension ni - ', trim(nf90_strerror(ierr))
           stop
        endif
        ierr = nf90_inq_dimid(ncid, "nj", dimid)
        ierr = nf90_inquire_dimension(ncid, dimid, len = nj)
        if (ierr /= 0) then
           print *, 'error in dimension nj - ', trim(nf90_strerror(ierr))
           stop
        endif

        allocate(xc(ni,nj))
        allocate(yc(ni,nj))
        allocate(mask(ni,nj))
        ierr = nf90_inq_varid(ncid, 'xc', varid)
        ierr = nf90_get_var(ncid, varid, xc)
        ierr = nf90_inq_varid(ncid, 'yc', varid)
        ierr = nf90_get_var(ncid, varid, yc)
        ierr = nf90_inq_varid(ncid, 'mask', varid)
        ierr = nf90_get_var(ncid, varid, mask)
        ierr = nf90_close(ncid)

        allocate(lonc(ni*nj))
        allocate(latc(ni*nj))
        allocate(gindx(ni*nj))
        ng = 0
        do i=1,ni
          do j=1,nj
            if (mask(ni,nj) ==1 ) then
              ng = ng + 1
              lonc(ng) = xc(ni,nj)
              latc(ng) = yc(ni,nj)
              gindx(ng) = ng
            end if
          end do
        end do

        ! domain round-robin
        ng_res = mod(ng, np)
        ng_mpi = (ng-ng_res)/np
        begg = myid*ng_mpi + 1 + min(myid+1, ng_res)
        endg = (myid+1)*ng_mpi + min(myid+1, ng_res)
     endif


      met_type = 'GSWP3_daymet4'
     !get met-grid lat/lon information, zone mappings of cpl_bypass_full met. datasets

     fname = trim(forcdir) // '/zone_mappings.txt'
     inquire(file=fname, exist=has_zonefile)
     if (has_zonefile) then
         open(unit=13, file=fname)
     else if (myid==0) then
       print*, 'NO zone file: ', trim(fname)
     end if

     ng = 0     !number of points
     do v=1,500000
         read(13,*, end=10) longxy(v), latixy(v), zone_map(v), grid_map(v)
         ng = ng + 1
     end do
10   continue
     close(unit=13)

    !-------------------------------------------------------------------------------------------
    t1  = MPI_Wtime()
    write (100+myid, *) 't1: ', t1-t0,' @rank ', myid

    allocate(timelen                            (1:14))        ; timelen                       (:)   = ival_int
    allocate(timelen_spinup                     (1:14))        ; timelen_spinup                (:)   = ival_int
    allocate(tindex               (begg:endg,1:14,1:2))        ; tindex                    (:,:,:)   = ival_int
    allocate(npf                                (1:14))        ; npf                           (:)   = ival
    allocate(add_offsets                        (1:14))        ; add_offsets                   (:)   = ival_float
    allocate(scale_factors                      (1:14))        ; scale_factors                 (:)   = ival_float
    allocate(timeres                            (1:14))        ; timeres                       (:)   = ival
    allocate(var_offset              (14,begg:endg,12))        ; var_offset                (:,:,:)   = ival
    allocate(var_mult                (14,begg:endg,12))        ; var_mult                  (:,:,:)   = ival
    !
    !---------------------------------------------------------------------------

    thisng = endg - begg + 1
    step_size = 1800.0d0
    do g = begg,endg
         i = 1 + (g - begg)
       
  !-----------------------------------Meteorological forcing  -----------------------------------

          met_nvars=7
 
          metvars(1) = 'TBOT'
          metvars(2) = 'PSRF'
          metvars(3) = 'QBOT'
          metvars(4) = 'FSDS'
          metvars(5) = 'PRECTmms'
          metvars(6) = 'WIND'
          metvars(7) = 'FLDS'

          !set defaults
          startyear_met       = 1980
          endyear_met         = 2014
          nyears              = endyear_met - startyear_met + 1

          !Figure out the closest point and which zone file to open
          mindist=99999
          do g3 = 1,ng
              ! in CPL_BYPASS met dataset, longitude is in format of 0-360, but 'lonc(g)' may or may not.
              if (lonc(g) .lt. 0) then
                if (longxy(g3) >= 180) longxy(g3) = longxy(g3)-360.0d0
              else if (lonc(g) .ge. 180) then
                if (longxy(g3) < 0) longxy(g3) = longxy(g3) + 360.0d0
              end if
              thisdist = 100*((latixy(g3) - latc(g))**2 + &
                              (longxy(g3) - lonc(g))**2)**0.5
              if (thisdist .lt. mindist) then 
                mindist = thisdist
                ztoget = zone_map(g3)
                gtoget = grid_map(g3)
              end if
          end do

          do v=1,met_nvars
            write(zst, '(I3)') 100+ztoget
            !daymet v4 with GSWP3 v2 for NA with user-defined zone-mappings.txt
            metdata_fname = 'GSWP3_daymet4_' // trim(metvars(v)) // '_1980-2014_z' // zst(2:3) // '.nc'
            fname = trim(forcdir) // '/' // trim(metdata_fname)
            ierr = nf90_open(fname, NF90_NOWRITE, met_ncids(v))
            if (myid==0 .and. i==1) print *, 'open met file: ', trim(fname), ' - ', trim(nf90_strerror(ierr))
       
            !get timestep information
            ierr = nf90_inq_dimid(met_ncids(v), 'DTIME', dimid)
            ierr = nf90_Inquire_Dimension(met_ncids(v), dimid, len = timelen(v))

            starti(1) = 1
            counti(1) = 2
            ierr = nf90_inq_varid(met_ncids(v), 'DTIME', varid)
            ierr = nf90_get_var(met_ncids(v), varid, timetemp, starti(1:1), counti(1:1))
            timeres(v)        = (timetemp(2)-timetemp(1))*24.0d0
            npf(v)            = 86400d0*(timetemp(2)-timetemp(1))/step_size
            timelen_spinup(v) = nyears*(365*nint(24./timeres(v)))
    
            ierr = nf90_inq_varid(met_ncids(v), trim(metvars(v)), varid)

            !get the conversion factors
            ierr = nf90_get_att(met_ncids(v), varid, 'scale_factor', scale_factors(v))
            if (ierr .ne. 0) scale_factors(v) = 1.0d0

            ierr = nf90_get_att(met_ncids(v), varid, 'add_offset', add_offsets(v))
            if (ierr .ne. 0) add_offsets(v) = 0.0d0

            !get the met data	     
            starti(1) = 1
            starti(2) = gtoget
            counti(1) = timelen_spinup(v)
            counti(2) = 1

            if (i == 1 .and. v == 1)  then 
              allocate(atm_input       (met_nvars,begg:endg,1,1:counti(1)))
            end if 

            ierr = nf90_get_var(met_ncids(v), varid, atm_input(v,g:g,1,1:counti(1)), starti(1:2), counti(1:2))
            ierr = nf90_close(met_ncids(v))
          end do

    end do

    t2 = MPI_Wtime()
    write (100+myid, *) 't2: ', t2-t1,' @rank ', myid

    deallocate(xc)
    deallocate(yc)
    deallocate(lonc)
    deallocate(latc)
    deallocate(gindx)

    deallocate(timelen)
    deallocate(timelen_spinup)
    deallocate(tindex)
    deallocate(npf)
    deallocate(add_offsets)
    deallocate(scale_factors)
    deallocate(timeres)
    deallocate(var_offset)
    deallocate(var_mult)
    deallocate(atm_input)

    call MPI_Finalize(ierr)

  !===============================================================================

end program cplbypass_metload



