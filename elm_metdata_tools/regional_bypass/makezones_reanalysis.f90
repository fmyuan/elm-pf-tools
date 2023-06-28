program makezones


use netcdf
implicit none
include 'mpif.h'
!include 'netcdf.h'

integer myid, np  ! NOTE: np must be <=7, with 1 variable on 1 rank ideally

integer res
integer v, n, i, z, y,m,j
character(len=4) yst, startyrst, endyrst
character(len=4) mst, myidst
character(len=4) zst
character(len=1) rst
character(len=150) metvars, myforcing, mydomain, forcdir, myres, odir
character(len=300) fname, filename_base
character(len=200) met_prefix
character(len=20) met_type

integer ierr, ncid, varid, dimid1
! met data zoning parameters
integer, parameter:: ngzx = 50000           ! max. grids in a N-S strip (i.e. zone)
integer nzx, nz, ng                         ! min. and actual numbers of zones (N-S strips), total masked grids
integer, pointer :: ngi(:), zng(:)          ! actual masked grids in each x/longitude interval and each zone
integer nzg                                 ! max. zone i-axis intervals
integer zoffset, zoffset_n                  ! when processing multiple Daymet Tiles, it's setting an offset so that zone numbering/counting is continuous
integer xindx, yindx                        ! 2-D daymet x/y index

real, pointer :: data_in(:,:,:)      !248
integer, pointer :: data_zone(:,:)   !2920
integer, pointer :: temp_zone(:,:)   !248
logical :: NEW_NC
real, pointer ::  longxy(:,:), latixy(:,:)
real, pointer ::  longxy_out(:,:), latixy_out(:,:)
integer, pointer ::  count_zone(:), ncid_out(:)
integer, pointer :: varids_out(:,:)
integer starti(3), counti(3), dimid(2), starti_out, starti_out_year



! domain/mask
integer ni, nj, nv
integer, pointer:: mask(:,:)
integer, pointer:: zstarti(:), zcounti(:)
character(len=150) fstmetfile, mask_varname
double precision, pointer :: xc(:,:), xc1d(:), yc(:,:), yc1d(:), xv(:,:,:), yv(:,:,:)
double precision :: xmin, xmax, ymin, ymax
character(len=20) site
logical :: UNSTRUCTURED

integer nty, ntm, ntd
integer startyear, endyear
double precision, pointer :: dtime(:)    ! timing in a year (1460 for 6-hourly)

real add_offsets(7), scale_factors(7), data_ranges(7,2)
integer ndaysm(12)
data ndaysm /31,28,31,30,31,30,31,31,30,31,30,31/
integer, parameter:: NHRSY = 8760        ! hours in a year


!--------------------------------------------------------------------------------


call MPI_init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

! -------------------------------------------------------------
! USER-SPECIFIED forcing data input directory, file header, domain, and period

met_type = 'ESM_daymet4' !'ESM_daymet4' !'GSWP3_daymet4'!'GSWP3'! 'CRUJRA' ! 'CRU-NCEP'

!Set the input data path
!forcdir = '/lustre/or-scratch/cades-ccsi/proj-shared/project_acme/e3sm_inputdata/atm/datm7/CRUJRAv2'

!forcdir = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/atm/datm7' // &
!           '/atm_forcing.datm7.GSWP3_daymet.56x24pt_kougarok-GRID'
!forcdir = '/nfs/data/ccsi/f9y/GSWP3_daymet' // &
!            '/TILE13867'

!forcdir = '/gpfs/alpine/proj-shared/cli146/GSWP3_daymet' // &
!            '/TILE10835'


forcdir = '/gpfs/alpine/proj-shared/cli144/5v1/CMIP6_City' // &
          '/Knoxville' // &
          '/ACCESS-CM2_ssp585_r1i1p1f1'
!          '/ACCESS-CM2_ssp585_r1i1p1f1_DBCCA_Daymet_2020_2059'
!          '/ACCESS-CM2_ssp585_r1i1p1f1_DBCCA_Daymet_2060_2099'
!          '/BCC-CSM2-MR_ssp585_r1i1p1f1_DBCCA_Daymet_1980_2019'
!          '/BCC-CSM2-MR_ssp585_r1i1p1f1_DBCCA_Daymet_2020_2059'
!          '/BCC-CSM2-MR_ssp585_r1i1p1f1_DBCCA_Daymet_2060_2099'
!          '/CNRM-ESM2-1_ssp585_r1i1p1f2_DBCCA_Daymet_1980_2019'
!          '/CNRM-ESM2-1_ssp585_r1i1p1f2_DBCCA_Daymet_2020_2059'
!          '/CNRM-ESM2-1_ssp585_r1i1p1f2_DBCCA_Daymet_2060_2099'
!          '/MPI-ESM1-2-HR_ssp585_r1i1p1f1_DBCCA_Daymet_1980_2019'
!          '/MPI-ESM1-2-HR_ssp585_r1i1p1f1_DBCCA_Daymet_2020_2059'
!          '/MPI-ESM1-2-HR_ssp585_r1i1p1f1_DBCCA_Daymet_2060_2099'
!          '/MRI-ESM2-0_ssp585_r1i1p1f1_DBCCA_Daymet_1980_2019'
!          '/MRI-ESM2-0_ssp585_r1i1p1f1_DBCCA_Daymet_2020_2059'
!          '/MRI-ESM2-0_ssp585_r1i1p1f1_DBCCA_Daymet_2060_2099'
!          '/NorESM2-MM_ssp585_r1i1p1f1_DBCCA_Daymet_1980_2019'
!          '/NorESM2-MM_ssp585_r1i1p1f1_DBCCA_Daymet_2020_2059'
!          '/NorESM2-MM_ssp585_r1i1p1f1_DBCCA_Daymet_2060_2099'


!
odir = './'
zoffset = 0   ! zone index counting previously generated data set (e.g. DAYMET tile already processed)
zoffset_n = 0 ! zone n-index counting previously generated data set

!Set the file header of the forcing, excluding 'clmforc.'
!myforcing = 'cruncep.V8.c2017'
!myforcing = 'CRUJRAV2.3.c2023.0.5x0.5'
!myforcing = 'CRUJRAV1.1.c2019.0.5x0.5'
!myforcing = 'GSWP3_daymet.56x24pt_kougarok-GRID.'
!myforcing = 'Daymet4.1km.'
myforcing = ''

! domain/mask
UNSTRUCTURED = .false.
!mydomain = 'domain.lnd.GSWP3_daymet.56x24pt_kougarok-GRID'
mydomain = ' '
!fstmetfile = 'TPHWL3Hrly/clmforc.Daymet4.1km.TPQWL.1980-01.nc' ! have to get info on grid no. and their centroids
!fstmetfile = 'TPHWL6Hrly/clmforc.CRUJRAV2.3.c2023.0.5x0.5.TPQWL.1901-01.nc' ! have to get info on grid no. and their centroids
fstmetfile = 'TPHWL1Hrly/clmforc.TPQWL.1980-01.nc' ! have to get info on grid no. and their centroids
mask_varname = 'TBOT'

! if setting the following as non-'-9999.0' value, mask will cut off by them
site = ''
xmin = -9999.0  ! xmin must be in 0~360 format
xmax = -9999.0  ! xmax must be in 0~360 format
ymin = -9999.0
ymax = -9999.0
! AK SewPenn NGEE sites: 
!   (1) Kougarok Road m64: -164.845~-164.800/65.175~65.149
!site = '-KM64'
!xmin = -164.845+360.0
!xmax = -164.800+360.0
!ymin = 65.149
!ymax = 65.175

!xmin = -110.0628+360.0
!xmax = -110.0428+360.0
!ymin = 31.7305
!ymax = 31.7505

!Set the date range and time resolution
!startyear = 1980   ! defaut 1980
!endyear   = 2014   ! defaut 2014
!res       = 3      ! Timestep in hours

!startyear = 1901
!endyear   = 2017
!res       = 6      !Timestep in hours
startyear = 1980   ! defaut 1980
endyear   = 2099   ! defaut 2014
res       = 1      ! Timestep in hours

!!-----------------------------------------------------------------

data_ranges(1,1) =-0.04
data_ranges(1,2) = 0.04   !Precip
data_ranges(2,1) = -20.
data_ranges(2,2) = 2000.  !FSDS
data_ranges(3,1) = 175.
data_ranges(3,2) = 350.   !TBOT
data_ranges(4,1) = 0.
data_ranges(4,2) = 0.10   !QBOT
data_ranges(5,1) = 0.
data_ranges(5,2) = 1000.  !FLDS
data_ranges(6,1) = 20000.
data_ranges(6,2) = 120000.  !PBOT
data_ranges(7,1) = -1.
data_ranges(7,2) = 100.     !WIND

do v=1,7
   add_offsets(v) = (data_ranges(v,2)+data_ranges(v,1))/2.
   scale_factors(v) = (data_ranges(v,2)-data_ranges(v,1))*1.1/2**15
end do

if (trim(mydomain) /= '') then
    fname = trim(forcdir) // '/' // trim(mydomain) // '.nc'
    ierr = nf90_open(trim(fname), NF90_NOWRITE, ncid)
    if (myid==0) print*, 'domain: ', trim(fname), ' - ', trim(nf90_strerror(ierr))
    ierr = nf90_inq_dimid(ncid, "ni", dimid1)
    ierr = nf90_inquire_dimension(ncid, dimid1, len = ni)
    if (ierr /= 0) then
       print *, 'error in dimension ni - ', trim(nf90_strerror(ierr))
       stop
    endif
    ierr = nf90_inq_dimid(ncid, "nj", dimid1)
    ierr = nf90_inquire_dimension(ncid, dimid1, len = nj)
    if (ierr /= 0) then
       print *, 'error in dimension nj - ', trim(nf90_strerror(ierr))
       stop
    endif
    ierr = nf90_inq_dimid(ncid, "nv", dimid1)
    ierr = nf90_inquire_dimension(ncid, dimid1, len = nv)
    if (ierr /= 0) then
       print *, 'error in dimension nv ', trim(nf90_strerror(ierr))
       stop
    endif

    if (myid==0) print *, 'ni, nj, nv, ng', ni, nj, nv, ni*nj

    allocate(xc(ni,nj))
    allocate(yc(ni,nj))
    allocate(xc1d(ni))
    allocate(yc1d(nj))
    allocate(mask(ni,nj))
    allocate(xv(ni,nj,nv))
    allocate(yv(ni,nj,nv))
    ierr = nf90_inq_varid(ncid, 'xc', varid)
    ierr = nf90_get_var(ncid, varid, xc)
    ierr = nf90_inq_varid(ncid, 'yc', varid)
    ierr = nf90_get_var(ncid, varid, yc)
    ierr = nf90_inq_varid(ncid, 'xv', varid)
    ierr = nf90_get_var(ncid, varid, xv)
    ierr = nf90_inq_varid(ncid, 'yv', varid)
    ierr = nf90_get_var(ncid, varid, yv)
    ierr = nf90_inq_varid(ncid, 'mask', varid)
    ierr = nf90_get_var(ncid, varid, mask)
    if(ierr /= 0) mask(:,:) = 1
    ierr = nf90_close(ncid)

    xc1d = xc(:,1)
    yc1d = yc(1,:)

! NO domain file, so have to build one from the first metdata
elseif(trim(fstmetfile) /= '') then

    fname = trim(forcdir) // '/' // trim(fstmetfile)
    ierr = nf90_open(trim(fname), NF90_NOWRITE, ncid)
    if (myid==0) print*, 'domain info from: ', trim(fname),  ' - ',  trim(nf90_strerror(ierr))

    ! dim: ni
    ierr = nf90_inq_dimid(ncid, "x", dimid1)
    if (ierr/=0) ierr = nf90_inq_dimid(ncid, "lon", dimid1)
    if (ierr/=0) then
       ierr = nf90_inq_dimid(ncid, "grd", dimid1) ! 1-D unstructured
       if (ierr==0) UNSTRUCTURED=.true.
    endif
    if (ierr/=0) then
       print *, 'STOPPED: CANNOT get info on dimension ni - ', trim(nf90_strerror(ierr))
       stop
       !
    else
    
       ierr = nf90_inquire_dimension(ncid, dimid1, len = ni)
       if (ierr /= 0) then
          print *, 'error in dimension ni -', trim(nf90_strerror(ierr))
          stop
       endif
    endif
    
    !dim: nj
    ierr = nf90_inq_dimid(ncid, "y", dimid1)
    if (ierr/=0) ierr = nf90_inq_dimid(ncid, "lat", dimid1)
    if (ierr==0) then
       ierr = nf90_inquire_dimension(ncid, dimid1, len = nj)
       if (ierr /= 0) then
          print *, 'error in dimension nj - ', trim(nf90_strerror(ierr))
          stop
       endif
    else
       ! in case that 'grd' as 1D unstructed dimension, nj=1
       ierr = nf90_inq_dimid(ncid, "grd", dimid1)
       if (ierr==0) then
          nj = 1
       else
          print *, 'CANNOT get info on dimension nj - ', trim(nf90_strerror(ierr))
          stop
       endif

    endif

    if (myid==0) print *, ' domain info - ni, nj, ng', ni, nj, ni*nj


    allocate(xc(ni,nj))
    allocate(yc(ni,nj))
    allocate(xc1d(ni))
    if (UNSTRUCTURED) then
      allocate(yc1d(ni))
    else
      allocate(yc1d(nj))
    endif
    allocate(mask(ni,nj))

    if (trim(met_type(1:12)) == 'GSWP3_daymet') then
      ierr = nf90_inq_varid(ncid, 'lon', varid)
      ierr = nf90_get_var(ncid, varid, xc)
      ierr = nf90_inq_varid(ncid, 'lat', varid)
      ierr = nf90_get_var(ncid, varid, yc)

      ierr = nf90_inq_varid(ncid, 'x', varid)  ! daymet geox axis in meters
      ierr = nf90_get_var(ncid, varid, xc1d)

      ierr = nf90_inq_varid(ncid, 'y', varid)  ! daymet geoy axis in meters
      ierr = nf90_get_var(ncid, varid, yc1d)
    else

      if (UNSTRUCTURED) then
        ierr = nf90_inq_varid(ncid, 'lon', varid)
        ierr = nf90_get_var(ncid, varid, xc)
        ierr = nf90_inq_varid(ncid, 'lat', varid)
        ierr = nf90_get_var(ncid, varid, yc)
        ierr = nf90_inq_varid(ncid, 'x', varid)  ! daymet geox axis in meters
        ierr = nf90_get_var(ncid, varid, xc1d)
        ierr = nf90_inq_varid(ncid, 'y', varid)  ! daymet geoy axis in meters
        ierr = nf90_get_var(ncid, varid, yc1d)

      else
        ! 2-D grids
        ierr = nf90_inq_varid(ncid, 'lon', varid)  ! lon axis in degree
        if (ierr==0) then
          ! lon/lat as grid mesh axis (middle)
          ierr = nf90_get_var(ncid, varid, xc1d)
          ierr = nf90_inq_varid(ncid, 'lat', varid)  ! lat axis in degree
          ierr = nf90_get_var(ncid, varid, yc1d)
      
          do j=1, nj
            do i=1,ni
              xc(i,j) = xc1d(i)
              yc(i,j) = yc1d(j)
            end do
          end do

        else
          ! don't have axis but full list of 2D grids' locations
          ierr = nf90_inq_varid(ncid, 'LONGXY', varid)  ! centroid of a grid
          ierr = nf90_get_var(ncid, varid, xc)
          ierr = nf90_inq_varid(ncid, 'LATIXY', varid)  ! centroid of a grid
          ierr = nf90_get_var(ncid, varid, yc)
          xc1d = xc(:,1)
          yc1d = yc(1,:)
        end if

      endif

    end if

    ! mask based on if metvar values are filled or NaN
    ierr = nf90_inq_varid(ncid, trim(mask_varname), varid)
    starti(1:3)  = 1
    counti(1)    = ni
    counti(2)    = nj
    counti(3)    = 1 ! only needs 1 count of time
    allocate(data_in(ni,nj,1))
    ierr = nf90_get_var(ncid, varid, data_in, starti, counti)
    !data_in = merge(1,0,data_in<=1e9)
    !mask = reshape(data_in, shape(mask))
    ! the above seems ok with gnu 9, but not with gnu 6 on CADES,
    ! so the following is good for all compilers
    if (xmin/=-9999.0 .and. xmin<0.0) xmin=xmin+360.0
    if (xmax/=-9999.0 .and. xmax<0.0) xmax=xmax+360.0
    do i=1,ni
      do j=1,nj
        if (xc(i,j)/=-9999.0 .and. xc(i,j)<0.0) xc(i,j)=xc(i,j)+360.0

        if ( (data_in(i,j,1) .le. 1e9) .or. &
             (data_in(i,j,1) .gt. huge(data_in(i,j,1))) .or. &
             (.not.(data_in(i,j,1) /= data_in(i,j,1)) ) ) then
          mask(i,j) = 1

          ! cutoff mask, if boxed x/y ranges given
          if(xmin/=-9999.0) then
            if(xc(i,j)<0.0 .and. (xc(i,j)+360.0)<xmin) then ! xmin must be in 0~360 format
              mask(i,j) = 0
            elseif (xc(i,j)>=0.0 .and. xc(i,j)<xmin) then
              mask(i,j) = 0
            endif
          endif
          if(xmax/=-9999.0 .and. mask(i,j)==1) then
            if(xc(i,j)<0.0 .and. (xc(i,j)+360.0)>xmax) then ! xmax must be in 0~360 format
              mask(i,j) = 0
            elseif (xc(i,j)>=0.0 .and. xc(i,j)>xmax) then
              mask(i,j) = 0
            endif
          endif
          if(ymin/=-9999.0 .and. mask(i,j)==1) then
            if (yc(i,j)<ymin) mask(i,j) = 0
          endif
          if(ymax/=-9999.0 .and. mask(i,j)==1) then
            if (yc(i,j)>ymax) mask(i,j) = 0
          endif

        else
          mask(i,j) = 0
        endif
      enddo
    enddo
    deallocate(data_in)


    ierr = nf90_close(ncid)

endif

if (myid .eq. 0) open(unit=8, file=trim(odir) // '/cpl_bypass_full' // trim(site) // '/' // 'zone_mappings.txt')
if (myid .eq. 0 .and. index(trim(met_type),'daymet') .gt. 0) then
  open(unit=9, file=trim(odir) // '/cpl_bypass_full' // trim(site) // '/' // 'daymet_elm_mappings.txt')
  write(9,*) '     lon          lat           geox            geoy        z     i     j     g '
else
  open(unit=9, file=trim(odir) // '/cpl_bypass_full' // trim(site) // '/' // 'cplbypass_mappings.txt')
  write(9,*) '     lon          lat        grid-x       grid-y       z     i     j     g '

endif
! possible zone number and its grids
! the following is for regular zones - may have remarkable variations of masked grids among zones
!nzx  = min(1, (ni-1)/ngzx+1)
!nzg = min(ngzx, (ni-1)/nzx+1)
!nzr = ni-nzx*nzg  ! should be non-positive integer
allocate(ngi(ni))
ngi = sum(mask,DIM=2)
ng  = sum(ngi)
nzx = (ng-1)/(ngzx-nj)+1
nzg = min(ngzx-nj, (ng-1)/nzx+1)       ! this is temporary (average) number (integer)

allocate(zng(nzx+1))                   ! an extra length just in case
allocate(zstarti(nzx+1))               ! zone starting index along x/lon dimension
allocate(zcounti(nzx+1))               ! zone grid count number along x/lon dimension
zng(:)     = 0
zstarti(:) = 0
zcounti(:) = 0

! do zone division
nz  = 1
zstarti(1) = 1
do i=1,ni-1
  zcounti(nz) = zcounti(nz) + 1
  zng(nz)     = zng(nz)+ngi(i)
  if ((zng(nz)+ngi(i+1))>(nzg+nj/2)) then ! a new zone, if adding next count of i-index over (nzg+nj/2)
    if (myid==0) then
      print *, 'zone, starti, endi, counti, grids: ', &
        nz, zstarti(nz), zstarti(nz)+zcounti(nz)-1, zcounti(nz),zng(nz)
    endif
    nz = nz + 1
    zstarti(nz) = i+1
  endif
  if(i==ni-1) then  ! do-loop ends at last second i-index
    zng(nz)=zng(nz)+ngi(i+1)
    zcounti(nz)=zcounti(nz)+1
    if (myid==0) then
      print *, 'zone, starti, endi, counti, grids: ', &
        nz, zstarti(nz), zstarti(nz)+zcounti(nz)-1, zcounti(nz),zng(nz)
    endif
  endif
enddo
nzg = maxval(zcounti) ! the max i count for all zones (for max. ni when read in data_in(nzg, nj, ntm)
if (myid==0) then
  print *, 'max. counts of ni among all zones, zone number, total masked grids - ', &
    nzg, nz, ng
endif

!
nty = (NHRSY-1)/res+1               ! ntime in a year
ntd = 23/res+1                      ! ntime in a day
ntm = ntd*31                        ! ntime in a month(max. days)
allocate(longxy(ni,nj))
allocate(latixy(ni,nj))
allocate(data_in(nzg, nj, ntm))     ! monthly data read-in for a N-S zone
allocate(data_zone(nty,nzg*nj))     ! yearly data for a N-S zone (1D grid)
allocate(temp_zone(ntm,nzg*nj))     ! monthly data array for a zone (1D grid)
allocate(longxy_out(nz,ngzx))      ! ? 15000
allocate(latixy_out(nz,ngzx))
allocate(count_zone(nz))
allocate(ncid_out(nz))
allocate(varids_out(nz,10))         ! 10 variables at the most
allocate(dtime(nty))

do v=myid+1,7,np
 ng = 0  ! index of all (land-)masked grids
 do z=1,nz
   NEW_NC = .True.
   ! the following may change masked gridcell from variable/time to variable/time
   ! when 'mask(i,j)' upon variable's filling values later
   !mask(:,:)=0
   if (v .eq. 1) metvars='PRECTmms'
   if (v .eq. 2) metvars='FSDS'
   if (v .eq. 3) metvars='TBOT'
   if (v .eq. 4) metvars='QBOT'
   if (v .eq. 5) metvars='FLDS'
   if (v .eq. 6) metvars='PSRF'
   if (v .eq. 7) metvars='WIND'

   write(rst,'(I1)') res
   myres = trim(rst) // 'Hrly'
 
   ! GSWP3
   if(trim(met_type) == 'GSWP3' .or. trim(met_type(1:12)) == 'GSWP3_daymet') then
     if (v .eq. 1) met_prefix = trim(forcdir) // '/Precip' // trim(myres) // &
         '/clmforc.' // trim(myforcing) // 'Prec.'

     if (v .eq. 2) met_prefix = trim(forcdir) // '/Solar' // trim(myres) // &
         '/clmforc.' // trim(myforcing) // 'Solr.'
     if (v .ge. 3) met_prefix = trim(forcdir) // '/TPHWL' // trim(myres) // &
         '/clmforc.' // trim(myforcing) // 'TPQWL.'

   elseif(trim(met_type) == 'ESM_daymet4') then
     if (v .eq. 1) met_prefix = trim(forcdir) // '/Precip' // trim(myres) // &
         '/clmforc.' // trim(myforcing) // 'PRECTmms.'

     if (v .eq. 2) met_prefix = trim(forcdir) // '/Solar' // trim(myres) // &
         '/clmforc.' // trim(myforcing) // 'Solr.'
     if (v .ge. 3) met_prefix = trim(forcdir) // '/TPHWL' // trim(myres) // &
         '/clmforc.' // trim(myforcing) // 'TPQWL.'


   elseif(trim(met_type) == 'CRUJRA') then
   ! CRUJRA
     if (v .eq. 1) met_prefix = trim(forcdir)  // '/Precip' // trim(myres) // &
         '/clmforc.' // trim(myforcing) // '.Prec.'
     if (v .eq. 2) met_prefix = trim(forcdir) // '/Solar' // trim(myres) // &
         '/clmforc.' // trim(myforcing) // '.Solr.'
     if (v .ge. 3) met_prefix = trim(forcdir) // '/TPHWL' // trim(myres) // &
         '/clmforc.' // trim(myforcing) // '.TPQWL.'

   elseif(trim(met_type) == 'CRU-NCEP') then

   !CRU-NCEP v7
     if (v .eq. 1) met_prefix = trim(forcdir) // '/atm_forcing.datm7.' &
            // trim(myforcing) // '_qianFill.0.5d.V7.c160715' // '/Precip' // trim(myres) // &
           '/clmforc.cruncep.V7.c2016.0.5d.Prec.'
     if (v .eq. 2) met_prefix = trim(forcdir) // '/atm_forcing.datm7.' &
            // trim(myforcing) // '_qianFill.0.5d.V7.c160715' // '/Solar' // &
            trim(myres) // '/clmforc.cruncep.V7.c2016.0.5d.Solr.'
     if (v .ge. 3) met_prefix = trim(forcdir) // '/atm_forcing.datm7.' &
            // trim(myforcing) // '_qianFill.0.5d.V7.c160715' // '/TPHWL' // &
            trim(myres) // '/clmforc.cruncep.V7.c2016.0.5d.TPQWL.'
   
   endif

   data_in(:,:,:)=1e36
   write(startyrst,'(I4)') startyear
   write(endyrst,'(I4)') endyear
   starti_out = 1

   do y=startyear,endyear
      starti_out_year = 1
      do m=1,12
         write(mst,'(I4)') 1000+m
         print*, ' '
         print*, trim(metvars), y, m, z, 'on RANK ', myid
         count_zone(:)=0
         write(yst,'(I4)') y
         fname = trim(met_prefix) // yst // '-' // mst(3:4) // '.nc'
         ierr = nf90_open(trim(fname), NF90_NOWRITE, ncid)
         if(ierr/=0) then
           print*, 'Error in open: ', trim(fname), ' ', trim(nf90_strerror(ierr))
           stop
         endif

         ! if known grid-centroid lon/lat (full 2D)
         if (trim(fstmetfile) /= '') then
           longxy = xc
           latixy = yc
         else
           ! if not yet, read in grid-centroid lon/lat (note: make sure it's full 2D)
           ierr = nf90_inq_varid(ncid, 'LONGXY', varid)
           if (ierr==0) then
             ierr = nf90_get_var(ncid, varid, longxy)

             ierr = nf90_inq_varid(ncid, 'LATIXY', varid)
             if (ierr/=0) then
               print *, "invalid variable 'LATIXY'",trim(nf90_strerror(ierr))
               stop
             else           
               ierr = nf90_get_var(ncid, varid, latixy)
             endif
           else
             ierr = nf90_inq_varid(ncid, 'lon', varid)
             ierr = nf90_get_var(ncid, varid, longxy)           

             if(ierr==0) then
               ierr = nf90_inq_varid(ncid, 'lat', varid)
               if (ierr/=0) then      
                 print *, "invalid variable 'lon' ",trim(nf90_strerror(ierr))
                 stop
               else 
                 ierr = nf90_get_var(ncid, varid, latixy)
               endif
             endif
           endif
         
         endif

         ! get monthly metdata
         ierr = nf90_inq_varid(ncid, trim(metvars), varid)
         starti(1:3) = 1
         starti(1)   = zstarti(z)
         counti(1)   = zcounti(z)
         counti(2)   = nj
         counti(3)   = ndaysm(m)*(ntd)
   
         ierr = nf90_get_var(ncid, varid, data_in(1:counti(1),1:counti(2),1:counti(3)), starti, counti)
         if(ierr /= 0) then
           print*, 'error: ', ierr,  ' - ', trim(nf90_strerror(ierr)), &
                   ' In reading file: ',fname, 'on rank:', myid

           stop
         else
           print*, 'Successfully READ file: ', fname
         endif
         ierr = nf90_close(ncid)

         do i=1,counti(1)

            do j=1,nj
               ! the following may change masked gridcell from variable/time to variable/time
               !if ( (data_in(i,j,1) .le. 1e9) .or. &
               !     (.not.(data_in(i,j,1) /= data_in(i,j,1))) ) mask(i,j)=1

               ! better to let 'longitude' in format of 0~360
               if (longxy(starti(1)+i-1,j)<0.0) then
                   longxy(starti(1)+i-1,j) = longxy(starti(1)+i-1,j) + 360.0
               endif

               if (mask(starti(1)+i-1,j) == 1) then

                  count_zone(z)=count_zone(z)+1
                  !if (y .eq. startyear .and. m .eq. 1) then
                  if (NEW_NC) then

                     ng = ng + 1
                     if (myid .eq. 0) then 

                           write(8,'(f12.5,1x,f12.6,1x,I5,1x,I9)') longxy(starti(1)+i-1,j), &
                                latixy(starti(1)+i-1,j), &
                                z+zoffset, count_zone(z)+zoffset_n  ! offset only for writing zone_mappings.txt

                           ! xindex, yindex for 2D <--> 1-D transformat
                           if (zoffset_n<=0) then
                               xindx = starti(1)+i-1
                               yindx = j
                           else
                               xindx = -9999  ! when joint tiles, the x/y index need to be re-done
                               yindx = -9999
                           endif
                           if (index(trim(met_type),'daymet') .gt. 0) then
                             if (UNSTRUCTURED) then
                               write(9,'(f12.6,1x,f12.6,1x, 2(f20.6, 1x),3(I5,1x), I9)') longxy(starti(1)+i-1,j), &
                                latixy(starti(1)+i-1,j),                  &
                                xc1d(starti(1)+i-1), yc1d(starti(1)+i-1), &
                                z+zoffset, xindx, yindx, ng+zoffset_n
                             else
                               write(9,'(f12.6,1x,f12.6,1x, 2(f20.6, 1x),3(I5,1x), I9)') longxy(starti(1)+i-1,j), &
                                latixy(starti(1)+i-1,j),              &
                                xc1d(starti(1)+i-1), yc1d(j),         &
                                z+zoffset, xindx, yindx, ng+zoffset_n
                             endif
                           endif
                     end if
                     longxy_out(z, count_zone(z)) = longxy(starti(1)+i-1,j)
                     latixy_out(z, count_zone(z)) = latixy(starti(1)+i-1,j)
                  end if

                  temp_zone(1:ndaysm(m)*ntd,count_zone(z)) = &
                       nint((data_in(i,j,1:ndaysm(m)*ntd)-add_offsets(v))/scale_factors(v))

               end if
            end do
         end do
         ! infos
         if (myid == 0 .and. NEW_NC .and. v==1) then
            print *, 'Total valid grids - ', count_zone(z), ' OF ', counti(1)*nj,' In zone - ', z
         endif

         
         !do z=mod(myid,24)+1,24,np
            write(zst,'(I4)') 1000+z+zoffset
            !if (y .eq. startyear .and. m .eq. 1) then
            if (NEW_NC) then
               if(trim(met_type(1:5)) == 'GSWP3' .or. trim(met_type) == 'ESM_daymet4') then
                 fname = trim(odir) // '/cpl_bypass_full' // trim(site) // '/' // trim(met_type) &
                      // '_' // trim(metvars) // '_' // startyrst // '-' // endyrst // '_z' // &
                    zst(3:4) // '.nc'
               else
                 fname = trim(odir) // '/cpl_bypass_full' // trim(site) // '/' // trim(myforcing) &
                      // '_' // trim(metvars) // '_' // startyrst // '-' // endyrst // '_z' // &
                    zst(3:4) // '.nc'
               endif

               ierr = nf90_create(trim(fname),cmode=or(nf90_clobber,nf90_64bit_offset),ncid=ncid_out(z))
               ierr = nf90_def_dim(ncid_out(z), 'n', count_zone(z), dimid(2))
               ierr = nf90_def_dim(ncid_out(z), 'DTIME', (endyear-startyear+1)*(NHRSY/res), dimid(1))
               ierr = nf90_def_var(ncid_out(z), 'DTIME', NF90_DOUBLE, dimid(1), &
                    varids_out(z,1))
               ierr = nf90_put_att(ncid_out(z), varids_out(z,1), 'long_name', &
                    'Day of Year')
               ierr = nf90_put_att(ncid_out(z), varids_out(z,1), 'units', &
                    'Days since ' // startyrst // '-01-01 00:00')
               ierr = nf90_def_var(ncid_out(z), 'LONGXY', NF90_FLOAT, dimid(2), &
                    varids_out(z,2))
               ierr = nf90_def_var(ncid_out(z), 'LATIXY', NF90_FLOAT, dimid(2), &
                    varids_out(z,3))
               ierr = nf90_def_var(ncid_out(z), trim(metvars), NF90_SHORT, &
                    dimid(1:2), varids_out(z,4))
               ierr = nf90_put_att(ncid_out(z), varids_out(z,4), 'add_offset', &
                    add_offsets(v))
               ierr = nf90_put_att(ncid_out(z), varids_out(z,4), 'scale_factor', &
                    scale_factors(v))
               ierr = nf90_enddef(ncid_out(z))
               ierr = nf90_put_var(ncid_out(z), varids_out(z,2), &
                    longxy_out(z, 1:count_zone(z)))
               ierr = nf90_put_var(ncid_out(z), varids_out(z,3), &
                    latixy_out(z, 1:count_zone(z)))

               NEW_NC = .False.
            end if
            do i=1,ndaysm(m)*(24/res)
               dtime(i) = (starti_out+i-1)/(24/res*1.0)-(res/24)*0.5
            end do
            starti(1) = starti_out
            counti(1) = ndaysm(m)*(24/res)
            ierr = nf90_put_var(ncid_out(z), varids_out(z,1), &
                 dtime(1:counti(1)), starti(1:1), counti(1:1))
            starti(2) = 1
            counti(2) = count_zone(z)

            data_zone(starti_out_year:(starti_out_year+counti(1)-1),starti(2):(starti(2) &
                 +counti(2)-1)) = temp_zone(1:counti(1),1:counti(2))
         !end do  !Zone loop
         starti_out_year = starti_out_year+ndaysm(m)*(24/res)
         starti_out = starti_out+ndaysm(m)*(24/res)
         !
      end do    !month loop

      starti(1) = (y-startyear)*(8760/res)+1
      starti(2) = 1
      counti(1) = 8760/res
      counti(2) = count_zone(z)
      ierr = nf90_put_var(ncid_out(z), varids_out(z,4), &
                   data_zone(1:counti(1), 1:counti(2)), starti(1:2), counti(1:2))
      if (ierr/=0) then
        print*, 'Variable writing status:', ierr,  ' - ', trim(nf90_strerror(ierr))
        print*, 'yr:',y, 'zone:',z, 'var:',trim(metvars), counti(1:2), starti(1:2)
        print*, 'file: ', trim(fname)
        stop
      else
        print*, 'Writing - ', trim(nf90_strerror(ierr)), ' : ',trim(metvars), ' ',trim(fname)
      endif

   end do       !year loop
   ierr = nf90_close(ncid_out(z))
 end do  !zone loop
 if (myid .eq. 0) close(8)  
 if (myid .eq. 0) close(9)
end do   !Variable loop

call mpi_barrier(ierr)

deallocate(xc)
deallocate(xc1d)
deallocate(yc)
deallocate(yc1d)
!if(associated(xv)) deallocate(xv)
!if(associated(yv)) deallocate(yv)
!if(associated(mask)) deallocate(mask)

deallocate(longxy)
deallocate(latixy)
deallocate(data_in)
deallocate(data_zone)
deallocate(temp_zone)
deallocate(longxy_out)
deallocate(latixy_out)
deallocate(count_zone)
deallocate(ncid_out)
deallocate(varids_out)
deallocate(dtime)
deallocate(ngi)
deallocate(zng)
deallocate(zstarti)
deallocate(zcounti)

call MPI_Finalize(ierr)

end program makezones
