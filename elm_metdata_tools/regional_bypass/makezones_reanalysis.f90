program makezones


!Take the CLM CRU-NCEP data and aggregate by time over the spinup period (1901-1920)
!  and "chop" it up into 24 longitundinal zones

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
character(len=150) metvars, myforcing, mydomain, forcdir, myres
character(len=300) fname, filename_base
character(len=200) met_prefix
character(len=20) met_type

integer ierr, ncid, varid, dimid1
! met data
integer, parameter:: nzx = 24, nzgx = 30 ! max. zones, max. longitudal grids in a N-S zone
integer ng, nz, nzg                      ! actual grid numbers, zone numbers, E-W grids in a N-S zone

real, pointer :: data_in(:,:,:)      !248
integer, pointer :: data_zone(:,:)   !2920
integer, pointer :: temp_zone(:,:)   !248
real, pointer ::  longxy(:,:), latixy(:,:)
real, pointer ::  longxy_out(:,:), latixy_out(:,:)
integer, pointer ::  count_zone(:), ncid_out(:)
integer, pointer :: varids_out(:,:)
integer starti(3), counti(3), dimid(2), starti_out, starti_out_year



! domain/mask
integer ni, nj, nv
integer, pointer:: mask(:,:)
double precision, pointer :: xc(:,:), xc1d(:), yc(:,:), yc1d(:), xv(:,:,:), yv(:,:,:)
double precision :: xmin, xmax, ymin, ymax

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

met_type = 'GSWP3'! 'CRUJRA' ! 'CRU-NCEP'

!Set the input data path
!forcdir = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/ACME_inputdata/atm/datm7/atm_forcing.datm7.CRUJRA.0.5d.v1.c190604'
forcdir = '/lustre/or-hydra/cades-ccsi/proj-shared/project_acme/cesminput_ngee/atm/datm7' // &
           '/atm_forcing.datm7.GSWP3_daymet.56x24pt_kougarok-GRID'
!forcdir = '/Users/f9y/clm4_5_inputdata/atm/datm7' // &
!            '/atm_forcing.datm7.GSWP3_daymet.56x24pt_kougarok-GRID'

!Set the file header of the forcing, excluding 'clmforc.'
!myforcing = 'cruncep.V8.c2017'
!myforcing = 'CRUJRAV1.1.c2019.0.5x0.5'
myforcing = 'GSWP3_daymet.56x24pt_kougarok-GRID'


! domain/mask
mydomain = 'domain.lnd.GSWP3_daymet.56x24pt_kougarok-GRID'

!Set the date range and time resolution
!startyear = 1901
!endyear   = 2017
!res       = 6      !Timestep in hours
startyear = 1980
endyear   = 2010
res       = 3      !Timestep in hours

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

!if (myid == 0) then ! need to read in data for each rank, when np=7, i.e for individual of 7 vars on each rank
    fname = trim(forcdir) // '/' // trim(mydomain) // '.nc'
    ierr = nf90_open(trim(fname), NF90_NOWRITE, ncid)
    if (myid==0) print*, 'domain: ', trim(fname), ierr
    ierr = nf90_inq_dimid(ncid, "ni", dimid1)
    ierr = nf90_inquire_dimension(ncid, dimid1, len = ni)
    if (ierr /= 0) then
       print *, 'error in dimension ni', ierr
       stop
    endif
    ierr = nf90_inq_dimid(ncid, "nj", dimid1)
    ierr = nf90_inquire_dimension(ncid, dimid1, len = nj)
    if (ierr /= 0) then
       print *, 'error in dimension nj', ierr
       stop
    endif
    ierr = nf90_inq_dimid(ncid, "nv", dimid1)
    ierr = nf90_inquire_dimension(ncid, dimid1, len = nv)
    if (ierr /= 0) then
       print *, 'error in dimension nv', ierr
       stop
    endif

    ng = ni*nj
    if (myid==0) print *, 'ni, nj, nj, ng', ni, nj, nv, ng

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
    xmin = minval(xc1d)
    xmax = maxval(xc1d)
    yc1d = yc(1,:)
    ymin = minval(yc1d)
    ymax = maxval(yc1d)
    if (myid==0) print *, 'range of xc: ', xmin, xmax
    if (myid==0) print *, 'range of yc: ', ymin, ymax

!endif


print*, 'myid', myid
if (myid .eq. 0) open(unit=8, file=trim(forcdir) // '/cpl_bypass_full/zone_mappings.txt')

nz  = min(nzx, (ni-1)/nzgx+1)
nzg = min(nzgx, (ni-1)/nz+1)
print *, 'zones, zone grids', nz, nzg

nty = (NHRSY-1)/res+1               ! ntime in a year
ntd = 23/res+1                      ! ntime in a day
ntm = ntd*31                        ! ntime in a month(max. days)
allocate(longxy(ni,nj))
allocate(latixy(ni,nj))
allocate(data_in(nzg, nj, ntm))     ! monthly data read-in for a N-S zone
allocate(data_zone(nty,nzg*nj))     ! yearly data for a N-S zone (1D grid)
allocate(temp_zone(ntm,nzg*nj))     ! monthly data array for a zone (1D grid)
allocate(longxy_out(nz,15000))      ! ? 15000
allocate(latixy_out(nz,15000))
allocate(count_zone(nz))
allocate(ncid_out(nz))
allocate(varids_out(nz,10))         ! 10 variables at the most
allocate(dtime(nty))

do v=myid+1,7,np
 print*, v
 do z=1,nz
   mask(:,:)=0
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
   if(trim(met_type) == 'GSWP3') then
     if (v .eq. 1) met_prefix = trim(forcdir) // '/Precip' // trim(myres) // &
         '/clmforc.' // trim(myforcing) // '.Prec.'
     if (v .eq. 2) met_prefix = trim(forcdir) // '/Solar' // trim(myres) // &
         '/clmforc.' // trim(myforcing) // '.Solr.'
     if (v .ge. 3) met_prefix = trim(forcdir) // '/TPHWL' // trim(myres) // &
         '/clmforc.' // trim(myforcing) // '.TPQWL.'

   elseif(trim(met_type) == 'CRUJRA') then
   ! CRUJRA
     if (v .eq. 1) met_prefix = trim(forcdir) // '/clmforc.' // trim(myforcing) // '.Prec.'
     if (v .eq. 2) met_prefix = trim(forcdir) // '/clmforc.' // trim(myforcing) // '.Solr.'
     if (v .ge. 3) met_prefix = trim(forcdir) // '/clmforc.' // trim(myforcing) // '.TPQWL.'

   elseif(trim(met_type) == 'CRUJRA') then

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
         print*, trim(metvars), y, m
         count_zone(:)=0
         write(yst,'(I4)') y
         fname = trim(met_prefix) // yst // '-' // mst(3:4) // '.nc'
         ierr = nf90_open(trim(fname), NF90_NOWRITE, ncid)
         print*, trim(fname), ierr, v
         ierr = nf90_inq_varid(ncid, 'LONGXY', varid)
         ierr = nf90_get_var(ncid, varid, longxy)
         ierr = nf90_inq_varid(ncid, 'LATIXY', varid)
         ierr = nf90_get_var(ncid, varid, latixy)
         ierr = nf90_inq_varid(ncid, trim(metvars), varid)
         starti(1:3)  = 1
         starti(1)    = (z-1)*nzg+1
         counti(1)  = nzg
         counti(2)  = nj
         counti(3)  = ndaysm(m)*(ntd)
   
         ierr = nf90_get_var(ncid, varid, data_in(1:counti(1),1:counti(2),1:counti(3)), starti, counti)
         print*, 'READ', ierr
         ierr = nf90_close(ncid)
         do i=1,nzg
            do j=1,nj
               if (data_in(i,j,1) .le. 1e9) mask(i,j)=1
               if (mask(i,j) == 1) then 

                  count_zone(z)=count_zone(z)+1
                  if (y .eq. startyear .and. m .eq. 1) then

                     if (myid .eq. 0) then 
                        !if (longxy((z-1)*nzg+i,j) .ge. 0.25) then
                           write(8,'(f12.5,1x,f12.6,1x,2(I5,1x))') longxy((z-1)*nzg+i,j), latixy((z-1)*nzg+i,j), &
                                z, count_zone(z)
                        !else
                        !   write(8,'(2(f10.3,1x),2(I5,1x))') xmax, latixy((z-1)*nzg+i,j), &
                        !        z, count_zone(z)
                        !end if
                     end if
                     longxy_out(z, count_zone(z)) = longxy((z-1)*nzg+i,j)
                     latixy_out(z, count_zone(z)) = latixy((z-1)*nzg+i,j)
                  end if
                  !print*, i,j,z, count_zone(z) 
                  temp_zone(1:ndaysm(m)*ntd,count_zone(z)) = &
                       nint((data_in(i,j,1:ndaysm(m)*ntd)-add_offsets(v))/scale_factors(v))
               end if
            end do
         end do
         print*, 'ZONE', z, count_zone(z)
         !if (y .eq. startyear .and. m .eq. 1 .and. myid .eq. 0) close(8)  
         
         
         !do z=mod(myid,24)+1,24,np
            write(zst,'(I4)') 1000+z
            if (y .eq. startyear .and. m .eq. 1) then 
               fname = trim(forcdir) // '/cpl_bypass_full/' // trim(myforcing) &
                      // '_' // trim(metvars) // '_' // startyrst // '-' // endyrst // '_z' // &
                    zst(3:4) // '.nc'
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
            !ierr = nf90_put_var(ncid_out(z), varids_out(z,4), &
            !     temp_zone(1:counti(1), 1:count_zone(z)), starti(1:2), counti(1:2))
            !print*,'WRITEVAR', ierr
         !end do  !Zone loop
         starti_out_year = starti_out_year+ndaysm(m)*(24/res)
         starti_out = starti_out+ndaysm(m)*(24/res)
         !print*, starti_out
      end do    !month loop
      starti(1) = (y-startyear)*(8760/res)+1
      starti(2) = 1
      counti(1) = 8760/res
      counti(2) = count_zone(z)
      print*, y, z, v, starti(1:2)
      print*, z, counti(1:2)
      ierr = nf90_put_var(ncid_out(z), varids_out(z,4), &
                   data_zone(1:counti(1), 1:counti(2)), starti(1:2), counti(1:2))
   end do       !year loop
   ierr = nf90_close(ncid_out(z))
 end do  !zone loop
 if (myid .eq. 0) close(8)  
end do   !Variable loop

deallocate(xc)
deallocate(xc1d)
deallocate(yc)
deallocate(yc1d)
deallocate(xv)
deallocate(yv)
deallocate(mask)

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

call MPI_Finalize(ierr)

end program makezones
