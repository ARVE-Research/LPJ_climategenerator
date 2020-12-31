program calc_wetdays

! ifort -xHost -o calcwetdays calcwetdays.f90 -I/share1/netcdf/4.7.1-impi/include -L/share1/netcdf/4.7.1-impi/lib -lnetcdff -lnetcdf
! gfortran -o calcwetdays calcwetdays.f90 -I/home/public/easybuild/software/netCDF-Fortran/4.5.2-gompi-2020a/include -L/home/public/easybuild/software/netCDF-Fortran/4.5.2-gompi-2020a/lib -lnetcdff -lnetcdf

use netcdf
use iso_fortran_env, only : real32,real64

implicit none

integer, parameter :: sp = real32
integer, parameter :: dp = real64

real(sp), parameter, dimension(12) :: ndays_in_month = [ 31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31. ]

real(sp), parameter :: missing = -9999.

character(200) :: precfile
character(200) :: coeff_file
character(200) :: outfile

integer :: status
integer :: ifid
integer :: dimid
integer :: varid
integer :: ovarid
integer :: ofid
integer :: xlen
integer :: ylen
integer :: tlen
integer :: x,y,m,t

real(dp) :: exponent

real(sp), allocatable, dimension(:,:,:) :: k0
real(sp), allocatable, dimension(:,:,:) :: k1
real(sp), allocatable, dimension(:,:,:) :: prec
real(sp), allocatable, dimension(:,:,:) :: wet

real(sp), dimension(2) :: valid_range

real(sp) :: wetf

real(sp), parameter :: default_k0 = 0.0005
real(sp), parameter :: default_k1 = 0.5

real(dp), allocatable, dimension(:) :: lon
real(dp), allocatable, dimension(:) :: lat
real(dp), allocatable, dimension(:) :: time

real(dp), allocatable, dimension(:,:) :: lon_bnds
real(dp), allocatable, dimension(:,:) :: lat_bnds
real(dp), allocatable, dimension(:,:) :: time_bnds

!---------------------------------------------------------------------------------

call getarg(1,precfile)
call getarg(2,coeff_file)
call getarg(3,outfile)

write(0,*)'estimating wet days on the basis of total precip, month, and geographic location'

status = nf90_open(coeff_file,nf90_nowrite,ifid)
if (status /= nf90_noerr) call handle_err(status)

!--------------------------------------------------------------------------------- 

! inquire dimensions

status = nf90_inq_dimid(ifid,'lon',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=xlen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'lat',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=ylen)
if (status /= nf90_noerr) call handle_err(status)

!---------------------------------------------------------------------------------
!get the two coefficients k0 and k1 for each grid cell and month

write(0,*)'get coeffs',xlen,ylen

allocate(k0(xlen,ylen,12))
allocate(k1(xlen,ylen,12))

status = nf90_inq_varid(ifid,'k0',varid) !a
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,k0)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ifid,'k1',varid) !b
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,k1)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ifid)
if (status /= nf90_noerr) call handle_err(status)

!set k0 and k1 to default values for undefined

where (k0 == missing) k0 = default_k0
where (k1 == missing) k1 = default_k1

!---------------------------------------------------------------------------------

status = nf90_open(precfile,nf90_nowrite,ifid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_open(outfile,nf90_write,ofid)
if (status /= nf90_noerr) call handle_err(status)

!---------------------------------------------------------------------------------

status = nf90_inq_dimid(ifid,'lon',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=xlen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'lat',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=ylen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ifid,'time',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ifid,dimid,len=tlen)
if (status /= nf90_noerr) call handle_err(status)

allocate(lon(xlen))
allocate(lat(ylen))
allocate(time(tlen))

allocate(lon_bnds(2,xlen))
allocate(lat_bnds(2,ylen))
allocate(time_bnds(2,tlen))

!---------------------------------------------------------------------------------

status = nf90_inq_varid(ifid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,lon)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'lon',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,lon)
if (status /= nf90_noerr) call handle_err(status)

!---

status = nf90_inq_varid(ifid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,lat)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'lat',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,lat)
if (status /= nf90_noerr) call handle_err(status)

!---

status = nf90_inq_varid(ifid,'time',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,time)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'time',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,time)
if (status /= nf90_noerr) call handle_err(status)

!---

status = nf90_inq_varid(ifid,'lon_bnds',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,lon_bnds)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'lon_bnds',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,lon_bnds)
if (status /= nf90_noerr) call handle_err(status)

!---

status = nf90_inq_varid(ifid,'lat_bnds',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(ifid,varid,lat_bnds)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'lat_bnds',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,lat_bnds)
if (status /= nf90_noerr) call handle_err(status)

!---

! status = nf90_inq_varid(ifid,'time_bnds',varid)
! if (status /= nf90_noerr) call handle_err(status)
! 
! status = nf90_get_var(ifid,varid,time_bnds)
! if (status /= nf90_noerr) call handle_err(status)
! 
! status = nf90_inq_varid(ofid,'time_bnds',varid)
! if (status /= nf90_noerr) call handle_err(status)
! 
! status = nf90_put_var(ofid,varid,time_bnds)
! if (status /= nf90_noerr) call handle_err(status)

!---------------------------------------------------------------------------------
write(0,*)'get prec'

allocate(prec(xlen,ylen,12))
allocate(wet(xlen,ylen,12))

status = nf90_inq_varid(ifid,'prec',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'wetdays',ovarid)
if (status /= nf90_noerr) call handle_err(status)

do t = 1,tlen,12

  status = nf90_get_var(ifid,varid,prec,start=[1,1,t],count=[xlen,ylen,12])
  if (status /= nf90_noerr) call handle_err(status)
  
  !correct prec for negative values

  where (prec < 0.) prec = 0.

  !--------------------------------------------------------------------------------- 

  write(0,*)'calculating',t,xlen,ylen

  wet = missing 

  do x = 1, xlen
    do y = 1, ylen
      do m = 1, 12
    
        if (prec(x,y,m) >= 0.) then
      
          exponent = -1. * k0(x,y,m) * prec(x,y,m) 
     
          wetf = (1. - exp(exponent))**k1(x,y,m)
        
          wetf = max(min(wetf,1.),0.)
        
          wet(x,y,m) = wetf * ndays_in_month(m)
        
        else

          wet(x,y,m)  = missing

        end if

      end do
    end do
  end do
  
  status = nf90_put_var(ofid,ovarid,wet,start=[1,1,t],count=[xlen,ylen,12])
  if (status /= nf90_noerr) call handle_err(status)

end do

status = nf90_close(ifid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ofid)
if (status /= nf90_noerr) call handle_err(status)

!---------------------------------------------------------------------------------
  
contains

!----------------------------

subroutine handle_err(status)

! Internal subroutine - checks error status after each netcdf call,
! prints out text message each time an error code is returned. 

  integer, intent (in) :: status
    
  if(status /= nf90_noerr) then 
    write(0,*)trim(nf90_strerror(status))
    stop
  end if

end subroutine handle_err

!----------------------------
end program calc_wetdays 
