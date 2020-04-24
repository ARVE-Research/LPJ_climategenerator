program addanom

!program to generate climate from a baseline present-day climatology and paleoclimate anomalies from a GCM
!Jed Kaplan, January 2020

! ifort -xHost -o addanom addanom.f90 -I/usr/local/include -L/usr/local/lib -lnetcdff -lnetcdf
! ifort -xHost -o addanom addanom.f90 -I/share1/netcdf/4.7.1-impi/include -L/share1/netcdf/4.7.1-impi/lib -lnetcdff -lnetcdf

use netcdf

implicit none

integer, parameter :: i2 = selected_int_kind(4)
integer, parameter :: i4 = selected_int_kind(8)
integer, parameter :: i8 = selected_int_kind(13)
integer, parameter :: sp = selected_real_kind(4)
integer, parameter :: dp = selected_real_kind(13)

integer, parameter, dimension(12) :: ndaymon = [ 31,28,31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]  !number of days in each month

real(sp), parameter :: nsecday = 86400.  !number of seconds in 24hrs

character(200) :: basefile
character(200) :: jobfile
character(200) :: outfile

character(200) :: anompath
character(200) :: tmpfile
character(200) :: dtrfile
character(200) :: cldfile
character(200) :: wetfile
character(200) :: prefile
character(200) :: windfile
character(200) :: lghtfile

character(10) :: exp

character(100) :: note1
character(100) :: note2

integer :: status
integer :: bfid
integer :: afid
integer :: ofid
integer :: dimid
integer :: varid

integer :: xlen
integer :: ylen
integer :: tlen

integer(i2), dimension(2) :: valid_range

real(sp), allocatable, dimension(:) :: mvals

integer(i2), allocatable, dimension(:,:,:) :: varb
real(sp),    allocatable, dimension(:,:,:) :: vara
real(sp),    allocatable, dimension(:,:,:) :: tanom
integer(i2), allocatable, dimension(:,:,:) :: outv

real(sp)   :: scale_factor
real(sp)   :: add_offset
integer(i2) :: missing

integer :: x
integer :: y

logical, allocatable, dimension(:,:) :: land

real(sp), parameter :: rmissing = -9999.

real(sp) :: anom_missing

namelist / anomfiles / basefile,anompath,tmpfile,dtrfile,prefile,wetfile,windfile,cldfile,lghtfile

!-----------------------------------------------

call getarg(1,jobfile)
call getarg(2,outfile)

open(10,file=jobfile,status='old')
read(10,nml=anomfiles)
close(10)

!----------------------------
!setup

status = nf90_open(basefile,nf90_nowrite,bfid)
if (status /= nf90_noerr) call handle_err(status)

write(0,*) 'using basefile: ',trim(basefile)

status = nf90_open(outfile,nf90_write,ofid)
if (status /= nf90_noerr) call handle_err(status)

!---

status = nf90_inq_dimid(bfid,'lon',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(bfid,dimid,len=xlen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(bfid,'lat',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(bfid,dimid,len=ylen)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(bfid,'time',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(bfid,dimid,len=tlen)
if (status /= nf90_noerr) call handle_err(status)

write(0,*) 'inquired dimensions of baseline file'

allocate(mvals(tlen))
allocate(varb(xlen,ylen,tlen))
allocate(vara(xlen,ylen,tlen))
allocate(tanom(xlen,ylen,tlen))
allocate(outv(xlen,ylen,tlen))
allocate(land(xlen,ylen))

write(0,*) 'allocated variables'

!----------------------------------------------------------------------------------------------------------------------------------------------------
!1. temperature
! baseline

write(0,*) '===== temperature ====='

status = nf90_inq_varid(bfid,'tmp',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,varb)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

write(0,*) 'Got scale factor attribute'

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

write(0,*) 'Got offset attribute'

status = nf90_get_att(bfid,varid,'missing_value',missing)
if (status /= nf90_noerr) call handle_err(status)

write(0,*) 'Got missing value attribute'

!---
! anomalies

status = nf90_open(trim(anompath)//tmpfile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'tsurf',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,tanom)
if (status /= nf90_noerr) call handle_err(status)

!---

land = .false.
outv = missing

do y = 1,ylen
  do x = 1,xlen

    if (all(varb(x,y,:) /= missing) .and. all(tanom(x,y,:) /= rmissing)) then

      ! convert baseline temperature to actual temperature using scale factor and offset
      mvals = real(varb(x,y,:)) * scale_factor + add_offset

     	! add temperature anomaly to converted baseline temperature
      mvals = mvals + tanom(x,y,:)

      ! re-convert output temperature to integer using offset and scale factor
      outv(x,y,:) = nint((mvals - add_offset) / scale_factor)

      land(x,y) = .true.

    end if
  end do
end do

status = nf90_inq_varid(ofid,'tmp',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,outv)
if (status /= nf90_noerr) call handle_err(status)

valid_range(1) = minval(outv,mask = outv /= missing)
valid_range(2) = maxval(outv,mask = outv /= missing)

write(0,*)'output data range:',valid_range(1) * scale_factor + add_offset,valid_range(2) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!----------------------------------------------------------------------------------------------------------------------------------------------------
!2. precipitation

write(0,*) '===== precipitation  ====='

status = nf90_inq_varid(bfid,'pre',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,varb)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

!status = nf90_get_att(bfid,varid,'missing_value',missing)
!if (status /= nf90_noerr) call handle_err(status)

!---
! anomalies

status = nf90_open(trim(anompath)//prefile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'prec',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,vara)
if (status /= nf90_noerr) call handle_err(status)

!---

outv = missing

do y = 1,ylen
  do x = 1,xlen

    if (all(varb(x,y,:) /= missing) .and. all(tanom(x,y,:) /= rmissing)) then

      mvals = real(varb(x,y,:)) * scale_factor + add_offset

      mvals = max(mvals + vara(x,y,:),0.)  !limit to positive values

      outv(x,y,:) = nint((mvals - add_offset) / scale_factor)

    end if
  end do
end do

status = nf90_inq_varid(ofid,'pre',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,outv)
if (status /= nf90_noerr) call handle_err(status)

valid_range(1) = minval(outv,mask = outv /= missing)
valid_range(2) = maxval(outv,mask = outv /= missing)

write(0,*)'output data range:',valid_range(1) * scale_factor + add_offset,valid_range(2) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!----------------------------------------------------------------------------------------------------------------------------------------------------
!3. cloud fraction

write(0,*) '===== cloud fraction ====='

status = nf90_inq_varid(bfid,'cld',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,varb)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'missing_value',missing)
if (status /= nf90_noerr) call handle_err(status)

!---
! anomalies

status = nf90_open(trim(anompath)//cldfile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'pcldt',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,vara)
if (status /= nf90_noerr) call handle_err(status)

!---

outv = missing

do y = 1,ylen
  do x = 1,xlen

    if (all(varb(x,y,:) /= missing) .and. all(vara(x,y,:) /= rmissing)) then

      vara(x,y,:) = vara(x,y,:)

      mvals = real(varb(x,y,:)) * scale_factor + add_offset		! baseline variable is in percent

      mvals = min(max(mvals + vara(x,y,:),0.),100.)  !limit to values between 0 and 100%

      outv(x,y,:) = nint((mvals - add_offset) / scale_factor)

    end if
  end do
end do

status = nf90_inq_varid(ofid,'cld',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,outv)
if (status /= nf90_noerr) call handle_err(status)

valid_range(1) = minval(outv,mask = outv /= missing)
valid_range(2) = maxval(outv,mask = outv /= missing)

write(0,*)'output data range:',valid_range(1) * scale_factor + add_offset,valid_range(2) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!------------------
!4. wind speed

write(0,*) '===== wind speed ====='

status = nf90_inq_varid(bfid,'wnd',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,varb)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'missing_value',missing)
if (status /= nf90_noerr) call handle_err(status)

!---
! anomalies

status = nf90_open(trim(anompath)//windfile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'wsurf',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,vara)
if (status /= nf90_noerr) call handle_err(status)

!---

outv = missing

do y = 1,ylen
  do x = 1,xlen

    if (all(varb(x,y,:) /= missing) .and. all(vara(x,y,:) /= rmissing)) then

      mvals = real(varb(x,y,:)) * scale_factor + add_offset

      mvals = max(mvals + vara(x,y,:),0.)  !limit to positive values

      outv(x,y,:) = nint((mvals - add_offset) / scale_factor)

    end if
  end do
end do

status = nf90_inq_varid(ofid,'wnd',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,outv)
if (status /= nf90_noerr) call handle_err(status)

valid_range(1) = minval(outv,mask = outv /= missing)
valid_range(2) = maxval(outv,mask = outv /= missing)

write(0,*)'output data range:',valid_range(1) * scale_factor + add_offset,valid_range(2) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!------------------
!5. temperature range

write(0,*) '===== dtr ====='

status = nf90_inq_varid(bfid,'dtr',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,varb)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,varb)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'missing_value',missing)
if (status /= nf90_noerr) call handle_err(status)

!---
! anomalies

status = nf90_open(trim(anompath)//dtrfile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'dtdiurn',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,vara)
if (status /= nf90_noerr) call handle_err(status)

!---

land = .false.
outv = missing

do y = 1,ylen
  do x = 1,xlen

    if (all(varb(x,y,:) /= missing) .and. all(vara(x,y,:) /= rmissing)) then

      mvals = real(varb(x,y,:)) * scale_factor + add_offset

      mvals = max(mvals + vara(x,y,:),0.)  !limit to positive values

      outv(x,y,:) = nint(max(mvals - add_offset,0.) / scale_factor)

      land(x,y) = .true.

    end if
  end do
end do

status = nf90_inq_varid(ofid,'dtr',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,outv)
if (status /= nf90_noerr) call handle_err(status)

valid_range(1) = minval(outv,mask = outv /= missing)
valid_range(2) = maxval(outv,mask = outv /= missing)

write(0,*)'output data range:',valid_range(1) * scale_factor + add_offset,valid_range(2) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!------------------
!6. wet days

write(0,*) '===== wet days ====='

status = nf90_inq_varid(bfid,'wet',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,varb)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

!status = nf90_get_att(bfid,varid,'missing_value',missing)
!if (status /= nf90_noerr) call handle_err(status)

!---
! anomalies

status = nf90_open(trim(anompath)//wetfile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'wetdays',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,vara)
if (status /= nf90_noerr) call handle_err(status)

!---

land = .false.
outv = missing

do y = 1,ylen
  do x = 1,xlen

    if (all(varb(x,y,:) /= missing) .and. all(tanom(x,y,:) /= rmissing)) then

      mvals = real(varb(x,y,:)) * scale_factor + add_offset

      mvals = mvals + vara(x,y,:)

      mvals = max(min(mvals,real(ndaymon)),0.)   !limit to positive values not more than ndaymon

      outv(x,y,:) = nint((mvals - add_offset) / scale_factor)

      land(x,y) = .true.

    end if
  end do
end do

status = nf90_inq_varid(ofid,'wet',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,outv)
if (status /= nf90_noerr) call handle_err(status)

valid_range(1) = minval(outv,mask = outv /= missing)
valid_range(2) = maxval(outv,mask = outv /= missing)

write(0,*)'output data range:',valid_range(1) * scale_factor + add_offset,valid_range(2) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!------------------
!7. lightning

write(0,*) '===== lightning ====='


status = nf90_inq_varid(bfid,'lght',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,varb)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'missing_value',missing)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'avg_period',note1)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'source',note2)
if (status /= nf90_noerr) call handle_err(status)

!---
! anomalies

status = nf90_open(trim(anompath)//lghtfile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'L',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,vara)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(afid,varid,'_FillValue',anom_missing)
if (status /= nf90_noerr) call handle_err(status)

!---

outv = missing

do y = 1,ylen
  do x = 1,xlen

    if (all(varb(x,y,:) /= missing) .and. all(vara(x,y,:) /= rmissing)) then

      mvals = real(varb(x,y,:)) * scale_factor + add_offset

      mvals = max(mvals + vara(x,y,:),0.)  !limit to positive values

      outv(x,y,:) = nint((mvals - add_offset) / scale_factor)

    end if
  end do
end do

status = nf90_inq_varid(ofid,'lght',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,outv)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'avg_period',note1)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'source',note2)
if (status /= nf90_noerr) call handle_err(status)

valid_range(1) = minval(outv,mask = outv /= missing)
valid_range(2) = maxval(outv,mask = outv /= missing)

write(0,*)'output data range:',valid_range(1) * scale_factor + add_offset,valid_range(2) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!------------------
!clean up and close climate baseline and anomaly files

deallocate(vara)

status = nf90_close(bfid)
if (status /= nf90_noerr) call handle_err(status)

!------------------
!8. elevation

status = nf90_open(basefile,nf90_nowrite,bfid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(bfid,'elv',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,varb)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'elv',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,varb)
if (status /= nf90_noerr) call handle_err(status)

!------------

status = nf90_close(bfid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ofid)
if (status /= nf90_noerr) call handle_err(status)

write(0,*) 'Done!'

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

end program addanom
