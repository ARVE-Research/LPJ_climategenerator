program makeclimate

! gcc -Difort -c overprint.c
! ifort -xHost -o makeclimate randomdistmod.f90 makeclimate.f90 overprint.o -I/usr/local/include -L/usr/local/lib -lnetcdff -lnetcdf
! ifort -xHost -o makeclimate randomdistmod.f90 makeclimate.f90 overprint.o -I/share1/netcdf/4.7.1-impi/include -L/share1/netcdf/4.7.1-impi/lib -lnetcdff -lnetcdf
! gfortran -O3 -o makeclimate randomdistmod.f90 makeclimate.f90 overprint.o -I/home/public/easybuild/software/netCDF-Fortran/4.5.2-gompi-2020a/include -L/home/public/easybuild/software/netCDF-Fortran/4.5.2-gompi-2020a/lib -lnetcdff -lnetcdf

! ./makeclimate 

!program to generate a climatology of arbitrary length by randomly selecting ycyc-year blocks of climate from the
!20th century reanalysis dataset

use typesizes
use netcdf
use randomdistmod, only : randomstate,ran_seed,ranur

implicit none

integer, parameter :: i2 = selected_int_kind(4)
integer, parameter :: sp = selected_real_kind(4)
integer, parameter :: dp = selected_real_kind(13)

integer(i2), parameter :: missing = -32768
real(sp), parameter :: missing_lightning = -32768.

integer, parameter :: xlen = 720
integer, parameter :: ylen = 360

integer :: ycyc
integer :: tlen

character(200) :: jobfile

character(200) :: basefile
character(200) :: anompath
character(200) :: tmpfile
character(200) :: dtrfile
character(200) :: prefile
character(200) :: cldfile
character(200) :: wndfile
character(200) :: wetfile
character(200) :: capefile

character(200) :: outfile

integer :: a,b
integer :: i,j,k
integer :: x,y
integer :: status
integer :: afid,bfid,ofid
integer :: varid
integer :: yr
integer :: day
integer :: m

integer, dimension(12) :: nd

integer(2), dimension(2) :: valid_range
real(dp),   dimension(2) :: valid_range_lightning
real(dp),   dimension(2) :: actual_range

real(sp),    allocatable, dimension(:,:,:) :: anom
real(sp),    allocatable, dimension(:,:,:) :: rbase
integer(i2), allocatable, dimension(:,:,:) :: base
real(sp),    allocatable, dimension(:,:,:) :: base_lightning
integer(i2), allocatable, dimension(:,:,:) :: vout
real(sp),    allocatable, dimension(:,:,:) :: vout_lightning

real(sp) :: scale_factor

real(dp), allocatable, dimension(:) :: time

real(dp), dimension(12) :: cape
real(dp), dimension(12) :: lanom

real(sp) :: limit

real(sp) :: add_offset

real(sp), parameter :: urs = 4.    !upper range CAPE:lightning scale factor
real(sp), parameter :: lrs = 0.95  !lower range CAPE:lightning scale factor

type(randomstate) :: rndst

integer :: startyr =   0  !year (BP) of the first year of the output climatology
integer :: numyrs  = 100  !default value in case it is not specified in the namelist
integer :: numcyc

integer, allocatable, dimension(:) :: seg
integer, allocatable, dimension(:) :: offset

logical, dimension(12) :: used

integer :: val

!index of the starting value of the 30-year climate blocks

integer, parameter, dimension(12) :: offset30  = [ 1, 361, 721, 1081, 241, 601, 961, 1321, 121, 481, 841, 1201 ]

integer :: climmonths

integer :: seed

character(80) :: status_msg

!character(10) :: cval

character(30) :: datestring
character(30) :: calendar

logical :: timeslice = .true.

namelist  / joboptions / basefile,anompath,tmpfile,dtrfile,prefile,cldfile,wndfile,wetfile,capefile,startyr,numyrs,timeslice

!----------------------------

call getarg(1,jobfile)
call getarg(2,outfile)

open(10,file=jobfile,status='old')
read(10,nml=joboptions)
close(10)

if (timeslice) then

  ycyc = 30
  tlen = 12 * ycyc

  seed = -(32768 + startyr)

  call ran_seed(seed,rndst)   !initialize the random seed with the start year

  !------
  !select periods for construction of climatology

  numcyc = numyrs / 30

  if (mod(numyrs,30) /= 0) numcyc = numcyc + 1

  write(0,'(i5,a)')numyrs,' years of climate requested'
  write(0,'(a,2i5,a,i5,a)')' will generate',numcyc,ycyc,'-year climate cycles = ',numcyc*ycyc,' years of climate'
  write(0,*)'generating sequence order'

  !generate a pseudo-random sequence of 30 year climate blocks
  !avoid repeating some blocks frequently and never using others by setting a flag vector
  !to keep choosing random numbers until one that has not already been used is selected
  !once all values have been used once, reset the flag vector and start again

  open(20,file='sequence.dat',status='unknown')

  allocate(seg(numcyc))
  allocate(offset(numcyc))

  used = .false.

  do i = 1,numcyc

    do
      val = nint(1. + 11. * ranur(rndst))  !generates a random integer between 1 and 12

      if (.not.used(val)) exit
    end do

    seg(i) = val
    used(val) = .true.

    if (all(used)) used = .false.

    write(20,*)i,seg(i)

    offset(i) = offset30(seg(i))

    write(0,*)seg(i),offset(i)

  end do

  close(20)

  datestring = 'days since 0000-01-01 00:00:00'
  calendar   = '365_day'
  startyr    = -1

else

  !transient climatology; use the 1871-2010 transient values

  ycyc = 20
  tlen = 12 * ycyc

  numcyc = 7

  allocate(seg(numcyc))
  allocate(offset(numcyc))

  do i = 1,numcyc
    seg(i) = i
    offset(i) = 1 + tlen * (i - 1)

    write(0,*)seg(i),offset(i)

  end do

  datestring = 'days since 1871-01-01 00:00:00'
  calendar   = 'gregorian'
  startyr    = 1871

  write(0,'(a,2i5,a,i5,a)')' transient climatology:',numcyc,ycyc,'-year climate cycles = ',numcyc*ycyc,' years of climate'

end if

!------
!open files

status = nf90_open(basefile,nf90_nowrite,bfid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_open(outfile,nf90_write,ofid)
if (status /= nf90_noerr) call handle_err(status)

!inquire dimensions and allocate arrays

allocate(base(xlen,ylen,12))
allocate(base_lightning(xlen,ylen,12))
allocate(rbase(xlen,ylen,12))
allocate(anom(xlen,ylen,1680))
allocate(vout(xlen,ylen,tlen))
allocate(vout_lightning(xlen,ylen,tlen))

!----------------------------------------------------------------------------
!time

!NB this section should be adjusted to accommodate leap years!

write(0,*)
write(0,*)'writing time'

climmonths = tlen * numcyc

allocate(time(climmonths))

day = 0
i = 1

yr = startyr

do y = 1,ycyc * numcyc  !years
  nd = ndaymon(yr)
  do m = 1,12           !months

    time(i) = real(day)
    i = i + 1
    day = day + nd(m)

  end do

  if (.not.timeslice) yr = yr + 1

end do

status = nf90_inq_varid(ofid,'time',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,time)
if (status /= nf90_noerr) call handle_err(status)

actual_range(1) = minval(time)
actual_range(2) = maxval(time)

status = nf90_put_att(ofid,varid,'actual_range',actual_range)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'units',trim(datestring))
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ofid,varid,'calendar',trim(calendar))
if (status /= nf90_noerr) call handle_err(status)

!----------------------------------------------------------------------------
!elevation

write(0,*)
write(0,*)'reading elevation'

status = nf90_inq_varid(bfid,'elv',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,base(:,:,1))
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(ofid,'elv',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ofid,varid,base(:,:,1))
if (status /= nf90_noerr) call handle_err(status)

!----------------------------------------------------------------------------
!temperature

write(0,*)
write(0,*)'reading temperature'

status = nf90_inq_varid(bfid,'tmp',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,base)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_open(trim(anompath)//tmpfile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'air',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,anom)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!----

valid_range = 0

status = nf90_inq_varid(ofid,'tmp',varid)
if (status /= nf90_noerr) call handle_err(status)

rbase = real(base) * scale_factor

write(0,'(a,2f8.5)')' calculating',scale_factor,add_offset

vout = missing

!work in blocks of 30 years

write(0,*) 'Creating 30-year-blocks for tmp'

do j = 1,numcyc

  k = offset(j)

  do i = 1,tlen,12

    a = k + i - 1
    b = a + 11

    where (base /= missing) vout(:,:,i:11+i) = nint((rbase + anom(:,:,a:b)) / scale_factor)

  end do

  m = 1 + tlen * (j-1)

 write(status_msg,*)'writing block',j,seg(j),k,m
 call overprint(status_msg)

  status = nf90_put_var(ofid,varid,vout,start=[1,1,m],count=[xlen,ylen,tlen])
  if (status /= nf90_noerr) call handle_err(status)

  valid_range(1) = min(valid_range(1),minval(vout,mask = vout /= missing))
  valid_range(2) = max(valid_range(2),maxval(vout,mask = vout /= missing))

end do

write(0,*)
write(0,'(a,2f6.1)')' writing valid range',real(valid_range) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

!----------------------------------------------------------------------------
!DTR temperature range

write(0,*)
write(0,*)'reading DTR'

status = nf90_inq_varid(bfid,'dtr',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,base)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_open(trim(anompath)//dtrfile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'dtr',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,anom)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!----

valid_range = 0

status = nf90_inq_varid(ofid,'dtr',varid)
if (status /= nf90_noerr) call handle_err(status)

rbase = real(base) * scale_factor

write(0,'(a,2f8.5)')' calculating',scale_factor,add_offset

vout = missing

!work in blocks of 30 years

write(0,*) 'Creating 30-year-blocks for dtr'

do j = 1,numcyc

  k = offset(j)

  do i = 1,tlen,12

    a = k + i - 1
    b = a + 11

    where (base /= missing) vout(:,:,i:11+i) = nint(max(rbase + anom(:,:,a:b),0.) / scale_factor)

  end do

  m = 1 + tlen * (j-1)

 write(status_msg,*)'writing block',j,seg(j),k,m
 call overprint(status_msg)

  status = nf90_put_var(ofid,varid,vout,start=[1,1,m],count=[xlen,ylen,tlen])
  if (status /= nf90_noerr) call handle_err(status)

  valid_range(1) = min(valid_range(1),minval(vout,mask = vout /= missing))
  valid_range(2) = max(valid_range(2),maxval(vout,mask = vout /= missing))

end do

write(0,*)
write(0,'(a,2f6.1)')' writing valid range',real(valid_range) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

!------
!precipitation

write(0,*)'reading precipitation'

status = nf90_inq_varid(bfid,'pre',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,base)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_open(trim(anompath)//prefile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'prec',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,anom)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!----

valid_range = 0

status = nf90_inq_varid(ofid,'pre',varid)
if (status /= nf90_noerr) call handle_err(status)

rbase = real(base) * scale_factor + add_offset

write(0,'(a,f8.5,f8.1)')' calculating',scale_factor,add_offset

vout = missing

!work in blocks of 30 years

write(0,*) 'Creating 30-year-blocks for pre'

do j = 1,numcyc

  k = offset(j)

  do i = 1,tlen,12

    a = k + i - 1
    b = a + 11

    where (base /= missing) vout(:,:,i:11+i) = nint((max(0.,rbase + anom(:,:,a:b)) - add_offset) / scale_factor)

  end do

  m = 1 + tlen * (j-1)

 write(status_msg,*)'writing block',j,seg(j),k,m
 call overprint(status_msg)

  status = nf90_put_var(ofid,varid,vout,start=[1,1,m],count=[xlen,ylen,tlen])
  if (status /= nf90_noerr) call handle_err(status)

  valid_range(1) = min(valid_range(1),minval(vout,mask = vout /= missing))
  valid_range(2) = max(valid_range(2),maxval(vout,mask = vout /= missing))

end do

write(0,*)
write(0,'(a,2f8.1)')' writing valid range',real(valid_range) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

!----------------------------------------------------------------------------
!wet days

write(0,*)
write(0,*)'reading wet days'

status = nf90_inq_varid(bfid,'wet',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,base)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_open(trim(anompath)//wetfile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'wet',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,anom)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!----

valid_range = 0

status = nf90_inq_varid(ofid,'wet',varid)
if (status /= nf90_noerr) call handle_err(status)

rbase = real(base) * scale_factor

write(0,'(a,2f8.5)')' calculating',scale_factor,add_offset

vout = missing

!work in blocks of years

write(0,*) 'Creating 30-year-blocks for wet'

yr = startyr

do j = 1,numcyc

  k = offset(j)

  do i = 1,tlen,12

    a = k + i - 1
    b = a + 11

    nd = ndaymon(yr)

    do y = 1,ylen
      do x = 1,xlen

        if (base(x,y,1) == missing) cycle

        vout(x,y,i:11+i) = nint(max(0.,min(real(nd),rbase(x,y,:) + anom(x,y,a:b))) / scale_factor)

      end do
    end do

    if (.not.timeslice) yr = yr + 1

  end do

  m = 1 + tlen * (j-1)

 write(status_msg,*)'writing block',j,seg(j),k,m
 call overprint(status_msg)

  status = nf90_put_var(ofid,varid,vout,start=[1,1,m],count=[xlen,ylen,tlen])
  if (status /= nf90_noerr) call handle_err(status)

  valid_range(1) = min(valid_range(1),minval(vout,mask = vout /= missing))
  valid_range(2) = max(valid_range(2),maxval(vout,mask = vout /= missing))

end do

write(0,*)
write(0,'(a,2f6.1)')' writing valid range',real(valid_range) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

!----------------------------------------------------------------------------
!windspeed

write(0,*)
write(0,*)'reading wind'

status = nf90_inq_varid(bfid,'wnd',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,base)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_open(trim(anompath)//wndfile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'wspd',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,anom)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!----

valid_range = 0

status = nf90_inq_varid(ofid,'wnd',varid)
if (status /= nf90_noerr) call handle_err(status)

rbase = real(base) * scale_factor

write(0,'(a,2f8.5)')' calculating',scale_factor,add_offset

vout = missing

!work in blocks of 30 years

write(0,*) 'Creating 30-year-blocks for wnd'

do j = 1,numcyc

  k = offset(j)

  do i = 1,tlen,12

    a = k + i - 1
    b = a + 11

    where (base /= missing) vout(:,:,i:11+i) = nint(max(0.,rbase + anom(:,:,a:b)) / scale_factor)

  end do

  m = 1 + tlen * (j-1)

 write(status_msg,*)'writing block',j,seg(j),k,m
 call overprint(status_msg)

  status = nf90_put_var(ofid,varid,vout,start=[1,1,m],count=[xlen,ylen,tlen])
  if (status /= nf90_noerr) call handle_err(status)

  valid_range(1) = min(valid_range(1),minval(vout,mask = vout /= missing))
  valid_range(2) = max(valid_range(2),maxval(vout,mask = vout /= missing))

end do

write(0,*)
write(0,'(a,2f6.1)')' writing valid range',real(valid_range) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

!----------------------------------------------------------------------------
!cloud cover fraction

write(0,*)
write(0,*)'reading cloud'

status = nf90_inq_varid(bfid,'cld',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,base)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_open(trim(anompath)//cldfile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'tcdc',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,anom)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!----

valid_range = 0

status = nf90_inq_varid(ofid,'cld',varid)
if (status /= nf90_noerr) call handle_err(status)

rbase = real(base) * scale_factor

write(0,'(a,2f8.5)')' calculating',scale_factor,add_offset

vout = missing

!work in blocks of 30 years

write(0,*) 'Creating 30-year-blocks for cld'

do j = 1,numcyc

  k = offset(j)

  do i = 1,tlen,12

    a = k + i - 1
    b = a + 11

    where (base /= missing) vout(:,:,i:11+i) = nint(min(100.,max(0.,rbase + anom(:,:,a:b))) / scale_factor)
  end do

  m = 1 + tlen * (j-1)

 write(status_msg,*)'writing block',j,seg(j),k,m
 call overprint(status_msg)

  status = nf90_put_var(ofid,varid,vout,start=[1,1,m],count=[xlen,ylen,tlen])
  if (status /= nf90_noerr) call handle_err(status)

  valid_range(1) = min(valid_range(1),minval(vout,mask = vout /= missing))
  valid_range(2) = max(valid_range(2),maxval(vout,mask = vout /= missing))

end do

write(0,*)
write(0,'(a,2f6.1)')' writing valid range',real(valid_range) * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range)
if (status /= nf90_noerr) call handle_err(status)

!----------------------------------------------------------------------------
!lightning

write(0,*)
write(0,*)'reading lightning'

status = nf90_inq_varid(bfid,'lght',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(bfid,varid,base_lightning)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (status /= nf90_noerr) call handle_err(status)

!----

status = nf90_open(trim(anompath)//capefile,nf90_nowrite,afid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_varid(afid,'cape',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_get_var(afid,varid,anom)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(afid)
if (status /= nf90_noerr) call handle_err(status)

!----

valid_range_lightning = 0.

status = nf90_inq_varid(ofid,'lght',varid)
if (status /= nf90_noerr) call handle_err(status)

rbase = base_lightning * scale_factor + add_offset

write(0,'(a,2f8.5)')' calculating',scale_factor,add_offset

vout_lightning = missing_lightning

!work in blocks of 30 years

write(0,*) 'Creating 30-year-blocks for lght'

do j = 1,numcyc

  k = offset(j)

  do y = 1,ylen
    do x = 1,xlen

      if (any(base_lightning(x,y,:) /= missing_lightning)) then

        limit = max(abs(maxval(anom(x,y,:))),abs(minval(anom(x,y,:)))) !absolute maximum magnitude of CAPE anomaly

        do i = 1,tlen,12

          a = k + i - 1
          b = a + 11

          cape = anom(x,y,a:b) / limit       !12 values, normalized to -1,+1

          !if CAPE anomaly is positive, can be up to 5x mean value
          !if CAPE anomaly is negative, minimum lightning will be 0.05x mean value
          !this range of values around the mean are supported by qualitative inspection
          !of 25 years of data from the the Alaska lightning strike data set

          where (cape == 0.)
            lanom = rbase(x,y,:)
          elsewhere (cape > 0.)
            lanom = rbase(x,y,:) * (urs * cape + 1.)
          elsewhere
            lanom = rbase(x,y,:) * (lrs * cape + 1.)
          end where

          vout_lightning(x,y,i:11+i) = (max(0.,lanom) - add_offset) / scale_factor

        end do
      else

        vout_lightning(x,y,:) = (0. - add_offset) / scale_factor

      end if
    end do
  end do

  m = 1 + tlen * (j-1)

 write(status_msg,*)'writing block',j,seg(j),k,m
 call overprint(status_msg)

  status = nf90_put_var(ofid,varid,vout_lightning,start=[1,1,m],count=[xlen,ylen,tlen])
  if (status /= nf90_noerr) call handle_err(status)

  valid_range_lightning(1) = minval(vout_lightning,mask = vout_lightning /= missing_lightning)
  valid_range_lightning(2) = maxval(vout_lightning,mask = vout_lightning /= missing_lightning)

end do

write(0,*)
write(0,'(a,2f6.1)')' writing valid range',valid_range_lightning ! * scale_factor + add_offset

status = nf90_put_att(ofid,varid,'valid_range',valid_range_lightning)
if (status /= nf90_noerr) call handle_err(status)

!----------------------------------------------------------------------------
!finish

99 continue

status = nf90_close(bfid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ofid)
if (status /= nf90_noerr) call handle_err(status)

contains

!--------------------------------------------------------------------------------

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

function ndaymon(year)

implicit none

integer, parameter, dimension(12) :: nd0 = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]  !number of days in each month
integer, parameter, dimension(12) :: nd1 = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]  !number of days in each month (leap year)

integer, intent(in) :: year  !requires year AD as input

integer, dimension(12) :: ndaymon

!---

if (leapyear(year)) then
  ndaymon = nd1
else
  ndaymon = nd0
end if

end function ndaymon

!----------------------------

logical function leapyear(year)

implicit none

integer, intent(in) :: year  !requires year AD as input

!---

if ((mod(year,4) == 0 .and. mod(year,100) /= 0) .or. mod(year,400) == 0) then
  leapyear = .true.
else
  leapyear = .false.
end if

end function leapyear

!----------------------------

end program makeclimate
