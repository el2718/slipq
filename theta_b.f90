module theta_b0
integer(1), allocatable:: rboundary(:,:,:)
real, allocatable:: p_margin(:,:), rfs(:,:,:,:), rfe(:,:,:,:), theta_b(:,:,:)
integer:: n_lon, n_lat, n_r, i_margin_end
contains

subroutine cal_theta_b(i,j,k)
implicit none
integer:: i, j, k, t
real:: vp(0:2), vp_car(0:2), cos_arc, max_cos_arc
!----------------------------------------------------------------------------
if (rboundary(i,j,k) .eq. 21 .or. (rboundary(i,j,k) .eq. 12)) then
    if (rboundary(i,j,k) .eq. 12) vp=rfs(:,i,j,k)
    if (rboundary(i,j,k) .eq. 21) vp=rfe(:,i,j,k)
    vp_car=[cos(vp(1))*[cos(vp(0)), sin(vp(0))], sin(vp(1))]

    max_cos_arc=-1.
    do t=0, i_margin_end
        cos_arc= dot_product(vp_car, p_margin(:, t))
        if (cos_arc .gt. max_cos_arc) max_cos_arc=cos_arc
    enddo

    if (.not. (max_cos_arc .lt.  1.)) then
	    theta_b(i,j,k)= 0.
    else if   (max_cos_arc .le. -1.)  then
        theta_b(i,j,k)= 3.141592653589793
    else 
        theta_b(i,j,k)= acos(max_cos_arc)
    endif
else
    theta_b(i,j,k)=0.
endif

end subroutine cal_theta_b

end module theta_b0

program main
use theta_b0
implicit none
integer:: i, j, k, nthreads, OMP_GET_NUM_PROCS
logical:: margin_open
real:: vp(0:2)
character(len=1) :: str_aux
!----------------------------------------------------------------------------
open(unit=1, file='dimension.bin', access='stream', status='old')
read(1) nthreads, n_lon, n_lat, n_r
close(1, status='delete')

allocate(p_margin(0:2, 0:n_lon*n_lat-1))
allocate(rboundary(0:n_lon-1, 0:n_lat-1, 0:n_r-1))
allocate(rfs(0:2, 0:n_lon-1, 0:n_lat-1, 0:n_r-1))
allocate(rfe(0:2, 0:n_lon-1, 0:n_lat-1, 0:n_r-1))
allocate(theta_b(0:n_lon-1, 0:n_lat-1, 0:n_r-1))

open(unit=1, file='rboundary.bin', access='stream', status='old')
read(1) rboundary
close(1, status='delete')
open(unit=1, file='rFs.bin', access='stream', status='old')
read(1) rFs
close(1, status='delete')
open(unit=1, file='rFe.bin', access='stream', status='old')
read(1) rFe
close(1, status='delete')
!----------------------------------------------------------------------------
i_margin_end=-1
do j=0, n_lat-1
do i=0, n_lon-1
    margin_open=.false.
    if (rboundary(i,j,0) .eq. 21 .or. (rboundary(i,j,0) .eq. 12)) then
        if (i+1 .lt. n_lon) margin_open =                  (rboundary(i+1,j,0) .eq. 11)
        if (i-1 .ge. 0)     margin_open = margin_open .or. (rboundary(i-1,j,0) .eq. 11)
        if (j+1 .lt. n_lat) margin_open = margin_open .or. (rboundary(i,j+1,0) .eq. 11)
        if (j-1 .ge. 0)     margin_open = margin_open .or. (rboundary(i,j-1,0) .eq. 11)
        if (margin_open) then
            if (rboundary(i,j,0) .eq. 12) vp=rfs(:,i,j,0)
            if (rboundary(i,j,0) .eq. 21) vp=rfe(:,i,j,0)
            i_margin_end=i_margin_end+1
            p_margin(:,i_margin_end)=[cos(vp(1))*[cos(vp(0)), sin(vp(0))], sin(vp(1))]
        endif
    endif
enddo
enddo
!----------------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(i,j,k), schedule(DYNAMIC)
do i=0, n_lon-1
do j=0, n_lat-1
do k=0, n_r-1
    call  cal_theta_b(i,j,k)
enddo
enddo
enddo
!$OMP END PARALLEL DO
!----------------------------------------------------------------------------
open(unit=1, file='theta_b.bin', access='stream', status='replace')
write(1) theta_b
close(1)
deallocate(p_margin, rboundary, rfs, rfe, theta_b)

! If the pop-up window for fastqsl.exe cannot be closed automatically on some Windows systems, please uncomment this line
! call system('taskkill /im theta_b.exe /f')

end program main