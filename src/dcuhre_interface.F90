module dcuhre_interface

use multidim_integrate, only: dcuhre
use, intrinsic :: iso_c_binding

implicit none
private

real(8),     parameter :: Pi = dacos(-1.d0)
complex(8) , parameter :: Xi = dcmplx(0.d0,1.d0)
integer,     parameter :: max_num_threads = 8192 ! Max number of threads concurrently using this module

!real(c_double), bind(C,name='cfor_check_kx'), public :: cfor_check_kx(4)
!real(c_double), bind(C,name='cfor_check_ky'), public :: cfor_check_ky(4)
!integer(c_int), bind(C,name='cfor_check_counter'), public :: cfor_check_counter(4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Define an interface for a callback C function
abstract interface
    subroutine callback (k,ndim,res,fSize,par) bind(C)
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in), value :: fSize
      integer(c_int), intent(in), value :: ndim
      real(c_double), dimension(ndim) , intent(inout) :: k
      real(c_double), dimension(fSize), intent(inout) :: res
      type(c_ptr)   , intent(in), value :: par
    end subroutine callback
end interface

!!!!!!!!!!!!!!!!!!!! Interface for abort function
interface
  subroutine cfor_stop_all() bind(C, name='cfor_stop_all')
  end subroutine
end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IMPORTANT
!  module variables are NOT thread safe!
!
! Here we use a workaround, where we create an array of 
! datas, where each thread will save its data on its 
! reserved position datas[thread_id].
!
! Array of function pointers
type callback_f_ptr
    procedure(callback), pointer, nopass :: integrand
    type(c_ptr) :: param
end type
type(callback_f_ptr) :: thread_data(max_num_threads)

contains 

subroutine dcuhre_interface_2D(fSize,x1,x2,y1,y2,acc_rel,acc_abs,Res,cfunc,param,istat) bind(C,name='dcuhre_interface_2D')
    use, intrinsic :: iso_c_binding
    use omp_lib, only: omp_get_num_threads, omp_get_thread_num
    implicit none

    !Input-Variables
    type(c_funptr), intent(in), value :: cfunc !cfunc pointer
    type(c_ptr)   , intent(in), value :: param
    real(c_double), intent(in), value :: x1, x2, y1, y2, acc_rel, acc_abs
    integer(c_int), intent(in), value :: fSize
    !Output-Variables
    real(c_double), intent(inout), dimension(fSize) :: Res
    integer(c_int), intent(out) :: istat
    !Subroutine-Internal-Variables
    integer :: thread_id, total_threads
    integer :: key
    integer :: nf,ifail,neval,maxsub
    integer :: maxpts,num,ndim,numfun, nsub
    real(8) :: xmin(2), xmax(2)
    real(8), allocatable :: finest(:), absest(:)
    real(8), allocatable :: Re_tmp(:)
    !integer, pointer :: ptr_tmp

    thread_id = 1+omp_get_thread_num()
    total_threads = omp_get_num_threads()

    if( total_threads > max_num_threads ) then
        write(0,*) 'Total threads given greater than max allowed threads: ' , total_threads, max_num_threads
        stop 
    endif

    ! Translate function pointer
    call c_f_procpointer (cfunc, thread_data(thread_id)%integrand)
    thread_data(thread_id)%param = param

    !cfor_check_counter = 0

    allocate ( finest(fSize), absest(fSize), Re_tmp(fSize) )

    ! Adaptive integration method using DCUHRE library
    xmin(1) =  x1
    xmin(2) =  y1
    xmax(1) =  x2
    xmax(2) =  y2
    
    ! For these parameters see the documentation of dcuhre.f
    maxpts     =       100000000! Max number of function evaluations
    nf         =        fSize
    neval      =       50000
    num        =         65
    maxsub     = (maxpts-num)/(2*num) + 1 
    ndim       =         2
    numfun     =        fSize
    nsub       =         0
    key        =         0
    
    !write(6,'(a,i2,z16)') 'Set ptr: ', thread_id-1, loc(ptr_tmp)
    !call flush(6)
    ! compute real part
    call dcuhre(ndim,nf,xmin,xmax,0,maxpts,fortran_integrand,  &
          acc_rel,acc_abs, key, 0, finest, absest, neval,      &
          ifail, nsub)
    
    Re_tmp=finest ! Save result
    if(ifail.eq.1) then
      write(0,*) 'Dcuhre interface: cannot reach accuracy error ext: ', absest(:min(size(absest),4))
    else if(ifail.ne.0) then
      write(0,*) 'Dcuhre interface: cannot compute integration, error: ',ifail
      call flush(0)
      stop
      !stop
    endif
    istat = ifail
   
    ! Save result
    Res = Re_tmp

    !write(6,*) 'Fortran counter: ', cfor_check_counter(:)
    !deallocate ( w )
    deallocate ( finest )
    deallocate ( absest ) 
    deallocate ( Re_tmp )
end subroutine

! Integrand function
subroutine fortran_integrand(ndim,z,nfun,res)
    use OMP_LIB, only: omp_get_thread_num, omp_get_num_threads
    implicit none

    integer, parameter :: dp = SELECTED_REAL_KIND(15,60)
    integer, intent(in)  :: ndim, nfun
    real(dp), intent(in) :: z(:)
    real(dp), intent(out) :: res(:)
    real(c_double) :: res2(nfun), res_scalar, k(ndim)
    type(c_ptr) :: param
    integer  :: nfun2, tid

    tid = 1+omp_get_thread_num()

    !cfor_check_kx(tid) = z(1)
    !cfor_check_ky(tid) = z(2)
    !cfor_check_counter(tid) = cfor_check_counter(tid) + 1
 
    k   = z
    res = 0.0_dp   
    res2= 0.0_c_double
    nfun2= nfun
    param = thread_data(tid)%param

    call thread_data(tid)%integrand( k, ndim, res2, nfun2 , param)

    res = res2

end subroutine


end module
