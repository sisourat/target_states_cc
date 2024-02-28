module general
implicit none

integer, parameter :: lenmax = 200
double precision, parameter :: pi = dacos(-1d0)
double complex, parameter :: imag = dcmplx(0d0,1d0)

end module general

module mycsf_lib
implicit none

 double precision, dimension(:,:,:,:), allocatable :: mycsf_basis
 double precision, dimension(:,:), allocatable :: mycoef_det
 integer :: myndetmax, myncsfmax

end module mycsf_lib


