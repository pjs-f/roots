program ejemplo2
  use iso_fortran_env, only: wp => real64
  use roots, only: fzero
  implicit none (type, external)
  real(wp) :: a,b
  real(wp) :: rerror, aerror
  integer :: code

  a = 3.0_wp
  b = 4.0_wp
  rerror = 5.0e-8_wp
  aerror = 0.0_wp

  call fzero ( f, a, b, rerror, aerror, code )
  write(*,'(a,i0)') 'Code = ', code
  write(*,'(a,g0)') 'Root = ', a

contains

  real(wp) function f(x)
    real(wp), intent(in) :: x
    f = 1.0_wp/( x - 3.0_wp ) - 6.0_wp
  end function f

end program ejemplo2
