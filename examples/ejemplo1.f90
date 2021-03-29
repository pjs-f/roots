program ejemplo1
  use iso_fortran_env, only: wp => real64
  use roots, only: biseccion
  implicit none (type, external)

  real(wp) :: a, b, tol, raiz
  integer :: n, clave

  ! Datos de iniciales
  a = 0.6_wp
  b = 0.8_wp
  tol = 0.5e-8_wp
  n = 100

  ! Determinar la raíz
  call biseccion(f,a,b,n,tol,raiz,clave)
  if (clave == 0) then
     write(*,*) 'Raíz =', raiz
     write(*,*) 'Iteraciones realizadas =', n
  else
     write(*,*) 'Error =', clave
  endif

contains

  real(wp) function f(x)
    ! Función que define la ecuación
    real(wp), intent(in) :: x
    f = cos(x) - x
  end function f

end program ejemplo1
