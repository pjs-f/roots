module roots

  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none

  logical, private, parameter :: debug = .true.

  abstract interface
     real(wp) function func(x)
       import :: wp
       implicit none
       real(wp), intent(in) :: x
     end function func
  end interface

contains

  subroutine biseccion ( f, a, b, n, tol, raiz, clave )

    ! METODO DE BISECCION para encontrar una solución
    ! de f(x)=0 dada la función continua f en el intervalo
    ! [a,b] donde f(a) y f(b) tienen signos opuestos.

    ! Argumentos
    procedure(func)          :: f     ! Función que define la ecuación
    real(wp), intent(in)     :: a     ! Extremo izquierdo del intervalo inicial
    real(wp), intent(in)     :: b     ! Extremo derecho del intervalo inicial
    integer,  intent(in out) :: n     ! Límite de iteraciones/iteraciones realizadas
    real(wp), intent(in)     :: tol   ! Tolerancia para el error absoluto
    real(wp), intent(out)    :: raiz  ! Aproximación a la raiz
    integer,  intent(out)    :: clave ! Clave de éxito: 
                                      !   0 : éxito
                                      !  >0 : iteraciones excedidas
                                      !  <0 : no se puede proceder (f de 
                                      !       igual signo en a y b)
    ! Variables locales
    integer  :: i
    real(wp) :: xl, xr, error

    clave = 1
    xl    = a
    xr    = b  
    if ( sign( 1.0_wp , f(xl) ) * sign ( 1.0_wp, f(xr) ) > 0.0_wp ) then
       clave = -1
       return
    endif

    if ( debug ) write(*,'(a)') '    i          x_i'

    do i=1,n

       error = (xr-xl)*0.5_wp
       raiz  = xl + error

       if ( debug ) write(*,'(i5,2x,g0)') i, raiz

       if ( error <= tol ) then
          clave = 0
          n     = i
          exit
       endif

       if ( sign( 1.0_wp, f(xl) ) * sign( 1.0_wp, f(raiz) ) > 0.0_wp ) then
          xl = raiz
       else
          xr = raiz
       endif

    enddo

  end subroutine biseccion

  
  subroutine newton ( f, df, x0, n, tol, raiz, clave )

    ! Metodo DE NEWTON-RAPHSON para encontrar una 
    ! solución de f(x)=0 dada la función derivable
    ! f y una aproximación inicial x0.

    ! Argumentos
    procedure(func)          :: f     ! Función que define la ecuación
    procedure(func)          :: df    ! Derivada de la función que define la ecuación
    real(wp), intent(in)     :: x0    ! Aproximación inicial a la raíz
    integer,  intent(in out) :: n     ! Límite de iteraciones/iteraciones realizadas
    real(wp), intent(in)     :: tol   ! Tolerancia para el error relativo
    real(wp), intent(out)    :: raiz  ! Aproximación a la raiz
    integer,  intent(out)    :: clave ! Clave de éxito: 
                                      !   0 : éxito
                                      !  >0 : iteraciones excedidas
    ! Variables locales
    integer  ::  i
    real(wp) :: xx0

    clave = 1
    xx0   = x0

    if ( debug ) write(*,'(a)') '    i          x_i'

    do i=1,n
       raiz  = xx0 - f(xx0)/df(xx0)

       if ( debug ) write(*,'(i5,2x, g0)') i, raiz

       if ( abs( raiz-xx0 ) <= tol * abs( raiz ) ) then
          clave = 0
          n     = i
          exit
       endif

       xx0 = raiz

    end do

  end subroutine newton

  
  subroutine secante ( f, x0, x1, n, tol, raiz, clave )

    ! ALGORITMO DE LA SECANTE para encontrar una solución
    ! de f(x)=0, siendo f una función continua, dada las
    ! aproximaciones iniciales x0 y x1.

    ! Argumentos
    procedure(func)          :: f     ! Función que define la ecuación
    real(wp), intent(in)     :: x0,x1 ! Aproximaciones iniciales a la raíz
    integer,  intent(in out) :: n     ! Límite de iteraciones/iteraciones realizadas
    real(wp), intent(in)     :: tol   ! Tolerancia para el error relativo
    real(wp), intent(out)    :: raiz  ! Aproximación a la raiz
    integer,  intent(out)    :: clave ! Clave de éxito: 
                                      !   0 : éxito
                                      !  >0 : iteraciones excedidas

    ! Variables locales
    integer :: i
    real(wp):: xx0, xx1, fx0, fx1

    clave = 1
    xx0   = x0
    xx1   = x1  
    fx0   = f(x0)
    fx1   = f(x1)

    if ( debug ) write(*,'(a)') '    i          x_i'

    do i= 1,n

       raiz  = xx1 - fx1 * ( ( xx1-xx0 ) / ( fx1-fx0 ) )

       if ( debug ) write(*,'(i5,2x, g0)') i, raiz

       if ( abs( raiz-xx1 ) <= tol * abs( raiz ) ) then
          clave = 0
          n     = i
          exit
       endif

       xx0 = xx1
       fx0 = fx1
       xx1 = raiz
       fx1 = f(raiz)

    end do

  end subroutine secante

  
  subroutine birge_vieta ( a, x0, n, tol, raiz, clave )
  
    ! METODO DE BIRGE-VIETA para resolver ECUACIONES ALGEBRAICAS:
    ! P(x) = 0 donde P es un polinomio de grado m de coeficientes reales.
    ! El método se basa en el método de Newton-Raphson implementando el
    ! esquema de Horner para la evaluación del polinomio y su derivada.

    ! Argumentos
    real(wp), intent(in)     :: a(0:)  ! Vector de m+1 elementos conteniendo
                                       ! los coeficientes del polinomio
    real(wp), intent(in)     :: x0     ! Aproximación inicial a la raíz
    real(wp), intent(in)     :: tol    ! Tolerancia para el error relativo
    integer,  intent(in out) :: n      ! Límite de iteraciones/iteraciones realizadas
    real(wp), intent(out)    :: raiz   ! Aproximación a la raiz
    integer,  intent(out)    :: clave  ! Clave de éxito: 
                                       !   0 : éxito
                                       !  >0 : iteraciones excedidas
    ! Variables locales
    integer :: m, i, j
    real(wp):: xx0,b,c

    clave = 1
    xx0   = x0
    m     = size(a)-1

    if ( debug ) write(*,'(a)') '    i          x_i'

    do i=1,n

       ! Esquema de Horner
       b = a(m)
       c = a(m)    
       do j=m-1,1,-1
          b = b*xx0+a(j)
          c = c*xx0+b
       enddo
       b = b*xx0+a(0)

       ! Método de Newton
       raiz  = xx0 - b/c

       if ( debug ) write(*,'(i5,2x,g0)') i, raiz

       if ( abs( raiz-xx0 ) <= tol * abs( raiz ) ) then
          clave = 0
          n     = i
          exit
       end if
       xx0 = raiz
    end do

  end subroutine birge_vieta

  
  subroutine punto_fijo ( f, x0, n, tol, raiz, clave )

    ! ALGORITMO DE PUNTO FIJO o DE APROXIMACIONES SUCESIVAS
    ! para encontrar una solución de x=f(x) dada una
    ! aproximación inicial x0.

    ! Argumentos
    procedure(func)          :: f     ! Función que define la ecuación
    real(wp), intent(in)     :: x0    ! Aproximación inicial a la raíz
    integer,  intent(in out) :: n     ! Límite de iteraciones/iteraciones realizadas
    real(wp), intent(in)     :: tol   ! Tolerancia para el error relativo
    real(wp), intent(out)    :: raiz  ! Aproximación a la raiz
    integer,  intent(out)    :: clave ! Clave de éxito:
                                      !   0 : éxito
                                      !  >0 : iteraciones excedidas

    ! Variables locales
    integer  :: i
    real(wp) :: xx0

    clave = 1
    xx0   = x0

    if ( debug ) write(*,'(a)') '    i          x_i'

    do i=1,n
       raiz  = f(xx0)

       if ( debug ) write(*,'(i5,2x,g0)') i, raiz

       if ( abs( raiz-xx0 ) <= tol * abs( raiz ) ) then
          clave = 0
          n     = i
          exit
       endif
       xx0 = raiz
    end do

  end subroutine punto_fijo

  
  subroutine fzero( f, b, c, re, ae, iflag )

    ! Search for a zero of a function F(X) in a given interval
    ! (B,C).  It is designed primarily for problems where F(B)
    ! and F(C) have opposite signs.
    !
    ! Original author:
    !
    !  Shampine, L. F., (SNLA)
    !  Watts, H. A., (SNLA)
    !
    ! Modified by
    !
    !  Pablo Santamaría
    !
    ! Description
    !
    !     FZERO searches for a root of the nonlinear equation F(X) = 0
    !     (when F(X) is a real function of a single real variable X),
    !     between the given values B and C until the width
    !     of the interval (B,C) has collapsed to within a tolerance
    !     specified by the stopping criterion,
    !
    !                    ABS(B-C) <= 2*MAX(RW*ABS(B),AE).
    !
    !     The method used is an efficient combination of bisection and the
    !     secant rule and is due to T. J. Dekker.
    !
    ! Arguments
    !
    !   f     :EXT   - Name of the function subprogram defining F(X).
    !                  This subprogram should have the form
    !                   REAL(WP) FUNCTION F(X)
    !                      USE iso_fortan_env, ONLY: WP => REAL64
    !                      IMPLICIT NONE
    !                      REAL(WP), INTENT(IN) :: X
    !                      F = ...
    !                   END FUNCTION F
    !                  An explicit interface should be provided
    !                  in the calling program.
    !
    !   b     :INOUT - One end of the interval (B,C).  The
    !                  value returned for B usually is the better
    !                  approximation to a zero of F.
    !
    !   c     :INOUT - The other end of the interval (B,C)
    !
    !
    !   re    :IN    - Relative error used for RW in the stopping criterion.
    !                  If the requested RE is less than machine precision,
    !                  then RW is set to approximately machine precision.
    !
    !   ae    :IN    - Absolute error used in the stopping criterion.  If
    !                  the given interval (B,C) contains the origin, then a
    !                  nonzero value should be chosen for AE.
    !
    !   iflag :OUT   - A status code.  User must check IFLAG after each
    !                  call.  Control returns to the user from FZERO in all
    !                  cases.
    !
    !                0  B is within the requested tolerance of a zero.
    !                   The interval (B,C) collapsed to the requested
    !                   tolerance, the function changes sign in (B,C), and
    !                   F(X) decreased in magnitude as (B,C) collapsed.
    !
    !                1  F(B) = 0.  However, the interval (B,C) may not have
    !                   collapsed to the requested tolerance.
    !
    !               -1  B may be near a singular point of F(X).
    !                   The interval (B,C) collapsed to the requested tol-
    !                   erance and the function changes sign in (B,C), but
    !                   F(X) increased in magnitude as (B,C) collapsed, i.e.
    !                     ABS(F(B out)) > MAX(ABS(F(B in)),ABS(F(C in)))
    !
    !               -2  No change in sign of F(X) was found although the
    !                   interval (B,C) collapsed to the requested tolerance.
    !                   The user must examine this case and decide whether
    !                   B is near a local minimum of F(X), or B is near a
    !                   zero of even multiplicity, or neither of these.
    !
    !               -3  Too many (> 500) function evaluations used.
    !
    ! References
    !
    !      L. F. Shampine and H. A. Watts, FZERO, a root-solving
    !             code, Report SC-TM-70-631, Sandia Laboratories,
    !             September 1970.
    !      T. J. Dekker, Finding a zero by means of successive
    !             linear interpolation, Constructive Aspects of the
    !             Fundamental Theorem of Algebra, edited by B. Dejon
    !             and P. Henrici, Wiley-Interscience, 1969.
    !
    ! Notes
    !
    ! This implementation is a Fortran 90 refactoring from the original
    ! dfzero.f subroutine:
    !
    ! * http://www.netlib.org/slatec/src/dfzero.f
    !
    ! with some changes from:
    !
    ! * ftp://ftp.wiley.com/public/college/math/sapcodes/f90code/saplib.f90
    !
    
    ! Arguments
    procedure(func)          :: f
    real(wp), intent(in out) :: b
    real(wp), intent(in out) :: c
    real(wp), intent(in)     :: re
    real(wp), intent(in)     :: ae
    integer,  intent(out)    :: iflag

    ! Local variables
    real(wp) :: a, acbs, acmb, cmb
    real(wp) :: fa, fb, fc, fx
    real(wp) :: p, q
    real(wp) :: aw, rw, tol
    integer  :: ic, kount, max_feval
    logical  :: force_bisect
    character(7) :: method


    ! Set error tolerances
    rw = max( re, epsilon(1.0_wp) )
    aw = max( ae, 0.0_wp )

    ! Initialize
    max_feval = 500
    ic = 0

    fb = f( b )
    fc = f( c )
    kount = 2

    a = c
    fa = fc
    acbs = abs( b-c )
    fx = max( abs(fb), abs(fc) )

    if ( debug ) then
       write(*,'(a,1x,a,11x,a,19x,a,17x,a,15x,a)') 'F-count', &
            & 'method', 'b', 'c', 'f(b)', 'f(c)'
       method = 'initial'
    end if

    ! Iteration loop
    do

       ! Perform interchange such abs(f(b)) < abs(f(c))
       if ( abs( fc ) < abs( fb ) ) then
          a  = b
          fa = fb
          b  = c
          fb = fc
          c  = a
          fc = fa
       end if

       if ( debug ) write(*,'(i5,*(2x,g0))') kount, method, b, c, fb, fc

       cmb  = 0.5_wp * (c-b)
       acmb = abs( cmb )
       tol  = max( rw*abs(b), aw ) ! or tol = rw*abs(b) + aw

       ! Test stopping criterion and function count
       ! Set iflag appropriately
       if ( acmb <= tol ) then
          if ( sign( 1.0_wp, fb ) * sign( 1.0_wp, fc ) > 0.0_wp ) then
             iflag = -2
          else if ( abs( fb ) > fx ) then
             iflag = -1
          else
             iflag = 0
          end if
          exit
       endif

       if ( abs( fb ) <= 0.0_wp ) then
          iflag = 1
          exit
       end if

       if ( kount >= max_feval ) then
          iflag = -3
          exit
       end if

       ! Calculate new iterate implicitly as b+p/q, where we arrange
       ! p >= 0.  The implicit form is used to prevent overflow.
       p = (b-a) * fb
       q = fa - fb

       if ( p < 0.0_wp ) then
          p = -p
          q = -q
       end if

       ! Update a and check for satisfactory reduction in the size of the
       ! bracketing interval every four increments.
       ! If not, force bisection.
       a  = b
       fa = fb
       ic = ic + 1
       force_bisect = .false.

       if ( ic >= 4 ) then

          if ( 8.0_wp * acmb >= acbs ) then
             force_bisect = .true.
          else
             ic   = 0
             acbs = acmb
          end if

       end if

       ! Get new iterate b
       if ( force_bisect ) then
          b = b + cmb         ! Use bisection (c+b)/2
          if (debug) method = 'bisect*'
       else if ( p <=  abs(q) * tol ) then ! If to small, increment by tolerance
          b = b + sign( tol, cmb )
          if ( debug ) method = 'minimal'
       else if ( p < cmb * q ) then ! Root ought to be between b and (c+b)/2
          b = b + p/q  ! Use secant rule
          if ( debug ) method = 'secant'
       else
          b = b + cmb  ! Use bisection (c+b)/2
          if ( debug ) method = 'bisect'
       end if

       ! Have completed computation for new iterate b.
       fb = f( b )
       kount = kount + 1

       ! Decide whether next step is interpolation or extrapolation.
       if ( sign( 1.0_wp, fb ) * sign( 1.0_wp, fc ) > 0.0_wp ) then
          c  = a
          fc = fa
       end if

    end do

  end subroutine fzero

end module roots
