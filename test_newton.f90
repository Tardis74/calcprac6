program test_newton
    use precision_mod, only: dp
    use newton_module, only: newton, vecnorm2
    implicit none

    integer, parameter :: n = 2
    real(dp), dimension(n) :: x0, x
    real(dp)               :: tol
    integer                :: max_iter, info, unit, lin_method
    real(dp)               :: norm_f

    ! Начальное приближение (немного отклоняемся от решения)
    x0 = [2.1_dp, 1.1_dp]
    tol = 1.0e-12_dp
    max_iter = 100
    lin_method = 1   ! Гаусс с полным выбором

    call newton(my_f, x0, max_iter, tol, x, info, lin_method)

    write(*, '(A, I0)') 'Статус завершения: ', info
    norm_f = norm2(my_f(x))
    write(*, '(A, ES15.7)') 'Норма F(X): ', norm_f
    write(*, '(A)') 'Найденное решение X:'
    write(*, '(2ES15.7)') x
    write(*, '(A)') 'Точное решение (2.0, 1.0):'
    write(*, '(2F15.7)') 2.0_dp, 1.0_dp

    open(newunit=unit, file='result.dat', status='replace', action='write')
    write(unit, '(ES25.16)') x
    close(unit)

contains

    function my_f(x) result(y)
        use precision_mod, only: dp
        implicit none
        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(size(x))       :: y
        real(dp), parameter :: sin2 = sin(2.0_dp)
        real(dp), parameter :: cos1 = cos(1.0_dp)
        real(dp), parameter :: log3 = log(3.0_dp)

        y(1) = sin(x(1)) - sin2 + x(2) - 1.0_dp
        y(2) = cos(x(2)) - cos1 + log(x(1) + 1.0_dp) - log3
    end function my_f

end program test_newton
