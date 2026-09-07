module newton_module
	use precision_mod, only: dp
	use linear_solver, only: solve !решатель из 3го задания
	implicit none
	private
	
	public :: newton, vecnorm2, vector_function_interface
	
	abstract interface
		function vector_function_interface(x) result(y)
			use precision_mod, only: dp
			real(dp), dimension(:), intent(in) :: x
			real(dp), dimension(size(x)) :: y
		end function vector_function_interface
	end interface

contains

	!Евклидова норма вектора
	function vecnorm2(vec) result(nrm)
		use precision_mod, only: dp
		real(dp), dimension(:), intent(in) :: vec
		real(dp) :: nrm
		nrm = sqrt(sum(vec**2))
	end function vecnorm2
	
	!Вычисление матрицы Якоби
	function jacobian(f, x) result(J)
		use precision_mod, only: dp
		procedure(vector_function_interface) :: f
		real(dp), dimension(:), intent(in) :: x
		real(dp), dimension(size(x), size(x)) :: J
		integer :: n, k, i
		real(dp) :: h, eps
		real(dp), dimension(:), allocatable :: fx, x_plus, f_plus
		
		n = size(x)
		eps = epsilon(1.0_dp)
		allocate(fx(n))
		fx = f(x)
		
		!$OMP PARALLEL DO PRIVATE(k, h, x_plus, f_plus, i) SHARED(J, fx, x)
		do k = 1, n
			h = sqrt(eps) * (1.0_dp + abs(x(k)))
			
			allocate(x_plus(n), f_plus(n))
			x_plus = x
			x_plus(k) = x_plus(k) + h
			f_plus = f(x_plus)
			
			do i = 1, n
				J(i, k) = (f_plus(i) - fx(i)) / h
			end do
			
			deallocate(x_plus, f_plus)
		end do
		!$OMP END PARALLEL DO
	end function jacobian

	subroutine newton(f, x0, max_iter, tol, x, info, lin_method)
		procedure(vector_function_interface) :: f
		real(dp), dimension(:), intent(in) :: x0
		integer, intent(in) :: max_iter
		real(dp), intent(in) :: tol
		real(dp), dimension(:), intent(out) :: x
		integer, intent(out) :: info
		integer, intent(in) :: lin_method
		integer :: n, iter
		real(dp), dimension(:), allocatable :: x_old, f_val, dx
		real(dp), dimension(:,:), allocatable :: J
		
		n = size(x)
		allocate(x_old(n), f_val(n), dx(n), J(n,n))
		x = x0
		info = 0 !0 - успех, 1 - превышено число итераций
		
		do iter = 1, max_iter
			x_old = x
			f_val = f(x)
			J = jacobian(f, x)
			
			!Решение системы J*dx = -F(x)
			dx = solve(J, -f_val, lin_method)
			x = x_old + dx
			
			if (vecnorm2(dx) < tol) then
				info = 0
				return
			end if
		end do
		
		!цикл завершился - лимит итераций
		info = 1
	end subroutine newton
end module newton_module
