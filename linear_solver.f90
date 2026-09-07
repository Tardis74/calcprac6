! linear_solver.f90
module linear_solver
  	use precision_mod, only: dp
  	implicit none
contains

  	function solve(A, B, method) result(X)
    		real(dp), intent(in) :: A(:,:), B(:)
    		integer, intent(in) :: method
    		real(dp), allocatable :: X(:)

    		integer :: n, k, i, loc(2)
    		real(dp), allocatable :: M(:,:)      ! расширенная матрица
    		integer, allocatable :: col_perm(:)  ! перестановки столбцов
    		real(dp), allocatable :: temp_row(:), temp_col(:)
    		real(dp) ::pivot
    		real(dp), allocatable :: temp_arr(:)

    		n = size(A, 1)
    		allocate(M(n, n+1), col_perm(n))
    		allocate(temp_row(n+1), temp_col(n)) 
    		M(:, 1:n) = A
    		M(:, n+1) = B
    		col_perm = [(i, i = 1, n)]

    		select case (method)

	    		case (0)  ! Гаусс без выбора главного элемента
		      		do k = 1, n
					pivot = M(k, k)
					if (abs(pivot) < 1e-12_dp) &
			  			print *, "Warning: near-zero pivot at step", k, " value =", pivot
					! Нормализация строки k
					M(k, k:n+1) = M(k, k:n+1) / pivot
					! Исключение ниже
					!$OMP PARALLEL DO PRIVATE(i)
					do i = k+1, n
						M(i, k:n+1) = M(i, k:n+1) - M(i, k) * M(k, k:n+1)
					end do
					!$OMP END PARALLEL DO
		      		end do

				! Обратная подстановка с dot_product
		      		allocate(X(n))
		      		do i = n, 1, -1
					X(i) = M(i, n+1) - dot_product(M(i, i+1:n), X(i+1:n))
		      		end do

		    	case (1)  ! Гаусс с полным выбором
		      		do k = 1, n
					! Поиск максимального элемента в подматрице M(k:n, k:n)
					loc = maxloc(abs(M(k:n, k:n)))
					loc(1) = loc(1) + k - 1
					loc(2) = loc(2) + k - 1
					if (abs(M(loc(1), loc(2))) < 1e-12_dp) &
						print *, "Warning: near-zero pivot at step", k, " value =", M(loc(1), loc(2))
					! Перестановка строк (если нужно)
					if (loc(1) /= k) then
						temp_row = M(k, k:n+1)
						M(k, k:n+1) = M(loc(1), k:n+1)
						M(loc(1), k:n+1) = temp_row
					end if

					! Перестановка столбцов (если нужно)
					if (loc(2) /= k) then
						temp_col = M(:, k)
						M(:, k) = M(:, loc(2))
						M(:, loc(2)) = temp_col
						! Обновляем перестановку столбцов
						i = col_perm(k)
						col_perm(k) = col_perm(loc(2))
						col_perm(loc(2)) = i
					end if

					! Нормализация строки k
					M(k, k:n+1) = M(k, k:n+1) / M(k, k)

					! Исключение ниже
					!$OMP PARALLEL DO PRIVATE(i)
					do i = k+1, n
						M(i, k:n+1) = M(i, k:n+1) - M(i, k) * M(k, k:n+1)
					end do
					!$OMP END PARALLEL DO
			      	end do

				! Обратная подстановка
				allocate(X(n))
				do i = n, 1, -1
					X(i) = M(i, n+1) - dot_product(M(i, i+1:n), X(i+1:n))
				end do

				! Восстановление порядка переменных
				allocate(temp_arr(n))
			      	temp_arr = X
				do i = 1, n
					X(col_perm(i)) = temp_arr(i)
			      	end do
			      	deallocate(temp_arr)

			case (2)  ! Жордан
	      			do k = 1, n
					! Поиск максимального элемента
					loc = maxloc(abs(M(k:n, k:n)))
					loc(1) = loc(1) + k - 1
					loc(2) = loc(2) + k - 1
					if (abs(M(loc(1), loc(2))) < 1e-12_dp) &
		  				print *, "Warning: near-zero pivot at step", k, " value =", M(loc(1), loc(2))

					! Перестановка строк
					if (loc(1) /= k) then
		  				temp_row = M(k, k:n+1)
		  				M(k, k:n+1) = M(loc(1), k:n+1)
		  				M(loc(1), k:n+1) = temp_row
					end if

					! Перестановка столбцов
					if (loc(2) /= k) then
			  			temp_col = M(:, k)
			  			M(:, k) = M(:, loc(2))
			  			M(:, loc(2)) = temp_col
			  			i = col_perm(k)
			  			col_perm(k) = col_perm(loc(2))
			  			col_perm(loc(2)) = i
					end if

					! Нормализация строки k
					M(k, k:n+1) = M(k, k:n+1) / M(k, k)

					! Исключение во всех остальных строках
					!$OMP PARALLEL DO PRIVATE(i)
					do i = 1, n
		  				if (i == k) cycle
		  				M(i, k:n+1) = M(i, k:n+1) - M(i, k) * M(k, k:n+1)
					end do
					!$OMP END PARALLEL DO
	      			end do

	      			! Решение – последний столбец
	      			allocate(X(n))
	      			X = M(:, n+1)

	      			! Восстановление порядка переменных
	      			allocate(temp_arr(n))
	      			temp_arr = X
	      			do i = 1, n
					X(col_perm(i)) = temp_arr(i)
	      			end do
	      			deallocate(temp_arr)

	    		case default
	      			print *, "Error: unknown method", method
	      			allocate(X(0))
	      			return
	    	end select

    		deallocate(M, col_perm)
	end function solve

end module linear_solver
