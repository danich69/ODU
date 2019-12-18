module odu
	use diffur
    use non_linear_system
    implicit none
    	
    contains
    
        subroutine Runge_Kutta(x)
        
            real(8), dimension(:, :) :: x
            real(8), dimension(0:int(interval/h) + 1) :: t
            real(8), dimension(dims) :: k1, k2, k3, k4
            
            x(1, :) = x_0
            t = 0
            
            do i = 1, size(x)/dims - 1
                
                k1 = h * f(t(i), x(i, :))
                k2 = h * f(t(i) + h/2.0, x(i, :) + k1/2.0)
                k3 = h * f(t(i) + h/2.0, x(i, :) + k2/2.0)
                k4 = h * f(t(i) + h, x(i, :) + k3)

                x(i+1, :) = x(i, :) + (k1 + 2 * k2 + 2 * k3 + k4)/6.0
                t(i+1) = t(i) + h
            
            enddo
            
        end subroutine Runge_Kutta
        
        subroutine Adams_extrapolation(x)
        
            real(8), dimension(:, :) :: x
            real(8), dimension(0:int(interval/h) + 1) :: t
            real(8), dimension(0 : extrapolate_order-1) :: A
            real(8), dimension(dims) :: s
            integer :: i, j
            
            t = 0.0
            
            call Runge_Kutta(x(1 : extrapolate_order, :))
            
            do i = 0, extrapolate_order
                t(i+1) = t(i) + h                
            enddo
            
            A = A_coefficents(extrapolate_order)
            
            do i = extrapolate_order, size(x)/dims - 1
            
                s = 0
            
                do j = 0, extrapolate_order-1
                
                    s = A(j)*f(t(i-j-1), x(i - j, :)) + s
                
                enddo
                
                x(i+1, :) = h * s + x(i, :)
                t(i+1) = t(i) + h
            
            enddo
            
        
        end subroutine Adams_extrapolation
        
        function A_coefficents(n)
        
            real(8), dimension(0 : n-1) :: A_coefficents
            real(8), dimension(0 : n-1) :: powers, new
            real(8) :: a
            integer :: n, j, i
            
            do j = 0, n-1
                
                powers = 0.0
                powers(0) = 1.0
                
                do i = 0, n-1
                
                    if (i /= j) then
                    
                        new = powers * i
                        powers = eoshift(powers, -1)
                        powers = powers + new
                        
                    endif
                    
                enddo
                
                do i = 0, n-1
                    powers(i) = powers(i)/real(i+1)
                enddo
                
                A_coefficents(j) = (-1.0_8)**j/(fact(j)*fact(n - 1 - j) )* sum(powers)
                
            enddo
        
        end function A_coefficents
        
        subroutine Adams_interpolation(x)
        
            real(8), dimension(:, :) :: x
            real(8), dimension(0:int(interval/h) + 1) :: t
            real(8), dimension(-1 : interpolate_order-2) :: B
            real(8), dimension(dims) :: s
            integer :: i, j
            
            x(1, :) = x_0
            t = 0
            
            call Runge_Kutta(x(1 : interpolate_order-1, :))
            
            do i = 1, interpolate_order-1
                t(i+1) = t(i) + h                
            enddo
            
            B = B_coefficents(interpolate_order)
            
            do i = interpolate_order - 1, int(interval/h)
            
                s = 0
            
                do j = 0, interpolate_order-2
                
                    s = B(j)*f(t(i-j), x(i - j, :)) + s
                
                enddo
                
                t(i+1) = t(i) + h
                
                call nonlinear_system_solver(x(i,:), funct, x(i+1,:))
                
            
            enddo
            
            contains
                    
                function funct(y)
        
                    real(8), dimension(:) :: y
                    real(8), dimension(size(y)) :: funct
            
                    funct = x(i, :) + h * s + h * B(-1) * f(t(i+1), y) - y
            
                end function funct
        
        end subroutine Adams_interpolation
        
        function B_coefficents(n)
        
            real(8), dimension(-1 : n-2) :: B_coefficents
            real(8), dimension(-1 : n-2) :: powers, shifter
            real(8) :: a
            integer :: n, j, i
            
            do j = -1, n-2
                
                powers = 0
                powers(-1) = 1
                
                do i = -1, n-2
                
                    if (i /= j) then
                    
                        shifter = powers * i
                        powers = eoshift(powers, -1)
                        powers = powers + shifter
                        
                    endif
                    
                enddo
                
                do i = -1, n-2
                    powers(i) = powers(i)/(i+2)
                enddo
                
                B_coefficents(j) = (-1.0)**(j+1)/(fact(j+1)*fact(n - 2 - j) )* sum(powers)
                
            enddo
        
        end function B_coefficents
        
        function fact(n)
            
            integer, intent(in) :: n
            integer :: i, fact
            
            fact = 1.0
            do i = 2, n
                fact = fact * i
            enddo

        end function fact

end module odu
