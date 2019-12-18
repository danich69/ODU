module diffur

    implicit none
    
    integer, parameter :: dims = 2
    integer, parameter :: interpolate_order = 4
    integer, parameter :: extrapolate_order = 4
    integer :: i
    
    real(8), parameter :: interval = 10
    real(8), dimension(dims), parameter :: x_0 = (/(1.0, i=1, dims)/)
    real(8), parameter :: h = 1e-6

    contains
    
        function f(t, x)
        
            real(8), dimension(:) :: x
            real(8), dimension(size(x)) :: f
            real(8) :: t
            
             f(1) = x(1)
             f(2) = x(2)
        
        end function f

end module diffur
