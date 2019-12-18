program in_out
	use diffur
    use odu
    implicit none
    	
!	-----------------------------------------------------
!	|	Program integrates differential equation   		|
!	|	Available methods are:              			|
!	|			Runge-Kutta - rk key					|
!	|			Adams interpolation - ai key			|
!	|			Adams extrapolation - ae key			|
!	-----------------------------------------------------

    real(8), dimension(1:int(interval/h + 1), size(x_0)) :: x
    real(8) :: t
    character*2 :: method
    character*6 :: file_name
    
    call getarg(1, method)
    write (file_name,"(a2,'.dat')") method
    
    select case(method)
        case('rk')
            call Runge_Kutta(x)
        case('ae')
            call Adams_extrapolation(x)
        case('ai')
            call Adams_interpolation(x)
        case default
            write(*,*) 'Inkorrect method! Please enter rk (Runge-Kutta)), ae (Adams extrapolation) or ai (Adams interpolation)'
            read(*,*) method
    end select
    
    t = 0
    
    open(1, file = trim(file_name))
    do i = 1, int(interval/h + 1)
        
        write(1,*) t, x(i, :)
        t = t + h
    
    enddo    
    close(1)

end program in_out
