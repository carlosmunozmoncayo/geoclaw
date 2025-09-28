module root_finding
    implicit none

    contains
subroutine newtons_metod(f,jacobian,x0,x,max_iter)
    !Root finding for a 2x2 system of equations using Newton's method
    implicit none
    !Input and output
    double precision, intent(in) :: x0(2)
    integer, intent(in) :: max_iter
    double precision, intent(out) :: x(2)
    interface
        function f(x) result(y)
            double precision, intent(in) :: x(2)
            double precision :: y(2)
        end function f
        function jacobian(x) result(J)
            double precision, intent(in) :: x(2)
            double precision :: J(2,2)
        end function jacobian
    end interface

    !Local variables
    integer :: i
    double precision :: J(2,2), invJ(2,2), dx(2), x1(2), x2(2)
    double precision, parameter :: tol = 1.0e-9

    !Newton's method
    x1 = x0
    do i = 1, max_iter
        !Calculate the Jacobian
        J = jacobian(x1)
        !Invert the Jacobian (exactly since it is 2x2)
        invJ(1,1) = J(2,2)
        invJ(2,2) = J(1,1)
        invJ(1,2) = -J(1,2)
        invJ(2,1) = -J(2,1)
        invJ = 1.d0/(J(1,1)*J(2,2) - J(1,2)*J(2,1)) * invJ
        !Calculate the step dx
        dx = matmul(invJ,-f(x1))
        !Update x2
        x2 = x1 + dx
        !Check for convergence
        if (maxval(abs(f(x2))) < tol) exit
        !Update x1
        x1 = x2
    end do
    !Output
    x = x2
end subroutine newtons_metod

end module root_finding