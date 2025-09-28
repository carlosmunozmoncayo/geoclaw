module exact_integrators
    implicit none
    double precision, parameter :: gammaEDC = 1.d0
    contains
    subroutine one_exact_step_EDC(h, hw, hp, csq_space, dt, hw2, hp2)
        implicit none
        double precision, intent(in) :: h, csq_space, dt
        double precision, intent(in) :: hw, hp
        double precision, intent(out) :: hw2, hp2
        !Locals
        double precision :: p,w,expdtA(2,2),result(2,1)
        double precision :: b, posc
        w = hw/h
        p = hp/h
        posc = dt*2.d0*csq_space/h
        b = dt*gammaEDC/h
        !Calculate the exact solution for the half problem:
        !(w,p)_t = h^{-1}(gamma*p,-2*csq*w) = A(p,w),
        !where A= [[0,gamma],[-2*csq,0]]/h.
        !The exact solution is given by:
        !(w,p)(dt) = exp(dt*A)*(w,p)(0).
        expdtA(1,1) = cos(dsqrt(b*posc)) !cosh(i*a) = cos(a)
        expdtA(1,2) = dsqrt(b/posc)*sin(dsqrt(b*posc)) !sinh(i*a) = i*sin(a)
        expdta(2,1) = -dsqrt(posc/b)*sin(dsqrt(b*posc)) !sin(i*a) = i*sinh(a)
        expdtA(2,2) = cos(dsqrt(b*posc)) !cosh(i*a) = cos(a)
        !Compute w,p
        w = expdtA(1,1)*w + expdtA(1,2)*p
        p = expdtA(2,1)*w + expdtA(2,2)*p
        !Output
        hw2 = w*h
        hp2 = p*h
    end subroutine one_exact_step_EDC

end module exact_integrators
