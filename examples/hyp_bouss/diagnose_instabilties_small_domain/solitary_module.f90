module solitary_module
    contains
!Solitary wave for Saint-Marie's system presented in Bristeau's paper
subroutine solitary_wave_Bristeau(x, h0, a, d, gamma, grav, t, h, u, p, w, x0, direction)
    implicit none
    double precision,intent(in) :: x, h0, a, d, gamma, grav, t,direction
    double precision,intent(out) :: h, u, p, w

    double precision, parameter :: pi = acos(-1.0d0)
    double precision :: sech, sechp, sechpp 
    double precision :: l, c0, L_star, x0, x_translated

    l = sqrt(((h0**3) / a + h0**2) * 2.0d0 / gamma)
    c0 = direction*sqrt(gamma / 2.0d0) * (l / d) * sqrt(grav * h0**3 / (0.5d0 * gamma * l**2 - h0**2))
    L_star = sqrt(4.0d0 / (3.0d0 * a)) * acosh(1.0d0 / 0.05d0)
    !x0 = L_star - 5.0d0 * c0
    x_translated = x - c0 * t - x0
    sech = 1.0d0 / cosh(x_translated / l)
    sechp = -sech * tanh(x_translated / l)
    sechpp = -sechp * tanh(x_translated / l) - sech**3
    h = h0 + a * sech**2
    u = c0 * (1.0d0 - d / h)
    w = -(a * c0 * d / (l * h)) * sech * sechp
    p = ((a * c0**2 * d**2) / (2.0d0 * l**2 * h**2)) * &
      (((2.0d0 * h0 - h) * sechp)**2 + h * sech * sechpp)
end subroutine solitary_wave_Bristeau
end module