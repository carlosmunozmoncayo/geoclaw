subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
      
    use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
    use geoclaw_module, only: sea_level
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: manning_coefficient
    use geoclaw_module, only: manning_break, num_manning
    use geoclaw_module, only: spherical_distance, coordinate_system
    use geoclaw_module, only: RAD2DEG, pi, dry_tolerance, DEG2RAD
    use geoclaw_module, only: rho_air
    use geoclaw_module, only: earth_radius, sphere_source
      
    use storm_module, only: wind_forcing, pressure_forcing, wind_drag
    use storm_module, only: wind_index, pressure_index
    use storm_module, only: storm_direction, storm_location

    use friction_module, only: variable_friction, friction_index

    !Parameters for Hyperbolic Relaxation 
    use bouss_module, only: boussEquations
    use bouss_module, only: S_EDC => Source_HypRel_EDC !Source term function
    use bouss_module, only: S_BBBD => Source_HypRel_BBBD !Source term function
    use bouss_module, only: c_sq => boussEDCcsq !Reference hyperbolic relaxation parameter
    use bouss_module, only: gammaEDC => boussEDCgamma !Approximate SGN (3/2) or Sainte-Marie equations (2)
    use bouss_module, only: csq_index => bouss_csq_index !Index of the hyperbolic relaxation parameter
    use bouss_module, only: dr_index => boussDecayRate_index !Index of the decay rate
    use bouss_module, only: small_tol=> small_tol_bouss
    use bouss_module, only: large_tol=> large_tol_bouss

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j, nman
    real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2), coeff
    real(kind=8) :: xm, xc, xp, ym, yc, yp, dx_meters, dy_meters
    real(kind=8) :: u, v, hu0, hv0
    real(kind=8) :: tau, wind_speed, theta, phi, psi, P_gradient(2), S(2)
    real(kind=8) :: Ddt, sloc(2)
    real(kind=8) :: tanyR, huv, huu, hvv
    real(kind=8) :: eta

    ! Algorithm parameters

    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    ! Physics
    ! Nominal density of water
    real(kind=8), parameter :: rho = 1025.d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Local EDC Hyperbolic Relaxation
    double precision :: bath !Local bathymetry
    double precision :: csq_space !Local hyperbolic relaxation parameter
    double precision :: hp, hw
    ! Local BBBD Hyperbolic Relaxation
    double precision :: hsg, hpb
    ! Source term substepping
    double precision :: dt_LSC !Linear stability condition for the source term
    double precision :: t_substep, dt_remainder
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ----------------------------------------------------------------
    ! Spherical geometry source term(s)
    !
    ! These should be included for shallow water on the sphere, 
    ! at least in the mass term, but were only added in v5.9.0 as an option.
    ! rundata.geo_data.sphere_source can now be set in setrun.py
    ! Set sphere_source = 0 to omit source terms for backward compatibility
    ! sphere_source = 1 should become the default?

    if ((coordinate_system == 2) .and. (sphere_source > 0)) then
        ! add in spherical source term in mass equation 
        ! if sphere_source in [1,2],
        ! and also in momentum equations if sphere_source == 2
        do j=1,my
            y = ylower + (j - 0.5d0) * dy
            tanyR = tan(y*DEG2RAD) / earth_radius
            do i=1,mx
                if (q(1,i,j) > dry_tolerance) then
                    ! source term in mass equation:
                    q(1,i,j) = q(1,i,j) + dt * tanyR * q(3,i,j)

                    if (sphere_source == 2) then
                        ! Momentum source terms that drop out if linearized:
                        ! These seem to have very little effect for
                        ! practical problems
                        huv = q(2,i,j)*q(3,i,j)/q(1,i,j)
                        huu = q(2,i,j)*q(2,i,j)/q(1,i,j)
                        hvv = q(3,i,j)*q(3,i,j)/q(1,i,j)
                        q(2,i,j) = q(2,i,j) + dt * tanyR * 2.d0*huv
                        q(3,i,j) = q(3,i,j) + dt * tanyR * (hvv - huu)
                    endif
                endif
            enddo
        enddo
    endif
                
    ! ----------------------------------------------------------------                
    ! Friction source term
    if (friction_forcing) then
        do j=1,my
            do i=1,mx
                ! Extract appropriate momentum
                if (q(1,i,j) < depth_tolerance) then
                    q(2:3,i,j) = 0.d0
                else
                    ! Apply friction source term only if in shallower water
                    if (q(1,i,j) <= friction_depth) then
                        if (.not.variable_friction) then
                            do nman = num_manning, 1, -1
                                if (aux(1,i,j) .lt. manning_break(nman)) then
                                    coeff = manning_coefficient(nman)
                                endif
                            enddo
                        else
                            coeff = aux(friction_index,i,j)
                        endif
                        
                        ! Calculate source term
                        gamma = sqrt(q(2,i,j)**2 + q(3,i,j)**2) * g     &   
                              * coeff**2 / (q(1,i,j)**(7.d0/3.d0))
                        dgamma = 1.d0 + dt * gamma
                        q(2, i, j) = q(2, i, j) / dgamma
                        q(3, i, j) = q(3, i, j) / dgamma
                    endif
                endif
            enddo
        enddo
    endif
    ! End of friction source term

    ! Coriolis source term
    ! TODO: May want to remove the internal calls to coriolis as this could 
    !       lead to slow downs.
    if (coriolis_forcing) then
        do j=1,my
            y = ylower + (j - 0.5d0) * dy
            fdt = coriolis(y) * dt ! Calculate f dependent on coordinate system

            ! Calculate matrix components
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)

            do i=1,mx
                q(2,i,j) = q(2, i, j) * a(1,1) + q(3, i, j) * a(1,2)
                q(3,i,j) = q(2, i, j) * a(2,1) + q(3, i, j) * a(2,2)
            enddo
        enddo
    endif
    ! End of coriolis source term

    ! wind -----------------------------------------------------------
    if (wind_forcing) then
        ! Need storm location and direction for sector based wind drag
        sloc = storm_location(t)
        theta = storm_direction(t)
        do j=1,my
            yc = ylower + (j - 0.5d0) * dy
            do i=1,mx
                xc = xlower + (i - 0.5d0) * dx
                if (q(1,i,j) > dry_tolerance) then
                    psi = atan2(yc - sloc(2), xc - sloc(1))
                    if (theta > psi) then
                        phi = (2.d0 * pi - theta + psi) * RAD2DEG
                    else
                        phi = (psi - theta) * RAD2DEG 
                    endif
                    wind_speed = sqrt(aux(wind_index,i,j)**2        &
                                    + aux(wind_index+1,i,j)**2)
                    tau = wind_drag(wind_speed, phi) * rho_air * wind_speed / rho
                    q(2,i,j) = q(2,i,j) + dt * tau * aux(wind_index,i,j)
                    q(3,i,j) = q(3,i,j) + dt * tau * aux(wind_index+1,i,j)
                endif
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------

    ! Atmosphere Pressure --------------------------------------------
    ! Handled in Riemann solver
    ! if (pressure_forcing) then
    !     do j=1,my  
    !         ym = ylower + (j - 1.d0) * dy
    !         yc = ylower + (j - 0.5d0) * dy
    !         yp = ylower + j * dy
    !         do i=1,mx  
    !             xm = xlower + (i - 1.d0) * dx
    !             xc = xlower + (i - 0.5d0) * dx
    !             xp = xlower + i * dx
                
    !             if (coordinate_system == 2) then
    !                 ! Convert distance in lat-long to meters
    !                 dx_meters = spherical_distance(xp,yc,xm,yc)
    !                 dy_meters = spherical_distance(xc,yp,xc,ym)
    !             else
    !                 dx_meters = dx
    !                 dy_meters = dy
    !             endif

    !             ! Extract depths
    !             h = q(1,i,j)

    !             ! Calculate gradient of Pressure
    !             P_gradient(1) = (aux(pressure_index,i+1,j) &
    !                            - aux(pressure_index,i-1,j)) / (2.d0 * dx_meters)
    !             P_gradient(2) = (aux(pressure_index,i,j+1) &
    !                            - aux(pressure_index,i,j-1)) / (2.d0 * dy_meters)

    !                 ! Modify momentum in each layer
    !             if (h > dry_tolerance) then
    !                 q(2, i, j) = q(2, i, j) - dt * h * P_gradient(1) / rho
    !                 q(3, i, j) = q(3, i, j) - dt * h * P_gradient(2) / rho
    !             end if
    !         enddo
    !     enddo
    ! endif

    if (boussEquations == -1) then
        do j=1,my
            do i=1,mx
                if (aux(dr_index,i,j)>small_tol) then
                    h = q(1,i,j)
                    csq_space = aux(csq_index,i,j) 
                    ! call one_exact_step_EDC(h, q(4,i,j), q(5,i,j), csq_space, dt)
                    ! call one_step_RK4_EDC(h, q(4,i,j), q(5,i,j), csq_space, dt)
                    call one_step_BE_EDC(h, q(4,i,j), q(5,i,j), csq_space, dt)
                end if
            enddo
        enddo
    endif

    if (boussEquations == -2) then
        do j=1,my
            do i=1,mx
                if (aux(dr_index,i,j)>small_tol) then
                    csq_space = aux(csq_index,i,j)
                    ! call one_step_RK4_BBBD(q(1,i,j), q(4,i,j), q(5,i,j), q(6,i,j), q(7,i,j),&
                    !                  csq_space, dt)
                    ! call one_step_RKGL4_BBBD(q(1,i,j), q(4,i,j), q(5,i,j), q(6,i,j), q(7,i,j),&
                    !                  csq_space, dt)
                    call one_exact_step_BBBD(q(1,i,j), q(4,i,j), q(5,i,j), q(6,i,j), q(7,i,j),&
                                     csq_space, dt)
                end if
            enddo
        enddo
    endif

    contains
        subroutine one_step_RK4_EDC(h, hw, hp, csq_space, dt)
            double precision, intent(in) :: h, csq_space, dt
            double precision, intent(inout) :: hw, hp
            double precision :: rk_EDC(4,2)
            rk_EDC(1,:) = S_EDC(hw,hp,h, csq_space)
            rk_EDC(2,:) = S_EDC(hw+0.5d0*dt*rk_EDC(1,1),hp+0.5d0*dt*rk_EDC(1,2), h, csq_space)
            rk_EDC(3,:) = S_EDC(hw+0.5d0*dt*rk_EDC(2,1),hp+0.5d0*dt*rk_EDC(2,2), h, csq_space)
            rk_EDC(4,:) = S_EDC(hw+dt*rk_EDC(3,1),hp+dt*rk_EDC(3,2), h, csq_space)
            hw = hw + (1.d0/6.d0)*dt*(rk_EDC(1,1)+2.d0*rk_EDC(2,1)+2.d0*rk_EDC(3,1)+rk_EDC(4,1))
            hp = hp + (1.d0/6.d0)*dt*(rk_EDC(1,2)+2.d0*rk_EDC(2,2)+2.d0*rk_EDC(3,2)+rk_EDC(4,2))
        end subroutine one_step_RK4_EDC

        subroutine one_step_BE_EDC(h, hw, hp, csq_space, dt)
            !A backward Euler step for the EDC equations
            implicit none
            double precision, intent(in) :: h, csq_space, dt
            double precision, intent(inout) :: hw, hp
            !Locals
            double precision :: p,w
            double precision ::denom,b,c
            w = hw/h
            p = hp/h
            denom = 1.d0 + (dt*dt*2.d0*csq_space*gammaEDC)/(h*h)
            b = dt*gammaEDC/h
            c = -2.d0*dt*csq_space/h
            !Output (Inverting the matrix explicitly and multiplying back h)
            hw = h*(w + b*p)/denom 
            hp = h*(c*w + p)/denom
        end subroutine one_step_BE_EDC

        subroutine one_exact_step_EDC(h, hw, hp, csq_space, dt)
            implicit none
            double precision, intent(in) :: h, csq_space, dt
            double precision, intent(inout) :: hw, hp
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
            result(:,1) = (/w,p/)
            result = matmul(expdtA,result)
            !Compute w,p
            ! w = expdtA(1,1)*w + expdtA(1,2)*p
            ! p = expdtA(2,1)*w + expdtA(2,2)*p
            !Output
            hw = result(1,1)*h
            hp = result(2,1)*h
        end subroutine one_exact_step_EDC
        
        subroutine one_step_RK4_BBBD(h, hw, hsg, hp, hpb, csq_space, dt)
            double precision, intent(in) :: h, csq_space, dt
            double precision, intent(inout) :: hw, hsg, hp, hpb
            double precision :: rk_BBBD(4,4)
            rk_BBBD(1,:) = S_BBBD(h,hw,hsg,hp,hpb,csq_space)
            rk_BBBD(2,:) = S_BBBD(h,hw+0.5d0*dt*rk_BBBD(1,1),hsg+0.5d0*dt*rk_BBBD(1,2), &
                                hp+0.5d0*dt*rk_BBBD(1,3),hpb+0.5d0*dt*rk_BBBD(1,4),&
                                csq_space)
            rk_BBBD(3,:) = S_BBBD(h,hw+0.5d0*dt*rk_BBBD(2,1),hsg+0.5d0*dt*rk_BBBD(2,2), &
                                hp+0.5d0*dt*rk_BBBD(2,3),hpb+0.5d0*dt*rk_BBBD(2,4),&
                                csq_space)
            rk_BBBD(4,:) = S_BBBD(h,hw+dt*rk_BBBD(3,1),hsg+dt*rk_BBBD(3,2), &
                                hp+dt*rk_BBBD(3,3),hpb+dt*rk_BBBD(3,4),&
                                csq_space)
            hw = hw + (1.d0/6.d0)*dt*(rk_BBBD(1,1)+2.d0*rk_BBBD(2,1)+ &
                                2.d0*rk_BBBD(3,1)+rk_BBBD(4,1))
            hsg = hsg + (1.d0/6.d0)*dt*(rk_BBBD(1,2)+2.d0*rk_BBBD(2,2)+ &
                                2.d0*rk_BBBD(3,2)+rk_BBBD(4,2))
            hp = hp + (1.d0/6.d0)*dt*(rk_BBBD(1,3)+ 2.d0*rk_BBBD(2,3)+ &
                                2.d0*rk_BBBD(3,3)+rk_BBBD(4,3))
            hpb = hpb + (1.d0/6.d0)*dt*(rk_BBBD(1,4)+2.d0*rk_BBBD(2,4)+ &
                                2.d0*rk_BBBD(3,4)+rk_BBBD(4,4))
        end subroutine one_step_RK4_BBBD 

        subroutine one_step_RKGL4_BBBD(h, hw, hsg, hp, hpb, csq_space, dt)
            !Take one step using the Gauss-Legendre 2 stage 4th order method
            !it is an implicit method, so we need some matrix inversion
            !Intended to solve q_t = A(q) with A(q) = A*q
            !And q = (w,sg,p,pb), 
            implicit none
            double precision, intent(in) :: h, csq_space, dt
            double precision, intent(inout) :: hw, hsg, hp, hpb
            !Matrix handling
            double precision, dimension(4,4) :: A, B1,B2,invB1,invB2,AinvB2A
            double precision :: q(4,1)
            double precision :: rhs1(4,1), rhs2(4,1)
            !Butcher tableau
            double precision :: a11, a12, a21, a22, bRK1, bRK2
            !Stages
            double precision :: k1(4,1), k2(4,1)
            !Locals
            integer :: i

            q(:,1) = (/hw,hsg,hp,hpb/)
            q = q/h
            A = 0.d0
            A(1,4) = 1.d0/h
            A(2,3) = 12.d0/h
            A(2,4) = -6.d0/h
            A(3,2) = -2.d0*csq_space/h
            A(4,1) = -6.d0*csq_space/h
            A(4,2) = 3.d0/h
            !Butcher tableau
            a11 = 0.25d0
            a12 = 0.25d0 -dsqrt(3.d0)/6.d0
            a21 = 0.25d0 + dsqrt(3.d0)/6.d0
            a22 = 0.25d0
            bRK1 = 0.5d0
            bRK2 = 0.5d0
            !B2 = I-dt*a22*A
            B2 = 0.d0
            do i=1,4
                B2(i,i) = 1.d0
            end do
            B2 = B2 - dt*a22*A
            invB2 = matinv4(B2)
            AinvB2A = matmul(A,matmul(invB2,A))
            
            !B1 = I-dt*a11*A
            B1 = 0.d0
            do i=1,4
                B1(i,i) = 1.d0
            end do
            B1 = B1 - dt*a11*A

            !RHS1
            rhs1 = matmul(A,q)+dt*a12*matmul(AinvB2A,q)
            invB1 = matinv4(B1-dt*dt*a12*a21*AinvB2A)
            k1 = matmul(invB1,rhs1)
            !RHS2
            rhs2 = matmul(A,q+dt*a21*k1)
            k2 = matmul(invB2,rhs2)
            !Update
            q = q + dt*(bRK1*k1+bRK2*k2)
            hw = q(1,1)*h
            hsg = q(2,1)*h
            hp = q(3,1)*h
            hpb = q(4,1)*h
        end subroutine one_step_RKGL4_BBBD


        subroutine one_exact_step_BBBD(h, hw, hsg, hp, hpb, csq_space, dt)
            implicit none

            double precision, intent(in) :: h, csq_space, dt
            double precision, intent(inout) :: hw, hsg, hp, hpb
            !Locals
            complex(kind=8) :: q(4,1)
            complex(kind=8) :: R(4,4), Diag(4), Rinv(4,4)
            complex(kind=8) :: z,sq7,a1,a2,a3,a4,a5,a6,a7,a8,c
            complex(kind=8) :: img = cmplx(0.d0,1.d0)
            integer :: i,j

            c = sqrt(csq_space)
            sq7 = sqrt(7.d0)
            a1 = 3.d0 + sq7
            a2 = -3.d0 + sq7
            a3 = 2.d0 + sq7
            a4 = -2.d0 + sq7
            a5 = 1.d0 - sq7
            a6 = 1.d0 + sq7
            a7 = sqrt(2.d0 - 2.d0*sq7/3.d0)
            a8 = sqrt(1.d0 - sq7/3.d0)

            R(1,:) = (1.d0/sqrt(6.d0))*(/img/sqrt(a1),-img/sqrt(a1),img/sqrt(-a2),1.d0/sqrt(a2)/)
            R(2,:) = (/ -img*a3*a8, img*a3*a8, -img*a4*a7/a2, img*a4*a7/a2 /)
            R(3,:) = (c/6.d0)*(/ a5, a5, a6, a6 /)
            R(4,:) = (/ c,c,c,c/)

            Rinv = complmatinv4(R)
            Diag = (dt*sqrt(6.d0)*c/h)*(/ -img*sqrt(a1),img*sqrt(a1),-sqrt(a2),sqrt(a2) /)

            !Initial conditions
            q(:,1) = (/hw,hsg,hp,hpb/)/h
            !More efficient matmul
            q = matmul(Rinv,q)
            do i=1,4
                q(i,1)=exp(Diag(i))*q(i,1)
            end do
            q = matmul(R,q)
            hw  = RealPart(q(1,1))*h
            hsg = RealPart(q(2,1))*h
            hp  = RealPart(q(3,1))*h
            hpb = RealPart(q(4,1))*h
        end subroutine one_exact_step_BBBD

        pure function matinv4(A) result(B)
            !! Performs a direct calculation of the inverse of a 4x4 matrix.
            double precision, intent(in) :: A(4,4)   !! Matrix
            double precision            :: B(4,4)   !! Inverse matrix
            double precision            :: detinv
        
            ! Calculate the inverse determinant of the matrix
            detinv = &
            1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)&
            -A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
            - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)&
            -A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
            + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)&
            -A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
            - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)&
            -A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))
        
            ! Calculate the inverse of the matrix
            B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))&
            +A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
            B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))&
            +A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
            B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))&
            +A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
            B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))&
            +A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
            B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))&
            +A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
            B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))&
            +A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
            B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))&
            +A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
            B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))&
            +A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
            B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))&
            +A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
            B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))&
            +A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
            B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))&
            +A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
            B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))&
            +A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
            B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))&
            +A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
            B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))&
            +A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
            B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))&
            +A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
            B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))&
            +A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
        end function
        pure function complmatinv4(A) result(B)
            !! Performs a direct calculation of the inverse of a 4x4 matrix.
            complex(kind=8), intent(in) :: A(4,4)   !! Matrix
            complex(kind=8)            :: B(4,4)   !! Inverse matrix
            complex(kind=8)          :: detinv
        
            ! Calculate the inverse determinant of the matrix
            detinv = &
            1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)&
            -A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
            - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)&
            -A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
            + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)&
            -A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
            - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)&
            -A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))
        
            ! Calculate the inverse of the matrix
            B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))&
            +A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
            B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))&
            +A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
            B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))&
            +A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
            B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))&
            +A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
            B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))&
            +A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
            B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))&
            +A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
            B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))&
            +A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
            B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))&
            +A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
            B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))&
            +A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
            B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))&
            +A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
            B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))&
            +A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
            B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))&
            +A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
            B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))&
            +A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
            B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))&
            +A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
            B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))&
            +A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
            B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))&
            +A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
        end function



end subroutine src2
