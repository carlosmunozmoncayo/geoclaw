! This routine should be a simplified version of src2
! which applies source terms for a 1-d slice of data along the
! edge of a grid.  This is called only from qad where the conservative
! fix-up is applied and is used to apply source terms over partial
! time steps to the coarse grid cell values used in solving Riemann
! problems at the interface between coarse and fine grids.
subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)
      
    use geoclaw_module, only: g => grav, rho, coriolis_forcing, coriolis
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: omega, coordinate_system, manning_coefficient
    use geoclaw_module, only: manning_break, num_manning, dry_tolerance, rho_air
    
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

    ! Input
    integer, intent(in) :: meqn, mbc, mx1d, maux
    real(kind=8), intent(in) :: t, dt
    real(kind=8), intent(inout) :: q1d(meqn, mx1d), aux1d(maux, mx1d)

    ! Local storage
    integer :: i, nman
    logical :: found
    real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2), coeff, tau
    real(kind=8) :: wind_speed, theta, P_atmos_x, P_atmos_y

    ! Algorithm parameters
    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

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
    
    ! Friction forcing
    if (friction_forcing) then

        do i=1,mx1d

            ! Extract depths
            h = q1d(1,i)
            hu = q1d(2,i)
            hv = q1d(3,i)

            ! If depth is near-zero, set momentum to zero
            if (h < depth_tolerance) then
                q1d(2:3,i) = 0.d0
                cycle 
            endif
            
            ! Apply friction source term only if in shallower water
            if (h <= friction_depth) then
                if (.not.variable_friction) then
                    do nman = num_manning, 1, -1
                        if (aux1d(1,i) .lt. manning_break(nman)) then
                            coeff = manning_coefficient(nman)
                        endif
                    enddo
                else
                    coeff = aux1d(friction_index, i)
                end if

                ! Calculate source term
                gamma = sqrt(hu**2 + hv**2) * (g * coeff**2) / h**(7.d0/3.d0)
                dgamma = 1.d0 + dt * gamma
                q1d(2, i) = q1d(2, i) / dgamma
                q1d(3, i) = q1d(3, i) / dgamma
            endif
        enddo
    endif
    
    ! Only lat-long coordinate system supported here right now
    if (coriolis_forcing .and. coordinate_system == 2) then

        do i=1,mx1d
            ! aux(3,:,:) stores the y coordinates multiplied by deg2rad
            fdt = 2.d0 * omega * sin(aux1d(3,i)) * dt

            ! Calculate matrix components
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)
    
            q1d(2,i) = q1d(2,i) * a(1,1) + q1d(3,i) * a(1,2)
            q1d(3,i) = q1d(2,i) * a(2,1) + q1d(3,i) * a(2,2)
        enddo
    endif

    ! = Wind Forcing =========================================================
    if (wind_forcing) then
        ! Cannot use sector based wind mappings here due to lack of information
        ! sloc = storm_location(t)
        ! theta = storm_direction(t)
        do i=1,mx1d
            if (q1d(1,i) > dry_tolerance) then
                wind_speed = sqrt(aux1d(wind_index,i)**2 &
                                + aux1d(wind_index+1,i)**2)
                tau = wind_drag(wind_speed,theta) * rho_air * wind_speed / rho(1)
                q1d(2,i) = q1d(2,i) + dt * tau * aux1d(wind_index,i)
                q1d(3,i) = q1d(3,i) + dt * tau * aux1d(wind_index+1,i)
            endif
        enddo
    endif
    ! ========================================================================
    
    ! == Pressure Forcing ====================================================
    ! Handled in Riemann solver

    ! Need to add dx and dy to calling signature from qad in order for this to 
    ! work, can probably get it from amr_module but need to know which grid we
    ! are working on
    ! if (.false.) then
    !     stop "Not sure how to proceed, need direction and the right dx or dy."
    ! endif

        if (boussEquations == -1) then
        do i=1,mx1d
                if (aux1d(dr_index,i)>small_tol) then
                    h = q1d(1,i)
                    csq_space = aux1d(csq_index,i)
                    call one_step_BE_EDC(h, q1d(4,i), q1d(5,i), csq_space, dt)
                end if
        enddo
    endif

    if (boussEquations == -2) then
        do i=1,mx1d
                if (aux1d(dr_index,i)>small_tol) then
                    csq_space = aux1d(csq_index,i)
                    call one_exact_step_BBBD(q1d(1,i), q1d(4,i), q1d(5,i), q1d(6,i), q1d(7,i),&
                                     csq_space, dt)
                end if
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

end subroutine src1d
