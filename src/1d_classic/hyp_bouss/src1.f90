subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
    ! Called to update q by solving source term equation
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This default version integrates manning friction or other friction terms if present
    
    ! Version from 1D branch of GeoClaw originally from Dave George

    ! Also handles radial source term, if desired.
    ! Note: assumes radial about lower boundary, wall bc should be imposed

    
    use geoclaw_module, only: dry_tolerance, grav, DEG2RAD
    use geoclaw_module, only: friction_forcing
    use geoclaw_module, only: frictioncoeff => friction_coefficient
    use geoclaw_module, only: earth_radius, coordinate_system
    use grid_module, only: xcell
    use bouss_module, only: boussEquations, alpha, cm1, c01, cp1, &
                      solve_tridiag_ms, build_tridiag_sgn, &
                      solve_tridiag_sgn, useBouss
    !Parameters for Hyperbolic Relaxation EDC
    use bouss_module, only: c_sq => boussEDCcsq !Reference hyperbolic relaxation parameter
    use bouss_module, only: gammaEDC => boussEDCgamma !Approximate SGN (3/2) or Sainte-Marie equations (2)
    use bouss_module, only: S_BBBD => Source_HypRel_BBBD !Source term function for BBBD

    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc)

    ! Locals
    real(kind=8) :: eta(0:mx+1)
    real(kind=8) :: gamma, rcell, tanxR, etax
    real(kind=8) :: rk_stage(1:mx,4), delt
    integer ::  i,k,ii,rk_order
    real(kind=8)  q0(meqn,1-mbc:mx+mbc)
    real(kind=8) psi(mx+2)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local EDC Hyperbolic Relaxation
    !The 2 above should be passed by bouss_module
    double precision, dimension(4,2) :: rk_EDC
    double precision :: bath !Local bathymetry
    double precision :: csq_space !Local hyperbolic relaxation parameter
    double precision :: h, hp, hu, hw, u, p, w
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Local BBBD
    double precision :: hwout, hsgout, hpout, hpbout

    if (frictioncoeff.gt.0.d0 .and. friction_forcing) then
          ! integrate source term based on Manning formula
            do i=1,mx
               if (q(1,i)<=dry_tolerance) then
                  q(2,i) = 0.0
               else
                  gamma= dsqrt(q(2,i)**2)*(grav*frictioncoeff**2)/(q(1,i)**(7.0/3.0))
                  q(2,i)= q(2,i)/(1.d0 + dt*gamma)
              endif
            enddo
    endif

!      ----------------------------------------------------------------
!Forcing terms for radially symmetric equations
    if (coordinate_system == -1 .and. boussEquations>=0) then
        ! radial source term for SWE:
        do i=1,mx
            if (q(1,i) .gt. dry_tolerance) then
                ! x is radial coordinate in meters, x>=0, 
                ! u = radial velocity
                q(1,i) = q(1,i) - dt/xcell(i) * q(2,i)
                u = q(2,i)/q(1,i)
                q(2,i) = q(2,i) - dt/xcell(i) * q(1,i)*u**2
            end if
        end do
    end if

    if (coordinate_system == -1 .and. boussEquations==-1) then
        do i=1,mx
            if (q(1,i) .gt. dry_tolerance) then
                ! x is radial coordinate in meters, x>=0, 
                ! u = radial velocity
                h = q(1,i)
                q(1,i) = h - dt/xcell(i) * q(2,i)
                u = q(2,i)/h
                q(2,i) = q(2,i) - dt/xcell(i) * h*u**2
                if (useBouss(i)) then
                    q(3,i) = q(3,i) - dt/xcell(i) * u*q(3,i)
                    q(4,i) = q(4,i) - dt/xcell(i) * (u*q(4,i)+h*u*c_sq)
                end if
            end if
        enddo
    end if

    if (coordinate_system == -1 .and. boussEquations==-2) then
        do i=1,mx
            if (q(1,i) .gt. dry_tolerance) then
                ! x is radial coordinate in meters, x>=0, 
                ! u = radial velocity
                h = q(1,i)
                q(1,i) = h - dt/xcell(i) * q(2,i)
                u = q(2,i)/h
                q(2,i) = q(2,i) - dt/xcell(i) * h*u**2
                if (useBouss(i)) then
                    q(3,i) = q(3,i) - dt/xcell(i) * u*q(3,i)
                    q(4,i) = q(4,i) - dt/xcell(i) * u*q(4,i)
                    q(5,i) = q(5,i) - dt/xcell(i) * (u*q(5,i)+h*u*c_sq)
                    q(6,i) = q(6,i) - dt/xcell(i) * u*q(6,i)
                end if
            end if
        enddo
    end if

!      ----------------------------------------------------------------

    if (coordinate_system == 2) then
        ! source term for x = latitude in degrees -90 <= x <= 90,
        ! u = velocity in latitude direction (m/s) on sphere:
        do i=1,mx
            if (q(1,i) .gt. dry_tolerance) then
                tanxR = tan(xcell(i)*DEG2RAD) / earth_radius
                q(1,i) = q(1,i) + dt * tanxR * q(2,i)
                u = q(2,i)/q(1,i)
                q(2,i) = q(2,i) + dt * tanxR * q(1,i)*u**2
            endif
         enddo
     endif

    ! -------------------------------------------------
    if (boussEquations > 0) then
    ! Boussinesq terms for MS and SGN 

        rk_order = 1  ! 1 for Forward Euler, 2 for second-order RK
        delt = dt / rk_order

        if (boussEquations == 2) then
            ! for SGN, need to factor matrix each step
            call build_tridiag_sgn(meqn,mbc,mx,xlower,dx,q,maux,aux)
        endif
                
        q0  = q
          
        do i=0,mx+1
           if (q(1,i) > dry_tolerance) then
               eta(i) = q(1,i)+aux(1,i)
           else
               eta(i) = 0.d0
           endif
        enddo

        !-----------------------
        ! First stage (only stage for rk_order == 1):
        
        if (boussEquations == 1) then
          call solve_tridiag_ms(mx,meqn,mbc,dx,q0,maux,aux,psi)
          ! returns solution psi = source term for Madsen
        else if (boussEquations == 2) then
          call solve_tridiag_sgn(mx,meqn,mbc,dx,q0,maux,aux,psi)
          ! modify solution psi for source term of SGN:
          do i=1,mx
              if (useBouss(i)) then
                  etax = cm1(i)*eta(i-1) + c01(i)*eta(i) + cp1(i)*eta(i+1)
                  psi(i+1) = q(1,i) * (grav/alpha * etax - psi(i+1))
              else
                  psi(i+1) = 0.d0
              endif
           enddo
        endif
              
        ! Forward Euler update to momentum q(2,:):
        ! Note psi(1) used for BC, so psi(i+1) updates q(2,i):
        q0(2,1:mx) = q0(2,1:mx) + delt*psi(2:mx+1)
          

        !-----------------------
        if (rk_order == 1) then
            q = q0 ! and we are done
              
        else if (rk_order == 2) then

            ! Second stage for 2-stage R-K, solve for psi at midpoint
            ! in time based on q0 computed in first stage:
            
            if (boussEquations == 1) then
                call solve_tridiag_ms(mx,meqn,mbc,dx,q0,maux,aux,psi)
                ! returns solution psi = source term for Madsen
            else if (boussEquations == 2) then
                call solve_tridiag_sgn(mx,meqn,mbc,dx,q0,maux,aux,psi)
                ! modify solution psi for source term of SGN:
                do i=1,mx
                    ! note that h,eta were not changed by first stage
                    if (useBouss(i)) then
                        etax = cm1(i)*eta(i-1) + c01(i)*eta(i) + cp1(i)*eta(i+1)
                        psi(i+1) = q(1,i) * (grav/alpha * etax - psi(i+1))
                    else
                        psi(i+1) = 0.d0
                    endif
                enddo
            endif
              
            ! Second stage is midpoint method dt=delt*2:
            ! Note psi(1) used for BC, so psi(i+1) updates q(2,i):
            q(2,1:mx) = q(2,1:mx) + 2.d0*delt*psi(2:mx+1)

        endif ! rk_order == 2

    endif ! end of Bouss terms

    ! -------------------------------------------------
    if (boussEquations == -1) then
    !Boussinesq terms for HypRel EDC
        do i=2-mbc,mx+mbc-1 !i=1,mx+mbc
            if (useBouss(i)) then
                call one_exact_step_EDC(q(1,i), q(3,i), q(4,i), c_sq, dt)
            endif
        end do
    end if
    if (boussEquations == -2) then
    !Boussinesq terms for HypRel BBBD 
        do i=2-mbc,mx+mbc-1 !i=1,mx+mbc
            if (useBouss(i)) then
                ! call one_exact_step_BBBD(q(1,i), q(3,i), q(4,i), q(5,i), q(6,i), c_sq, dt)
                call one_step_RKGL4_BBBD(q(1,i), q(3,i), q(4,i), q(5,i), q(6,i), c_sq, dt)
            endif
        end do
    end if

    contains
        !---------------------------------------------------------------------
        !Exact integration of Source terms for EDC and BBBD systems
        !---------------------------------------------------------------------
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
            implicit none
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
            implicit none
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
            A(3,2) = -csq_space/h
            A(4,1) = -6.d0*csq_space/h
            A(4,2) = 3.d0*csq_space/h
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
        !---------------------------------------------------------------------
        !---------------------------------------------------------------------

end subroutine src1

