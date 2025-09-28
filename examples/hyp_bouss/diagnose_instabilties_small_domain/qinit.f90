
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use geoclaw_module, only: sea_level
    use amr_module, only: t0
    use qinit_module, only: qinit_type,add_perturbation
    use qinit_module, only: variable_eta_init
    use qinit_module, only: force_dry,use_force_dry,mx_fdry, my_fdry
    use qinit_module, only: xlow_fdry, ylow_fdry, xhi_fdry, yhi_fdry
    use qinit_module, only: dx_fdry, dy_fdry
    use qinit_module, only: tend_force_dry


    use geoclaw_module, only: grav
    use geoclaw_module, only: dry_tolerance
    use bouss_module, only: boussEquations
    use bouss_module, only : c_sq => boussEDCcsq !Reference hyperbolic relaxation parameter
    use bouss_module, only: bouss_csq_index
    use bouss_module, only: gamma => boussEDCgamma
    ! use solitary_module, only: solitary_wave_Bristeau
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j,m, ii,jj
    real(kind=8) :: x,y
    real(kind=8) :: veta(1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8) :: ddxy
    real(kind=8) :: radius, uhat, u, v, p, w, x0,y0, direction
    real(kind=8) :: h, Hcap, A, d
    real(kind=8) :: xcell, ycell, inner_radius
    real(kind=8) :: scale = 125000.0d0*10.0d0
    real(kind=8) :: random_noise1, random_noise2
    
    q = 0.d0   ! initialize all elements to 0
    
    if (variable_eta_init) then
        ! Set initial surface eta based on eta_init
        call set_eta_init(mbc,mx,my,xlower,ylower,dx,dy,t0,veta)
      else
        veta = sea_level  ! same value everywhere
      endif

    forall(i=1:mx, j=1:my)
        q(1,i,j) = max(0.d0, veta(i,j) - aux(1,i,j))
    end forall

    if (use_force_dry .and. (t0 <= tend_force_dry)) then
     ! only use the force_dry if it specified on a grid that matches the 
     ! resolution of this patch, since we only check the cell center:
     ddxy = max(abs(dx-dx_fdry), abs(dy-dy_fdry))
     if (ddxy < 0.01d0*min(dx_fdry,dy_fdry)) then
       do i=1,mx
          x = xlower + (i-0.5d0)*dx
          ii = int((x - xlow_fdry + 1d-7) / dx_fdry)
          do j=1,my
              y = ylower + (j-0.5d0)*dy
              jj = int((y - ylow_fdry + 1d-7) / dy_fdry)
              jj = my_fdry - jj  ! since index 1 corresponds to north edge
              if ((ii>=1) .and. (ii<=mx_fdry) .and. &
                  (jj>=1) .and. (jj<=my_fdry)) then
                  ! grid cell lies in region covered by force_dry,
                  ! check if this cell is forced to be dry 
                  ! Otherwise don't change value set above:                  
                  if (force_dry(ii,jj) == 1) then
                      q(1,i,j) = 0.d0
                      endif
                  endif
          enddo ! loop on j
       enddo ! loop on i
       endif ! dx and dy agree with dx_fdry, dy_fdry
    endif ! use_force_dry

    
    ! Add perturbation to initial conditions
    if (qinit_type > 0) then
        call add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    endif

    !Generate a radial solitary wave
    Hcap = 4.d3
    A = 1.d0
    d = Hcap
    x0 = -90.6!-100.d0
    y0 = -1.2!-40.d0
    direction = 1.d0
    inner_radius = 0.05*scale
    do i=1,mx
        xcell = xlower + (i-0.5d0)*dx
        do j=1,my
            ycell = ylower + (j-0.5d0)*dy 

            radius = dsqrt((xcell-x0)**2+(ycell-y0)**2)
            q(1,i,j) = q(1,i,j)+1.d0*exp(-100.d0*radius**2)
            q(2:meqn,i,j) = 0.d0

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Adding random noise to the velocity components
            ! call random_number(random_noise1)
            ! call random_number(random_noise2)
            ! q(2,i,j) = q(1,i,j)*(random_noise1-random_noise2)!*1.e-45

            ! call random_number(random_noise1)
            ! call random_number(random_noise2)
            ! q(3,i,j) = q(1,i,j)*(random_noise1-random_noise2)!*1.e-45
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            

            ! radius = sqrt((xcell-x0)**2+(ycell-y0)**2)*scale

            ! call solitary_wave_Bristeau(radius, Hcap, A, d,&
            !     gamma, grav, 0.d0, h, uhat, p, w, inner_radius, direction)  
            
            ! u = uhat*(xcell-x0)*scale/radius
            ! v = uhat*(ycell-y0)*scale/radius

            ! q(1,i,j) = q(1,i,j)+h-Hcap
            ! q(2,i,j) = h*u 
            ! q(3,i,j) = h*v
            ! q(4:meqn,i,j) = 0.d0
            ! if (boussEquations==-1) then
            !     q(4,i,j) = q(1,i,j)*w
            !     q(5,i,j) = q(1,i,j)*p
            ! else if (boussEquations==-2) then
            !     q(4,i,j) = q(1,i,j)*w
            !     q(5,i,j) = 2*q(1,i,j)*w
            !     q(6,i,j) = q(1,i,j)*p
            !     q(7,i,j) = gamma*q(6,i,j)
            ! end if
        enddo
    enddo
    
end subroutine qinit
