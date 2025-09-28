real(kind=8) pure function get_max_speed(val,mitot,mjtot,nvar,aux,naux,nghost,hx,hy)

    use geoclaw_module, only: dry_tolerance, coordinate_system
    use geoclaw_module, only: grav, earth_radius, DEG2RAD

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
    
    ! Arguments
    integer, intent(in) :: mitot,mjtot,nvar,naux,nghost
    real(kind=8), intent(in) :: hx,hy
    real(kind=8), intent(in) :: val(nvar,mitot,mjtot), aux(naux,mitot,mjtot)
    
    ! Locals
    integer :: i,j
    real(kind=8) :: ymetric,hyphys,xmetric,hxphys,u,v,sig,sp_over_h

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local EDC Hyperbolic Relaxation
    double precision :: bath !Local bathymetry
    double precision :: csq_space !Local hyperbolic relaxation parameter
    double precision :: p, w
    ! Local BBBD Hyperbolic Relaxation
    double precision :: sg, pb
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    sp_over_h = 0.d0   ! compute max speed over h, since dx may not equal dy
    if (coordinate_system == 2) then
        do j = nghost+1, mjtot-nghost
            ymetric = earth_radius*deg2rad
            hyphys = ymetric*hy

            do i = nghost+1, mitot-nghost
                xmetric = cos(aux(3,i,j)) * earth_radius * DEG2RAD
                hxphys = xmetric * hx
                if (val(1,i,j) > dry_tolerance) then
                    u  = val(2,i,j) / val(1,i,j)
                    v  = val(3,i,j) / val(1,i,j)
                    if (boussEquations == 0) then
                        sig = sqrt(grav*val(1,i,j))
                    else if (boussEquations == -1) then
                        p = abs(val(5,i,j) / val(1,i,j))
                        csq_space = aux(csq_index,i,j) 
                        sig = sqrt(grav*val(1,i,j)+p+csq_space)
                    else if (boussEquations == -2) then
                        p = abs(val(6,i,j) / val(1,i,j))
                        csq_space = aux(csq_index,i,j) 
                        sig = sqrt(grav*val(1,i,j)+p+csq_space)
                    end if   
                endif
                sp_over_h = max((abs(u)+sig)/hxphys,(abs(v)+sig)/hyphys,sp_over_h)
            end do
        end do
    else  ! speeds in cartesian coords, no metrics needed
        do j = nghost+1, mjtot-nghost
            do i = nghost+1, mitot-nghost
                if (val(1,i,j) > dry_tolerance) then
                    u  = val(2,i,j) / val(1,i,j)
                    v  = val(3,i,j) / val(1,i,j)
                    if (boussEquations == 0) then
                        sig = sqrt(grav*val(1,i,j))
                    else if (boussEquations == -1) then
                        p = abs(val(5,i,j) / val(1,i,j))
                        csq_space = aux(csq_index,i,j) 
                        sig = sqrt(grav*val(1,i,j)+p+csq_space)
                    else if (boussEquations == -2) then
                        p = abs(val(6,i,j) / val(1,i,j))
                        csq_space = aux(csq_index,i,j) 
                        sig = sqrt(grav*val(1,i,j)+p+csq_space)
                    end if   
                    sp_over_h = max((abs(u)+sig)/hx,(abs(v)+sig)/hy,sp_over_h)
                endif
            end do
        end do
    endif
      
    get_max_speed = sp_over_h

end function get_max_speed
