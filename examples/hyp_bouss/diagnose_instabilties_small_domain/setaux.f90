! Set auxiliary arrays
!
! Routine that sets the values found in the aux array including topography and
! other fields that are used in the source terms and other things.
!
! In the default routine sets the following fields depending on the input
! parameters present:
!  (1) topography
!  (2) capacity (set if coordinate_system == 2)
!  (3) length ratio of edge (set if coordinate_system == 2)
!  (frition_index). Location of manning's N (variable_friction)
!  (wind_index, wind_index + 1).  Location of x and y wind speeds (wind_forcing)
!  (pressure_index).  Location of pressure field (pressure_forcing)
!
! There are a couple of additional things of note:
!
!  - This routine is called anytime a grid is created during re-gridding to fill
!    in the aux arrays.
!  - Built-in to this routine is an ability to copy aux values from grids from
!    the same level.  If you do not want this behavior you will need to modify
!    this routine.
!  - If you are using periodic BCs then you must handle this here.  Look to the
!    setting of topography for an example of how this may need to be done.
!  - If your aux arrays are time dependent and set somewhere else it may be
!    prudent to set a default value here as is done with the wind and pressure.
!
subroutine setaux(mbc,mx,my,xlow,ylow,dx,dy,maux,aux)

    use amr_module, only: mcapa, xupper, yupper, xlower, ylower, NEEDS_TO_BE_SET
    use amr_module, only: xperdom, yperdom

    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
    use geoclaw_module, only: sea_level, ambient_pressure

    use storm_module, only: wind_forcing, pressure_forcing
    use storm_module, only: wind_index, pressure_index, set_storm_fields

    use friction_module, only: variable_friction, friction_index
    use friction_module, only: set_friction_field

    use topo_module

    use adjoint_module, only : adjoint_flagging,innerprod_index

    use bouss_module, only : c_sq => boussEDCcsq !Reference hyperbolic relaxation parameter
    use bouss_module, only : boussEquations
    use bouss_module, only: dr_index => boussDecayRate_index
    use bouss_module, only: csq_index => bouss_csq_index
    use bouss_module, only: boussTransitionType !Using initial depth (0) or distance(1)
    use bouss_module, only: transition_type_fun => bouss_transition_type_fun
    use bouss_module, only: transition_low => bouss_transition_low !Min depth or distance
    use bouss_module, only: transition_up => bouss_transition_up !Max depth or distance
    use bouss_module, only: use_bouss !To decide whether to solve SWEs or BTEs
    
    
    use bouss_module, only: transition_function
    use bouss_regions_module, only: distance_to_bouss_regions
    use bouss_regions_module, only: num_bouss_rregions
    
   

    implicit none

    ! Arguments
    integer, intent(in) :: mbc,mx,my,maux
    real(kind=8), intent(in) :: xlow,ylow,dx,dy
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: ii,jj,m, iint,jint
    real(kind=8) :: x,y,xm,ym,xp,yp,topo_integral
    real(kind=8) :: xper,yper,xperm,yperm,xperp,yperp
    character(len=*), parameter :: aux_format = "(2i4,4d15.3)"
    integer :: skipcount,iaux,ilo,jlo
    logical :: exact_distance = .true.

    ! To store temporal data
    real(kind=8) :: h_at_rest, dist_to_regions
    integer :: neighborhood_idx(8,2),iidx,jidx, kk

    ! Lat-Long coordinate system in use, check input variables
    if (coordinate_system == 2) then
        if (mcapa /= 2 .or. maux < 3) then
            print *,'ERROR in setaux:  for coordinate_system==2'
            print *,'     need mcapa == 2 and maux >= 3'
            print *,'     have mcapa = ',mcapa,'  maux = ',maux
            stop
        endif
    endif

    ! Compute integer indices based off same corner of domain to reduce round
    ! off discrepancies
    ilo = floor((xlow - xlower + .05d0*dx)/dx)
    jlo = floor((ylow - ylower + .05d0*dy)/dy)

    ! Set geometry values
    if (coordinate_system == 1) then
        if (maux == 0 .and. mcapa > 0) then
            print *, "ERROR:  Capacity array requested but number of aux"
            print *, "variables is set to 0."
            stop
        end if
    else if (coordinate_system == 2) then
        do jj = 1 - mbc, my + mbc
            do ii = 1 - mbc, mx + mbc
        ! do jj = 1 , my
        !     do ii = 1 , mx
                ym = ylower + (jlo+jj-1.d0) * dy
                yp = ylower + (jlo+jj) * dy
                aux(2,ii,jj) = deg2rad * earth_radius**2                      &
                                * (sin(yp * deg2rad) - sin(ym * deg2rad)) / dy
                aux(3,ii,jj) = ym * deg2rad
            end do
        end do
    end if

    ! Storm fields if used
    if (wind_forcing) then
        aux(wind_index, :, :) = 0.d0
        aux(wind_index + 1, :, :) = 0.d0
    endif
    if (pressure_forcing) then
        aux(pressure_index, :, :) = ambient_pressure
    endif

    ! Innerproduct field if used
    if (adjoint_flagging) then
        do jj=1-mbc,my+mbc
            do ii=1-mbc,mx+mbc
                if (aux(1,ii,jj) .eq. NEEDS_TO_BE_SET) then
                    aux(innerprod_index,ii,jj) = 0.d0
                endif
            enddo
        enddo
    endif


    ! ==========================================================================
    !  Set Bathymetry
    ! Set analytical bathymetry here if requested
    if (test_topography > 0) then
        do jj = 1 - mbc, my + mbc
            ym = ylower + (jlo+jj-1.d0) * dy
            yp = ylower + (jlo+jj) * dy
            y = 0.5d0*(ym+yp)
            do ii = 1 - mbc, mx + mbc
                xm = xlower + (ilo+ii-1.d0) * dx
                xp = xlower + (ilo+ii) * dx
                x = 0.5d0*(xm+xp)
                aux(1,ii,jj) = test_topo(x)
            end do
        end do

    ! Set bathymetry based on input files
    else if (mtopofiles > 0) then

        do jj=1-mbc,my+mbc
            ym = ylower + (jlo+jj-1.d0) * dy
            yp = ylower + (jlo+jj) * dy
            y = 0.5d0*(ym+yp)


            do ii=1-mbc,mx+mbc
                xm = xlower + (ilo+ii-1.d0) * dx
                xp = xlower + (ilo+ii) * dx
                x = 0.5d0*(xm+xp)

                ! Parameter NEEDS_TO_BE_SET initialized in amr_module.f90
                ! saves time by otherwise copying instead of reinitializing
                if (aux(1,ii,jj) .ne. NEEDS_TO_BE_SET) then
                    cycle
                endif

                topo_integral = 0.d0
                if ((y>yupper).or.(y<ylower).or.(x>xupper).or.(x<xlower)) then
                    if (.not.(xperdom .or. yperdom)) then
                        ! Skip setting as this cell sticks out of the physical
                        ! domain and we are not setting periodic BCs
                        cycle

                    else
                        ! We evaluate the periodic BC topography by computing
                        ! the appropriate periodic grid cell coordinates and
                        ! again calling cellgridintegrate
                        call wrap_coords(x,y,xperm,xper,xperp,yperm,yper, &
                                         yperp,dx,dy)
                        call cellgridintegrate(topo_integral,  &
                                               xperm,xperp,yperm,yperp)
                    endif
                else
                    ! Cell does not extend outside of physical domain
                    call cellgridintegrate(topo_integral,xm,xp,ym,yp)
                endif

                ! Correct for geometry
                if (coordinate_system == 2) then
                    aux(1,ii,jj) = topo_integral / (dx * dy * aux(2,ii,jj))
                else
                    aux(1,ii,jj) = topo_integral / (dx * dy)
                endif
            enddo
        enddo
    else
        print *, "ERROR:  There is no way to set bathymetry!  Either "
        print *, "        provide topography files or request topography "
        print*,  "        defined by a function."
        stop
    end if

    ! Copy topo to ghost cells if outside physical domain and not periodic
    if (.not. yperdom) then
        do jj=1-mbc,my+mbc
            y = ylower + (jlo+jj-.5d0) * dy
            if ((y < ylower) .or. (y>yupper)) then
                do ii=1-mbc,mx+mbc
                    x = xlower + (ilo+ii-.5d0) * dx
                    iint = ii + max(0, ceiling((xlower-x)/dx)) &
                         - max(0, ceiling((x-xupper)/dx))
                    jint = jj + max(0, ceiling((ylower-y)/dy)) &
                         - max(0, ceiling((y-yupper)/dy))
                    aux(1,ii,jj) = aux(1,iint,jint)
                    !Copy also BTEs aux variables
                    if (boussEquations<0) then
                        aux(csq_index,ii,jj) = aux(csq_index,iint,jint)
                        aux(dr_index,ii,jj) = aux(dr_index,iint,jint)
                    endif
                enddo
            endif
        enddo
    endif
    if (.not. xperdom) then
        do ii=1-mbc,mx+mbc
            x =  xlower + (ilo+ii-.5d0) * dx
            if ((x < xlower) .or. (x > xupper)) then
                do jj=1-mbc,my+mbc
                    y = ylower + (jlo+jj-.5d0) * dy
                    iint = ii + max(0, ceiling((xlower-x)/dx)) &
                         - max(0, ceiling((x-xupper)/dx))
                    jint = jj + max(0, ceiling((ylower-y)/dy)) &
                         - max(0, ceiling((y-yupper)/dy))
                    aux(1,ii,jj) = aux(1,iint,jint)
                    !Copy also BTEs aux variables
                    if (boussEquations<0) then
                        aux(csq_index,ii,jj) = aux(csq_index,iint,jint)
                        aux(dr_index,ii,jj) = aux(dr_index,iint,jint)
                    endif
                enddo
            endif
        enddo
    endif



    ! If using a variable friction field initialize the coefficients to 0
    if (variable_friction) then
        call set_friction_field(mx, my, mbc, maux, xlow, ylow, dx, dy, aux)
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Set Hyperbolic relaxation parameter for EDC and BBBD systems
    if (boussEquations<0) then
        ! if (.not. use_bouss) then
        !     !Set aux to zero
            ! aux(dr_index,:,:) = 0.d0
            ! aux(csq_index,:,:) = 0.d0
        !     ! return !Dangerous if something else in setaux is modified below
        ! endif
        !Case statement
        select case(boussTransitionType)
        case(0) !Depth based transition
            !idxs for buffer zone
            neighborhood_idx(:,1) = [-1,-1,-1,0 ,0, 1,1,1]
            neighborhood_idx(:,2) = [-1, 0, 1,-1,1,-1,0,1]
            do jj=1,my
                do ii=1,mx
                    ! if (aux(dr_index,ii,jj) .ne. NEEDS_TO_BE_SET) then
                    !     !Do not overwrite values that have been set
                    !     cycle
                    ! endif
                    h_at_rest = sea_level - aux(1,ii,jj)
                    !Go around neighborhood of this cell
                    ! if(jj>2-mbc .and. jj<my+mbc-1 .and. ii>2-mbc .and. ii<mx+mbc-1) then
                        do kk=1,8
                            iidx = ii + neighborhood_idx(kk,1)
                            jidx = jj + neighborhood_idx(kk,2)
                            h_at_rest = min(h_at_rest,sea_level - aux(1,iidx,jidx))
                        end do
                        aux(dr_index,ii,jj) = transition_function(h_at_rest,transition_low,&
                                    transition_up,transition_type_fun)
                        aux(csq_index,ii,jj) = c_sq * aux(dr_index,ii,jj)
                    ! end if
                enddo
            enddo
        case(1) !Distance based transition
            if (0 == num_bouss_rregions) then
                print *, "ERROR:  No regions specified for Boussinesq equations"
                print *, "Distance based transition requires regions to be set"
                call abort
            endif
            do jj=1-mbc,my+mbc
                ym = ylower + (jlo+jj-1.d0) * dy
                yp = ylower + (jlo+jj) * dy
                y = 0.5d0*(ym+yp)
                do ii=1-mbc,mx+mbc
                    ! if (aux(dr_index,ii,jj) .ne. NEEDS_TO_BE_SET) then
                    !     cycle
                    ! endif
                    xm = xlower + (ilo+ii-1.d0) * dx
                    xp = xlower + (ilo+ii) * dx
                    x = 0.5d0*(xm+xp)
                    !To get the exact distance instead of an approximate one
                    !Set exact_distance = .true. 
                    ! exact_distance = .true.
                    dist_to_regions = distance_to_bouss_regions(x,y,exact_distance)
                    aux(dr_index,ii,jj) = transition_function(dist_to_regions,transition_low,&
                                    transition_up,transition_type_fun)
                    aux(csq_index,ii,jj) = c_sq * aux(dr_index,ii,jj)
                end do
            end do
        case default
            print *, "ERROR:  Unknown transition type for Boussinesq equations"
            call abort
        end select
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

    ! Provide wrapper function for providing periodic coordinates
    subroutine wrap_coords(x,y,xperm,xper,xperp,yperm,yper,yperp,dx,dy)

        use amr_module, only: xperdom, yperdom, xupper, yupper, xlower, ylower

        implicit none

        ! Arguments
        real(kind=8), intent(in) :: x,y,dx,dy
        real(kind=8), intent(out) :: xperm,xper,xperp,yperm,yper,yperp
        real(kind=8) :: xdomain, ydomain

        ! need size of domain for wrapping
        xdomain = xupper - xlower
        ydomain = yupper - ylower
        xper = x
        yper = y

        ! test for wrapping of coordinates
        if (x > xupper .and. xperdom) xper = xper - xdomain
        if (x < xlower .and. xperdom) xper = xper + xdomain
        if (y > yupper .and. yperdom) yper = yper - ydomain
        if (y < ylower .and. yperdom) yper = yper + ydomain

        ! adjust rest of variables
        xperm = xper - 0.5d0 * dx
        xperp = xper + 0.5d0 * dx
        yperm = yper - 0.5d0 * dy
        yperp = yper + 0.5d0 * dy

    end subroutine wrap_coords

end subroutine setaux

