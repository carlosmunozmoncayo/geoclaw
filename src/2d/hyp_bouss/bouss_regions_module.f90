module bouss_regions_module

    use regions_module, only: region_type, ruled_region_type
    use amr_module, only: parmunit
    !!!!!!!
    !Necessary stuff for coordinate_system==2:
    use geoclaw_module, only: DEG2RAD, earth_radius, spherical_distance
    use geoclaw_module, only: coordinate_system,latlon2xy, xy2latlon
    ! use bouss_module, only: projection_center
    !!!!!!!
    implicit none

    type point
        real(kind=8) :: x,y
    end type point

    integer :: num_bouss_rregions
    type(ruled_region_type), target, allocatable :: bouss_rregions(:)
    contains

    subroutine set_bouss_rregions(fname)

        ! Read new-style Ruled Rectangles from flagregions.data
        implicit none
      
        ! Function Arguments
        character(len=*), optional, intent(in) :: fname
      
        ! Locals
        integer, parameter :: unit = 7
        integer, parameter :: unit2 = 45
        integer :: i,j,nrules1, spatial_region_type
        logical :: foundFile
        type(ruled_region_type), pointer :: rr
        type(ruled_region_type), pointer :: rr_xy
        real(kind=8) :: rr_x1,rr_x2,rr_y1,rr_y2
        
        write(parmunit,*) ' '
        write(parmunit,*) '--------------------------------------------'
        write(parmunit,*) 'BOUSS RULED REGIONS:'
        write(parmunit,*) '-----------'

        if (present(fname)) then
            call opendatafile(unit,fname)
        else
            call opendatafile(unit,'bouss_flagregions.data')
        endif

        read(unit,"(i2)") num_bouss_rregions
        if (num_bouss_rregions == 0) then
            write(parmunit,*) '  No ruled regions specified for Boussinesq-type system'
            
        else
            ! BTEs region data
            allocate(bouss_rregions(num_bouss_rregions))
            ! Read ruled regions
            do i=1,num_bouss_rregions
                rr => bouss_rregions(i)
                read(unit,*) rr%name
                rr%name = trim(rr%name)
                read(unit,*) rr%min_level
                read(unit,*) rr%max_level
                read(unit,*) rr%t_low
                read(unit,*) rr%t_hi
                read(unit,*) spatial_region_type

                if (spatial_region_type == 1) then
                    ! read rectangle extent:
                    read(unit,*) rr_x1,rr_x2,rr_y1,rr_y2
                    ! turn into ruled rectangle:
                    rr%ixy = 1
                    rr%method = 0
                    rr%nrules = 2
                    allocate(rr%s(rr%nrules), rr%lower(rr%nrules), & 
                             rr%upper(rr%nrules))
                    rr%s(1) = rr_x1
                    rr%s(2) = rr_x2
                    rr%ds = rr_x2 - rr_x1
                    rr%lower(1) = rr_y1
                    rr%upper(1) = rr_y2
                    rr%lower(2) = rr_y1
                    rr%upper(2) = rr_y2
                
                else if (spatial_region_type == 2) then    
                    read(unit,*) rr%file_name
                    write(6,*) '+++ Ruled region name: ', rr%name
                    write(6,*) '+++ Ruled region file_name: ', rr%file_name
                    inquire(file=trim(rr%file_name),exist=foundFile)
                    if (.not. foundFile) then
                      write(*,*) 'Missing rregions file...'
                      write(*,*) 'Looking for: ',trim(rr%file_name)
                      stop
                      endif

                    open(unit=unit2,file=trim(rr%file_name),status='old')
                    read(unit2,*) rr%ixy
                    read(unit2,*) rr%method
                    read(unit2,*) rr%ds
                    read(unit2,*) rr%nrules
                    allocate(rr%s(rr%nrules), rr%lower(rr%nrules), & 
                             rr%upper(rr%nrules))
                    do j=1,rr%nrules
                        read(unit2,*) rr%s(j), rr%lower(j), rr%upper(j)
                    enddo
                        
                else
                    write(6,*) '*** Error: unexpected spatial_region_type'
                endif
                    
                ! compute bounding box:
                
                if (rr%method == 0) then
                    nrules1 = rr%nrules - 1
                else
                    nrules1 = rr%nrules
                endif
                    
                if (rr%ixy == 1) then
                    rr%x1bb = rr%s(1)
                    rr%x2bb = rr%s(rr%nrules)
                    rr%y1bb = minval(rr%lower(1:nrules1))
                    rr%y2bb = maxval(rr%upper(1:nrules1))
                else
                    rr%y1bb = rr%s(1)
                    rr%y2bb = rr%s(rr%nrules)
                    rr%x1bb = minval(rr%lower(1:nrules1))
                    rr%x2bb = maxval(rr%upper(1:nrules1))
                endif
                    
                write(6,*) '+++ rregion bounding box: '
                write(6,*) rr%x1bb,rr%x2bb,rr%y1bb,rr%y2bb
                
                write(6,*) '+++ i, rr%s(1), rr%ds: ',i, rr%s(1), rr%ds
            enddo
            nullify(rr)
        endif
        close(unit)


    end subroutine set_bouss_rregions

    function distance_to_bouss_regions(x,y,exact) result(dist)
        !Computes distance from point (x,y) to the nearest ruled rectangle
        !in the list of ruled regions bouss_rregions
        implicit none
        real(kind=8), intent(in) :: x,y
        logical, optional, intent(in) :: exact
        logical :: exact_flag
        real(kind=8) :: dist
        integer :: i
        dist = 1.e20
        if (present(exact)) then
            exact_flag = exact
        else
            exact_flag = .false. !default
        end if
        if (exact_flag) then
            do i=1,num_bouss_rregions
                dist = min(dist,distance_to_ruled_rectangle_exact(point(x,y),&
                bouss_rregions(i)))
            end do
        else
            do i=1,num_bouss_rregions
                if (bouss_rregions(i)%method == 1) then
                    dist = min(dist,distance_to_ruled_rectangle_method1(point(x,y),&
                    bouss_rregions(i)))
                else !bouss_rregions(i)%method == 0
                    dist = min(dist,distance_to_ruled_rectangle_method0(point(x,y),&
                    bouss_rregions(i)))
                end if
            end do
        end if
    end function distance_to_bouss_regions

    function distance_to_ruled_rectangle_method1(p,rr) result(dist)
        !Computes the distance from the point p to the ruled rectangle rr
        !Using an heuristic method. Not exact, but fast.
        implicit none
        type(point), intent(in) :: p
        type(ruled_region_type), intent(in) :: rr
        real(kind=8) :: dist

        !Locals (I'll define them as I need them)
        integer :: idx, nr
        real(kind=8) :: xlow,xup,ylow,yup
        type(point) :: point1,point2,point3,point4
        real(kind=8) :: d1,d2,d3

        nr = rr%nrules


        if (rr%ixy==1) then !s is an array in x
            idx = binary_search(p%x,rr%s) !s(idx-1) < p%x <= s(idx)
            !Edge cases
            if (1>=idx) then
                point1 = point(rr%s(1),rr%lower(1))
                point2 = point(rr%s(1),rr%upper(1))
                point3 = point(rr%s(2),rr%lower(2))
                point4 = point(rr%s(2),rr%upper(2))
                d1 = distance_to_segment(p,point1,point2)
                d2 = distance_to_segment(p,point1,point3)
                d3 = distance_to_segment(p,point2,point4)
                dist = min(d1,d2,d3)
                return
            end if
            if (nr<idx) then
                point1 = point(rr%s(nr),rr%lower(nr))
                point2 = point(rr%s(nr),rr%upper(nr))
                point3 = point(rr%s(nr-1),rr%lower(nr-1))
                point4 = point(rr%s(nr-1),rr%upper(nr-1))
                d1 = distance_to_segment(p,point1,point2)
                d2 = distance_to_segment(p,point1,point3)
                d3 = distance_to_segment(p,point2,point4)
                dist = min(d1,d2,d3)
                return
            end if
            !Interior cases
            ylow = get_y_line(p%x,point(rr%s(idx-1),rr%lower(idx-1)), &
                            point(rr%s(idx),rr%lower(idx)))
            yup = get_y_line(p%x,point(rr%s(idx-1),rr%upper(idx-1)), &
                            point(rr%s(idx),rr%upper(idx)))
            !Check if point is inside ruled rectangle
            if (ylow <= p%y .and. p%y <= yup) then
                dist = 0.
                return
            else if (ylow > p%y) then
                !Point is below the ruled rectangle
                if (idx == 2) then
                    !x between s(1) and s(2)
                    !Just compute distance to 2 segments
                    point1 = point(rr%s(1),rr%lower(1))
                    point2 = point(rr%s(2),rr%lower(2))
                    point3 = point(rr%s(min(3,nr)),rr%lower(min(3,nr)))
                    d1 = distance_to_segment(p,point1,point2)
                    d2 = distance_to_segment(p,point2,point3)
                    dist = min(d1,d2)
                    return
                else if (idx == nr) then
                    !x between s(nr-1) and s(nr)
                    !Just compute distance to 2 segments
                    point1 = point(rr%s(nr),rr%lower(nr))
                    point2 = point(rr%s(nr-1),rr%lower(nr-1))
                    point3 = point(rr%s(max(1,nr-2)),rr%lower(max(1,nr-2)))
                    d1 = distance_to_segment(p,point1,point2)
                    d2 = distance_to_segment(p,point2,point3)
                    dist = min(d1,d2)
                    return
                else
                    !x between s(idx-1) and s(idx)
                    !Compute distance to 3 segments connecting
                    !the points (s(i),lower(i)) for i=idx-2,idx+1
                    point1 = point(rr%s(idx-2),rr%lower(idx-2))
                    point2 = point(rr%s(idx-1),rr%lower(idx-1))
                    point3 = point(rr%s(idx),rr%lower(idx))
                    point4 = point(rr%s(idx+1),rr%lower(idx+1))
                    d1 = distance_to_segment(p,point1,point2)
                    d2 = distance_to_segment(p,point2,point3)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                end if
            else ! yup < p%y
                !Point is above the ruled rectangle
                if (idx == 2) then
                    !x between s(1) and s(2)
                    !Just compute distance to 2 segments
                    point1 = point(rr%s(1),rr%upper(1))
                    point2 = point(rr%s(2),rr%upper(2))
                    point3 = point(rr%s(min(3,nr)),rr%upper(min(3,nr)))
                    d1 = distance_to_segment(p,point1,point2)
                    d2 = distance_to_segment(p,point2,point3)
                    dist = min(d1,d2)
                    return
                else if (idx == nr) then
                    !x between s(nr-1) and s(nr)
                    !Just compute distance to 2 segments
                    point1 = point(rr%s(nr),rr%upper(nr))
                    point2 = point(rr%s(nr-1),rr%upper(nr-1))
                    point3 = point(rr%s(max(1,nr-2)),rr%upper(max(1,nr-2)))
                    d1 = distance_to_segment(p,point1,point2)
                    d2 = distance_to_segment(p,point2,point3)
                    dist = min(d1,d2)
                    return
                else
                    !x between s(idx-1) and s(idx)
                    !Compute distance to 3 segments connecting
                    !the points (s(i),upper(i)) for i=idx-2,idx+1
                    point1 = point(rr%s(idx-2),rr%upper(idx-2))
                    point2 = point(rr%s(idx-1),rr%upper(idx-1))
                    point3 = point(rr%s(idx),rr%upper(idx))
                    point4 = point(rr%s(idx+1),rr%upper(idx+1))
                    d1 = distance_to_segment(p,point1,point2)
                    d2 = distance_to_segment(p,point2,point3)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                end if
            end if

        else  !(rr%ixy==2) s is an array in y
            idx = binary_search(p%y,rr%s) !s(idx-1) < p%y <= s(idx)
            !Edge cases
            if (1>=idx) then
                point1 = point(rr%lower(1),rr%s(1))
                point2 = point(rr%upper(1),rr%s(1))
                point3 = point(rr%lower(2),rr%s(2))
                point4 = point(rr%upper(2),rr%s(2))
                d1 = distance_to_segment(p,point1,point2)
                d2 = distance_to_segment(p,point1,point3)
                d3 = distance_to_segment(p,point2,point4)
                dist = min(d1,d2,d3)
                return
            end if
            if (nr<idx) then
                point1 = point(rr%lower(nr),rr%s(nr))
                point2 = point(rr%upper(nr),rr%s(nr))
                point3 = point(rr%lower(nr-1),rr%s(nr-1))
                point4 = point(rr%upper(nr-1),rr%s(nr-1))
                d1 = distance_to_segment(p,point1,point2)
                d2 = distance_to_segment(p,point1,point3)
                d3 = distance_to_segment(p,point2,point4)
                dist = min(d1,d2,d3)
                return
            end if
            !Interior cases
            xlow = get_x_line(p%y,point(rr%lower(idx-1),rr%s(idx-1)), &
                            point(rr%lower(idx),rr%s(idx)))
            xup = get_x_line(p%y,point(rr%upper(idx-1),rr%s(idx-1)), &
                            point(rr%upper(idx),rr%s(idx)))
            !Check if point is inside ruled rectangle
            if (xlow <= p%x .and. p%x <= xup) then
                dist = 0.
                return
            else if (xlow > p%x) then
                !Point is to the left of the ruled rectangle
                if (idx == 2) then
                    !y between s(1) and s(2)
                    !Just compute distance to 2 segments
                    point1 = point(rr%lower(1),rr%s(1))
                    point2 = point(rr%lower(2),rr%s(2))
                    point3 = point(rr%lower(3),rr%s(3))
                    d1 = distance_to_segment(p,point1,point2)
                    d2 = distance_to_segment(p,point2,point3)
                    dist = min(d1,d2)
                    return
                else if (idx == nr) then
                    !y between s(nr-1) and s(nr)
                    !Just compute distance to 2 segments
                    point1 = point(rr%lower(nr),rr%s(nr))
                    point2 = point(rr%lower(nr-1),rr%s(nr-1))
                    point3 = point(rr%lower(nr-2),rr%s(nr-2))
                    d1 = distance_to_segment(p,point1,point2)
                    d2 = distance_to_segment(p,point2,point3)
                    dist = min(d1,d2)
                    return
                else
                    !y between s(idx-1) and s(idx)
                    !Compute distance to 3 segments connecting
                    !the points (lower(i),s(i)) for i=idx-2,idx+1
                    point1 = point(rr%lower(idx-2),rr%s(idx-2))
                    point2 = point(rr%lower(idx-1),rr%s(idx-1))
                    point3 = point(rr%lower(idx),rr%s(idx))
                    point4 = point(rr%lower(idx+1),rr%s(idx+1))
                    d1 = distance_to_segment(p,point1,point2)
                    d2 = distance_to_segment(p,point2,point3)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                end if
            else ! xup < p%x
                !Point is to the right of the ruled rectangle
                if (idx == 2) then
                    !y between s(1) and s(2)
                    !Just compute distance to 2 segments
                    point1 = point(rr%upper(1),rr%s(1))
                    point2 = point(rr%upper(2),rr%s(2))
                    point3 = point(rr%upper(3),rr%s(3))
                    d1 = distance_to_segment(p,point1,point2)
                    d2 = distance_to_segment(p,point2,point3)
                    dist = min(d1,d2)
                    return
                else if (idx == nr) then
                    !y between s(nr-1) and s(nr)
                    !Just compute distance to 2 segments
                    point1 = point(rr%upper(nr),rr%s(nr))
                    point2 = point(rr%upper(nr-1),rr%s(nr-1))
                    point3 = point(rr%upper(nr-2),rr%s(nr-2))
                    d1 = distance_to_segment(p,point1,point2)
                    d2 = distance_to_segment(p,point2,point3)
                    dist = min(d1,d2)
                    return
                else
                    !y between s(idx-1) and s(idx)
                    !Compute distance to 3 segments connecting
                    !the points (upper(i),s(i)) for i=idx-2,idx+1
                    point1 = point(rr%upper(idx-2),rr%s(idx-2))
                    point2 = point(rr%upper(idx-1),rr%s(idx-1))
                    point3 = point(rr%upper(idx),rr%s(idx))
                    point4 = point(rr%upper(idx+1),rr%s(idx+1))
                    d1 = distance_to_segment(p,point1,point2)
                    d2 = distance_to_segment(p,point2,point3)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                end if
            end if
        end if
    end function distance_to_ruled_rectangle_method1
    
    function distance_to_ruled_rectangle_method0(p,rr) result(dist)
        !Computes the distance from the point p to the ruled rectangle rr
        !Using an heuristic method. Not exact, but fast.
        implicit none
        type(point), intent(in) :: p
        type(ruled_region_type), intent(in) :: rr
        real(kind=8) :: dist

        !Locals (I'll define them as I need them)
        integer :: idx, nr
        real(kind=8) :: xlow,xup,ylow,yup
        type(point) :: point1,point2,point3,point4
        type(point) :: point5,point6
        real(kind=8) :: d1,d2,d3

        nr = rr%nrules

        if (1==rr%ixy) then
            idx = binary_search(p%x,rr%s) !s(idx-1) < p%x <= s(idx)
            !Edge cases
            if (1>=idx) then
                !Line 1 (left edge)
                point1 = point(rr%s(1),rr%lower(1))
                point2 = point(rr%s(1),rr%upper(1))
                !Line 2 (below)
                point3 = point(rr%s(2),rr%lower(1))
                point4 = point(rr%s(2),rr%lower(2))
                !Line 3 (above)
                point5 = point(rr%s(2),rr%upper(1))
                point6 = point(rr%s(2),rr%upper(2))
                d1 = distance_to_segment(p,point1,point2)
                d2 = distance_to_segment(p,point3,point4)
                d3 = distance_to_segment(p,point5,point6)
                dist = min(d1,d2,d3)
                return
            else if (nr<idx) then
                !Line 1 (right edge)
                point1 = point(rr%s(nr),rr%lower(nr-1))
                point2 = point(rr%s(nr),rr%upper(nr-1))
                !Line 2 (below)
                point3 = point(rr%s(nr-1),rr%lower(nr-1))
                point4 = point(rr%s(nr-1),rr%lower(nr-2))
                !Line 3 (above)
                point5 = point(rr%s(nr-1),rr%upper(nr-1))
                point6 = point(rr%s(nr-1),rr%upper(nr-2))
                d1 = distance_to_segment(p,point1,point2)
                d2 = distance_to_segment(p,point3,point4)
                d3 = distance_to_segment(p,point5,point6)
                dist = min(d1,d2,d3)
                return
            end if
            !Interior cases
            ylow = rr%lower(idx-1)
            yup = rr%upper(idx-1)
            !Check if point is inside ruled rectangle
            if (ylow<= p%y .and. p%y<= yup) then
                dist = 0.d0
                return
            else if (ylow > p%y) then
                d1 = ylow - p%y
                if (2==idx) then
                    !2 point buffer zone
                    point1 = point(rr%s(2),rr%lower(1))
                    point2 = point(rr%s(2),rr%lower(2))
                    point3 = point(rr%s(3),rr%lower(2))
                    point4 = point(rr%s(3),rr%lower(3))
                    d2 = distance_to_segment(p,point1,point2)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                else if (nr==idx) then
                    !2 point buffer zone
                    point1 = point(rr%s(nr-1),rr%lower(nr-1))
                    point2 = point(rr%s(nr-1),rr%lower(nr-2))
                    point3 = point(rr%s(nr-2),rr%lower(nr-2))
                    point4 = point(rr%s(nr-2),rr%lower(nr-3))
                    d2 = distance_to_segment(p,point1,point2)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                else
                    !1 point buffer zone
                    point1 = point(rr%s(idx),rr%lower(idx-1))
                    point2 = point(rr%s(idx),rr%lower(idx))
                    point3 = point(rr%s(idx-1),rr%lower(idx-1))
                    point4 = point(rr%s(idx-1),rr%lower(idx-2))
                    d2 = distance_to_segment(p,point1,point2)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                end if
            else ! yup < p%y
                d1 = p%y - yup
                if (2==idx) then
                    !2 point buffer zone
                    point1 = point(rr%s(2),rr%upper(1))
                    point2 = point(rr%s(2),rr%upper(2))
                    point3 = point(rr%s(3),rr%upper(2))
                    point4 = point(rr%s(3),rr%upper(3))
                    d2 = distance_to_segment(p,point1,point2)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                else if (nr==idx) then
                    !2 point buffer zone
                    point1 = point(rr%s(nr-1),rr%upper(nr-1))
                    point2 = point(rr%s(nr-1),rr%upper(nr-2))
                    point3 = point(rr%s(nr-2),rr%upper(nr-2))
                    point4 = point(rr%s(nr-2),rr%upper(nr-3))
                    d2 = distance_to_segment(p,point1,point2)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                else
                    !1 point buffer zone
                    point1 = point(rr%s(idx),rr%upper(idx-1))
                    point2 = point(rr%s(idx),rr%upper(idx))
                    point3 = point(rr%s(idx-1),rr%upper(idx-1))
                    point4 = point(rr%s(idx-1),rr%upper(idx-2))
                    d2 = distance_to_segment(p,point1,point2)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                end if
            end if
        else !2==rr%ixy
            idx = binary_search(p%y,rr%s) !s(idx-1) < p%y <= s(idx)
            !Edge cases
            if (1>=idx) then
                !Line 1 (left edge)
                point1 = point(rr%lower(1),rr%s(1))
                point2 = point(rr%upper(1),rr%s(1))
                !Line 2 (below)
                point3 = point(rr%lower(1),rr%s(2))
                point4 = point(rr%lower(2),rr%s(2))
                !Line 3 (above)
                point5 = point(rr%upper(1),rr%s(2))
                point6 = point(rr%upper(2),rr%s(2))
                d1 = distance_to_segment(p,point1,point2)
                d2 = distance_to_segment(p,point3,point4)
                d3 = distance_to_segment(p,point5,point6)
                dist = min(d1,d2,d3)
                return
            else if (nr<idx) then
                !Line 1 (right edge)
                point1 = point(rr%lower(nr-1),rr%s(nr))
                point2 = point(rr%upper(nr-1),rr%s(nr))
                !Line 2 (below)
                point3 = point(rr%lower(nr-1),rr%s(nr-1))
                point4 = point(rr%lower(nr-2),rr%s(nr-1))
                !Line 3 (above)
                point5 = point(rr%upper(nr-1),rr%s(nr-1))
                point6 = point(rr%upper(nr-2),rr%s(nr-1))
                d1 = distance_to_segment(p,point1,point2)
                d2 = distance_to_segment(p,point3,point4)
                d3 = distance_to_segment(p,point5,point6)
                dist = min(d1,d2,d3)
                return
            end if
            !Interior cases
            xlow = rr%lower(idx-1)
            xup = rr%upper(idx-1)
            !Check if point is inside ruled rectangle
            if (xlow<= p%x .and. p%x<= xup) then
                dist = 0.d0
                return
            else if (xlow > p%x) then
                d1 = xlow - p%x
                if (2==idx) then
                    !2 point buffer zone
                    point1 = point(rr%lower(1),rr%s(2))
                    point2 = point(rr%lower(2),rr%s(2))
                    point3 = point(rr%lower(2),rr%s(3))
                    point4 = point(rr%lower(3),rr%s(3))
                    d2 = distance_to_segment(p,point1,point2)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                else if (nr==idx) then
                    !2 point buffer zone
                    point1 = point(rr%lower(nr-1),rr%s(nr))
                    point2 = point(rr%lower(nr-1),rr%s(nr-1))
                    point3 = point(rr%lower(nr-2),rr%s(nr-1))
                    point4 = point(rr%lower(nr-2),rr%s(nr-2))
                    d2 = distance_to_segment(p,point1,point2)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                else
                    !1 point buffer zone
                    point1 = point(rr%lower(idx-1),rr%s(idx))
                    point2 = point(rr%lower(idx),rr%s(idx))
                    point3 = point(rr%lower(idx-1),rr%s(idx-1))
                    point4 = point(rr%lower(idx-2),rr%s(idx-1))
                    d2 = distance_to_segment(p,point1,point2)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                end if
            else ! xup < p%x
                d1 = p%x - xup
                if (2==idx) then
                    !2 point buffer zone
                    point1 = point(rr%upper(1),rr%s(2))
                    point2 = point(rr%upper(2),rr%s(2))
                    point3 = point(rr%upper(2),rr%s(3))
                    point4 = point(rr%upper(3),rr%s(3))
                    d2 = distance_to_segment(p,point1,point2)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                else if (nr==idx) then
                    !2 point buffer zone
                    point1 = point(rr%upper(nr-1),rr%s(nr))
                    point2 = point(rr%upper(nr-1),rr%s(nr-1))
                    point3 = point(rr%upper(nr-2),rr%s(nr-1))
                    point4 = point(rr%upper(nr-2),rr%s(nr-2))
                    d2 = distance_to_segment(p,point1,point2)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                else
                    !1 point buffer zone
                    point1 = point(rr%upper(idx-1),rr%s(idx))
                    point2 = point(rr%upper(idx),rr%s(idx))
                    point3 = point(rr%upper(idx-1),rr%s(idx-1))
                    point4 = point(rr%upper(idx-2),rr%s(idx-1))
                    d2 = distance_to_segment(p,point1,point2)
                    d3 = distance_to_segment(p,point3,point4)
                    dist = min(d1,d2,d3)
                    return
                end if
            end if
        end if

    end function distance_to_ruled_rectangle_method0
    
    function distance_to_ruled_rectangle_exact(p,rr) result(dist)
        !Computes the distance from the point p to the ruled rectangle rr
        !Using an exact method. Slow! Can be generalized to any 
        !polygon, but I'm sticking to the ruled rectangle data structure
        !for now.
        implicit none
        type(point), intent(in) :: p
        type(ruled_region_type), intent(in) :: rr
        real(kind=8) :: dist

        !Locals
        integer :: i, idx
        type(point) :: point1,point2
        real(kind=8) :: yup,ylow,xup,xlow

        dist = 1.e20
        if (1==rr%method) then
            if (1==rr%ixy) then !s is horizontal
                !Check if point is inside ruled rectangle
                !This also helps us halve the number of segments to check
                idx = binary_search(p%x,rr%s) 
                if (idx>1 .and. idx<=rr%nrules) then
                    !Could be inside ruled rectangle
                    yup = get_y_line(p%x,point(rr%s(idx-1),rr%upper(idx-1)), &
                                    point(rr%s(idx),rr%upper(idx)))
                    ylow = get_y_line(p%x,point(rr%s(idx-1),rr%lower(idx-1)), &
                                    point(rr%s(idx),rr%lower(idx)))
                    if (ylow <= p%y .and. p%y <= yup) then
                        dist = 0.d0 !Point is inside ruled rectangle
                        return
                    else if (ylow > p%y) then
                        !Iterate over lower segments
                        do i=1,rr%nrules-1
                            point1 = point(rr%s(i),rr%lower(i))
                            point2 = point(rr%s(i+1),rr%lower(i+1))
                            dist = min(dist,distance_to_segment(p,point1,point2))
                        end do
                        return
                    else ! yup < p%y
                        !Iterate over upper segments
                        do i=1,rr%nrules-1
                            point1 = point(rr%s(i),rr%upper(i))
                            point2 = point(rr%s(i+1),rr%upper(i+1))
                            dist = min(dist,distance_to_segment(p,point1,point2))
                        end do
                        return
                    end if
                end if
                !Not inside ruled rectangle, iterate over all segments
                !Start at bottom left corner and move clockwise
                !Go to top left corner
                point1 = point(rr%s(1),rr%lower(1))
                point2 = point(rr%s(1),rr%upper(1))
                dist = min(dist,distance_to_segment(p,point1,point2))
                !Go to top right corner
                do i=2,rr%nrules
                    point1 = point2
                    point2 = point(rr%s(i),rr%upper(i))
                    dist = min(dist,distance_to_segment(p,point1,point2))
                end do
                !Go to bottom right corner
                point1 = point2
                point2 = point(rr%s(rr%nrules),rr%lower(rr%nrules))
                dist = min(dist,distance_to_segment(p,point1,point2))
                !Go to bottom left corner
                do i=rr%nrules-1,2,-1
                    point1 = point2
                    point2 = point(rr%s(i-1),rr%lower(i-1))
                    dist = min(dist,distance_to_segment(p,point1,point2))
                end do
                return
            else !s is vertical
                !Check if point is inside ruled rectangle
                !This also helps us halve the number of segments to check
                idx = binary_search(p%y,rr%s) 
                if (idx>1 .and. idx<=rr%nrules) then
                    !Could be inside ruled rectangle
                    xup = get_x_line(p%y,point(rr%upper(idx-1),rr%s(idx-1)), &
                                    point(rr%upper(idx),rr%s(idx)))
                    xlow = get_x_line(p%y,point(rr%lower(idx-1),rr%s(idx-1)), &
                                    point(rr%lower(idx),rr%s(idx)))
                    if (xlow <= p%x .and. p%x <= xup) then
                        dist = 0.d0 !Point is inside ruled rectangle
                        return
                    else if (xlow > p%x) then
                        !Iterate over left segments
                        do i=1,rr%nrules-1
                            point1 = point(rr%lower(i),rr%s(i))
                            point2 = point(rr%lower(i+1),rr%s(i+1))
                            dist = min(dist,distance_to_segment(p,point1,point2))
                        end do
                        return
                    else ! xup < p%x
                        !Iterate over right segments
                        do i=1,rr%nrules-1
                            point1 = point(rr%upper(i),rr%s(i))
                            point2 = point(rr%upper(i+1),rr%s(i+1))
                            dist = min(dist,distance_to_segment(p,point1,point2))
                        end do
                        return
                    end if
                    
                end if
                !Not inside ruled rectangle, iterate over all segments
                !Start at bottom right corner and move clockwise
                !Go to bottom left corner
                point1 = point(rr%upper(1),rr%s(1))
                point2 = point(rr%lower(1),rr%s(1))
                dist = min(dist,distance_to_segment(p,point1,point2))
                !Go to top left corner
                do i=2,rr%nrules
                    point1 = point2
                    point2 = point(rr%lower(i),rr%s(i))
                    dist = min(dist,distance_to_segment(p,point1,point2))
                end do
                !Go to top right corner
                point1 = point2
                point2 = point(rr%upper(rr%nrules),rr%s(rr%nrules))
                dist = min(dist,distance_to_segment(p,point1,point2))
                !Go to bottom right corner
                do i=rr%nrules-1,2,-1
                    point1 = point2
                    point2 = point(rr%lower(i-1),rr%s(i-1))
                    dist = min(dist,distance_to_segment(p,point1,point2))
                end do
                return
            end if
        else !rr%method == 0
            if (rr%ixy==1) then !s is horizontal
                !Check if point is inside ruled rectangle
                !This also helps us halve the number of segments to check
                idx = binary_search(p%x,rr%s) 
                if (idx>1 .and. idx<=rr%nrules) then
                    !Could be inside ruled rectangle
                    yup = rr%upper(idx-1)
                    ylow = rr%lower(idx-1)
                    if (ylow <= p%y .and. p%y <= yup) then
                        dist = 0.d0 !Point is inside ruled rectangle
                        return
                    else if (ylow > p%y) then
                        !Iterate over lower segments
                        !Moving from left to right
                        do i=2,rr%nrules-1
                            !Two steps: go right, then go up/down
                            !Go right
                            point1 = point(rr%s(i-1),rr%lower(i-1))
                            point2 = point(rr%s(i),rr%lower(i-1))
                            dist = min(dist,distance_to_segment(p,point1,point2))
                            !Go up/down
                            point1 = point2
                            point2 = point(rr%s(i),rr%lower(i))
                            dist = min(dist,distance_to_segment(p,point1,point2))
                        end do
                        !Go to bottom right corner
                        point1 = point2
                        point2 = point(rr%s(rr%nrules),rr%lower(rr%nrules-1))
                        dist = min(dist,distance_to_segment(p,point1,point2))
                        return
                    else ! yup < p%y
                        !Iterate over upper segments
                        !Moving from left to right
                        do i=2,rr%nrules-1
                            !Two steps: go right, then go up/down
                            !Go right
                            point1 = point(rr%s(i-1),rr%upper(i-1))
                            point2 = point(rr%s(i),rr%upper(i-1))
                            dist = min(dist,distance_to_segment(p,point1,point2))
                            !Go up/down
                            point1 = point2
                            point2 = point(rr%s(i),rr%upper(i))
                            dist = min(dist,distance_to_segment(p,point1,point2))
                        end do
                        !Go to top right corner
                        point1 = point2
                        point2 = point(rr%s(rr%nrules),rr%upper(rr%nrules-1))
                        dist = min(dist,distance_to_segment(p,point1,point2))
                        return
                    end if
                end if
                !Not inside ruled rectangle, iterate over all segments
                !Start at bottom left corner and move clockwise
                !Go to top left corner
                point1 = point(rr%s(1),rr%lower(1))
                point2 = point(rr%s(1),rr%upper(1))
                dist = min(dist,distance_to_segment(p,point1,point2))
                !Start moving to the right
                do i=2,rr%nrules-1
                    !Two steps: go right, then go up/down
                    !Go right
                    point1 = point2
                    point2 = point(rr%s(i),rr%upper(i-1))
                    dist = min(dist,distance_to_segment(p,point1,point2))
                    !Go up/down
                    point1 = point2
                    point2 = point(rr%s(i),rr%upper(i))
                    dist = min(dist,distance_to_segment(p,point1,point2))
                end do
                !Go to top right corner
                point1 = point2
                point2 = point(rr%s(rr%nrules),rr%upper(rr%nrules-1))
                dist = min(dist,distance_to_segment(p,point1,point2))
                !Go to bottom right corner
                point1 = point2
                point2 = point(rr%s(rr%nrules),rr%lower(rr%nrules-1))
                dist = min(dist,distance_to_segment(p,point1,point2))
                !Start moving to the left
                do i=rr%nrules-1,2,-1
                    !Two steps: go left, then go up/down
                    !Go left
                    point1 = point2
                    point2 = point(rr%s(i),rr%lower(i))
                    dist = min(dist,distance_to_segment(p,point1,point2))
                    !Go up/down
                    point1 = point2
                    point2 = point(rr%s(i),rr%lower(i-1))
                    dist = min(dist,distance_to_segment(p,point1,point2))
                end do
                !Go to bottom left corner
                point1 = point2
                point2 = point(rr%s(1),rr%lower(1))
                dist = min(dist,distance_to_segment(p,point1,point2))
                return
            else !s is vertical
                !Check if point is inside ruled rectangle
                !This also helps us halve the number of segments to check
                idx = binary_search(p%y,rr%s)
                if(idx>1 .and. idx<=rr%nrules) then
                    !Could be inside ruled rectangle
                    xup = rr%upper(idx-1)
                    xlow = rr%lower(idx-1)
                    if (xlow <= p%x .and. p%x <= xup) then
                        dist = 0.d0 !Point is inside ruled rectangle
                        return
                    else if (xlow > p%x) then
                        !Iterate over left segments
                        !Moving from bottom to top
                        do i=1,rr%nrules-2
                            !Two steps: go up, then go left/right
                            !Go up
                            point1 = point(rr%lower(i),rr%s(i))
                            point2 = point(rr%lower(i),rr%s(i+1))
                            dist = min(dist,distance_to_segment(p,point1,point2))
                            !Go left/right
                            point1 = point2
                            point2 = point(rr%lower(i+1),rr%s(i+1))
                            dist = min(dist,distance_to_segment(p,point1,point2))
                        end do
                        !Go to top left corner
                        point1 = point2
                        point2 = point(rr%lower(rr%nrules-1),rr%s(rr%nrules))
                    else ! xup < p%x
                        !Iterate over right segments
                        !Moving from bottom to top
                        do i=1,rr%nrules-2
                            !Two steps: go up, then go left/right
                            !Go up
                            point1 = point(rr%upper(i),rr%s(i))
                            point2 = point(rr%upper(i),rr%s(i+1))
                            dist = min(dist,distance_to_segment(p,point1,point2))
                            !Go left/right
                            point1 = point2
                            point2 = point(rr%upper(i+1),rr%s(i+1))
                            dist = min(dist,distance_to_segment(p,point1,point2))
                        end do
                        !Go to top right corner
                        point1 = point2
                        point2 = point(rr%upper(rr%nrules-1),rr%s(rr%nrules))
                        dist = min(dist,distance_to_segment(p,point1,point2))
                    end if
                end if
                !Start at bottom right corner and move clockwise
                !Go to bottom left corner
                point1 = point(rr%upper(1),rr%s(1))
                point2 = point(rr%lower(1),rr%s(1))
                dist = min(dist,distance_to_segment(p,point1,point2))
                !Start moving up
                do i=2,rr%nrules-1
                    !Two steps: go up, then go left/right
                    !Go up
                    point1 = point2
                    point2 = point(rr%lower(i-1),rr%s(i))
                    dist = min(dist,distance_to_segment(p,point1,point2))
                    !Go left/right
                    point1 = point2
                    point2 = point(rr%lower(i),rr%s(i))
                    dist = min(dist,distance_to_segment(p,point1,point2))
                end do
                !Go to top left corner
                point1 = point2
                point2 = point(rr%lower(rr%nrules-1),rr%s(rr%nrules))
                dist = min(dist,distance_to_segment(p,point1,point2))
                !Go to top right corner
                point1 = point2
                point2 = point(rr%upper(rr%nrules-1),rr%s(rr%nrules))
                dist = min(dist,distance_to_segment(p,point1,point2))
                !Start moving down
                do i=rr%nrules-1,2,-1
                    !Two steps: go down, then go left/right
                    !Go down
                    point1 = point2
                    point2 = point(rr%upper(i),rr%s(i))
                    dist = min(dist,distance_to_segment(p,point1,point2))
                    !Go left/right
                    point1 = point2
                    point2 = point(rr%upper(i-1),rr%s(i))
                    dist = min(dist,distance_to_segment(p,point1,point2))
                end do
                !Go to bottom right corner
                point1 = point2
                point2 = point(rr%upper(1),rr%s(1))
                dist = min(dist,distance_to_segment(p,point1,point2))
                return
            end if
        end if

    end function distance_to_ruled_rectangle_exact

    pure function binary_search(id,arr) result(jloc)
        !! binary search of a sorted array with non-repeated
        !! elements (taking machine precision into account).
        !! Modified from the NASTRAN-95 code.
        !! returns the index jloc of an element such that
        !! arr(jloc-1)<id<=arr(jloc)
        !! if arr(1)>id, jloc=1
        !! if arr(n)<=id, jloc=n+1
        !! if an error occurs, jloc=-1
        implicit none

        double precision, intent(in) :: id
            !! key word to match in `arr`
        double precision, dimension(:),intent(in) :: arr
            !! array to search (it is
            !! assumed to be sorted)
        integer :: jloc
            !! the first matched index in 'arr'
            !! (if not found, 0 is returned)

        integer :: j,k,khi,klo,n
        integer,parameter :: iswtch = 16
        integer :: prevj

        n = size(arr)
        prevj = 0
        jloc = 0

        if (id<arr(1)) then
            jloc = 1
            return
        else if (id>=arr(n)) then
            jloc = n+1
            return
        end if

        if ( n<iswtch ) then
            ! sequential search more efficient
            do j = 1 , n
                jloc = j
                if (id<=arr(j)) then
                    return
                end if
            end do
            jloc = n+1
            return !Point is to the right of the array
        else

            klo = 1
            khi = n
            k = (klo+khi+1)/2
            do
                j = k
                if ( id<arr(j) ) then
                    khi = k
                else if ( id<=arr(j) .and. id>arr(j-1) ) then
                    jloc = j
                    return
                else
                    klo = k
                end if
                if ( khi-klo<1 ) then
                    jloc = -1
                    return ! error
                else if ( khi-klo==1 ) then
                    jloc = khi
                    return
                else
                    k = (klo+khi+1)/2
                end if
            end do

        end if

    end function binary_search

    pure function euclidean_distance(p1,p2) result(dist)
        ! Compute the Euclidean distance between two points
        implicit none
        type(point), intent(in) :: p1,p2
        real(kind=8) :: dist
        dist = sqrt((p1%x-p2%x)**2 + (p1%y-p2%y)**2)
    end function euclidean_distance

    pure function get_y_line(x,p1,p2) result(y)
        ! Compute the y-coordinate of a point on a line
        ! given two points on the line and the x-coordinate
        implicit none
        real(kind=8), intent(in) :: x
        type(point), intent(in) :: p1,p2
        real(kind=8) :: y
        y = p1%y + (x-p1%x)*(p2%y-p1%y)/(p2%x-p1%x)
    end function get_y_line

    pure function get_x_line(y,p1,p2) result(x)
        ! Compute the x-coordinate of a point on a line
        ! given two points on the line and the y-coordinate
        implicit none
        real(kind=8), intent(in) :: y
        type(point), intent(in) :: p1,p2
        real(kind=8) :: x
        x = p1%x + (y-p1%y)*(p2%x-p1%x)/(p2%y-p1%y)
    end function get_x_line

    function distance_to_segment(p,p1,p2) result(dist)
        ! Compute the distance between the point p and the line segment given
        ! by the points p1 and p2
        implicit none
        type(point), intent(in) :: p,p1,p2
        real(kind=8) :: dist

        !Locals
        type(point) :: p3
        !For coordinate_system == 2
        real(kind=8), dimension(2) :: latlon, xy
        type(point) :: pxy, p1xy, p2xy
        
        ! if (.FALSE.) then ! (coordinate_system==2) then
        !     !I am recomputing the xy coordinates of the ruled rectangle.
        !     !Far from optimal, but cleaner code for the moment.
        !     xy = latlon2xy([p%x,p%y],projection_center)
        !     pxy%x = xy(1)
        !     pxy%y = xy(2)
        !     xy = latlon2xy([p1%x,p1%y],projection_center)
        !     p1xy%x = xy(1)
        !     p1xy%y = xy(2)
        !     xy = latlon2xy([p2%x,p2%y],projection_center)
        !     p2xy%x = xy(1)
        !     p2xy%y = xy(2)
        !     p3 = closest_point_in_segment(pxy,p1xy,p2xy)
        !     latlon = xy2latlon([p3%x,p3%y],projection_center)
        !     dist = spherical_distance(p%x,p%y,latlon(1),latlon(2))
        ! else
            p3 = closest_point_in_segment(p,p1,p2)
            dist = euclidean_distance(p,p3)
        ! end if
        
    end function distance_to_segment

    pure function closest_point_in_segment(p,p1,p2) result(p3)
        ! Compute the point in the segment p1-p2 that is closest to p
        ! on the plane
        implicit none
        type(point), intent(in) :: p,p1,p2
        type(point) :: p3
        real(kind=8) :: px,py,norm,u,x,y
        px = p2%x-p1%x
        py = p2%y-p1%y
        norm = px*px + py*py
        u =  ((p%x - p1%x) * px + (p%y - p1%y) * py) / norm
        if (u > 1) then
            u = 1.0
        elseif (u < 0) then
            u = 0.0
        end if
        x = p1%x + u * px
        y = p1%y + u * py
        p3%x = x
        p3%y = y
    end function closest_point_in_segment

end module bouss_regions_module