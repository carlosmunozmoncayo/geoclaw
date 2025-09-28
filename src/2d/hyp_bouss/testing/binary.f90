pure function binary_search(id,arr) result(jloc)
    !! binary search of a sorted array with non-repeated
    !! elements (taking machine precision into account).
    !! returns the index jloc of an element such that
    !! arr(jloc-1)<id<=arr(jloc)
    !! if arr(1)>id, jloc=0
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