program main
  use fitpack_grid_surfaces
  use Utils
  implicit none
  integer :: L
  real(8) :: V0
  character(len=80) :: filename
  real(8), allocatable :: vol(:,:), grid(:,:)
  
  call parse_command_line(L, V0, filename)
  allocate(vol(2*L,L), grid(2*L,2*L))
  
  vol = interpolate(L, V0)

  grid = combine(vol)
  call write2file(filename, grid)
  deallocate(vol, grid)
contains

  recursive function interpolate(L, V0) result(vol)

    integer, intent(in) :: L
    real(8), intent(in) :: V0
    integer :: i, j, ierr, Lc
    real(8) :: Vinc
    real(8), allocatable :: vol(:,:), coarse(:,:)
    type(fitpack_grid_surface) :: surface
    Lc = L/2
    allocate(vol(2*L,L), coarse(2*Lc,2*Lc))
  
    if (L <= 50) then
      Vinc=V0/L
      vol=0.0d0

      do i=0, L-1
        vol(L+i, 1) = V0-Vinc*i
      end do

      call relax(vol, L)
    else
      Vinc=V0/Lc
      coarse=interpolate(Lc, V0)
      
      do i=2, Lc-2
        coarse(2:Lc-i, i) = coarse(Lc+1-i, Lc-1:i+1:-1)
      end do
      ierr = surface%new_fit([(real(i, kind=8)/(Lc-1), i=0, Lc-1)], [(real(i, kind=8)/(Lc-1/2), i=0, 2*Lc-1)], coarse, 0.0d0, 3)
      vol = surface%eval( [(real(i, kind=8)/(real(L, kind=8)-1.0), i=0, L-1)], &
                          [(real(i, kind=8)/(real(L, kind=8)-0.5), i=0, 2*L-1)], ierr)
      vol(:, [1,L]) = 0.0d0
      vol([1, 2*L], :) = 0.0d0
      do j=2,L-2
        do i=2, L-j
          vol(i,j) = 0.0d0
        end do
      end do
      Vinc=V0/L
      do i=0, L-1
        vol(L+i, 1) = V0-Vinc*i
      end do
      call relax(vol, L)
  end if

  end function interpolate

  subroutine relax(vol, L)
    real(8), intent(inout) :: vol(:, :)
    integer, intent(in) :: L
    integer :: i, j, k, n=0
    logical :: loop=.true.
    real(8) :: new

    do while (loop)
      loop=.false.
      do j=2,L-1
        do i=L+2-j, 2*L-1
          new = (vol(i-1,j) + vol(i+1,j) + vol(i,j-1) + vol(i,j+1)) / 4.0
          if (abs(new - vol(i,j)) > 1E-5) then
            loop=.true.
          end if
          vol(i,j) = new
        end do
        new = (2*vol(L+2-j,j) + 2*vol(L+1-j,j+1)) / 4.0
        if (abs(new - vol(L+1-j,j)) > 1E-5) then
          loop=.true.
        end if
        vol(L+1-j,j) = new
      end do
    end do
  end subroutine relax

end program main

