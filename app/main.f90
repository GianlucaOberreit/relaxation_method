program main
  use fitpack_grid_surfaces
  use Utils
  implicit none
  integer :: L, i
  real(8) :: V0, Vinc
  character(len=80) :: filename
  real(8), allocatable :: vol(:,:), grid(:,:)
  
  call parse_command_line(L, V0, filename)
  allocate(vol(2*L,L), grid(2*L,2*L))
  
  vol = interpolate(L)

  grid = combine(vol)
  call write2file(filename, grid)

  deallocate(vol, grid)
contains

  recursive function interpolate(L) result(vol)

    integer, intent(in) :: L
    integer :: i, Lc
    real(8) :: Vinc
    real(8), allocatable :: vol(:,:), coarse(:,:)
    Lc = L/2
    allocate(vol(2*L,L), coarse(2*Lc,2*Lc))
  
    if (L <= 20) then
      Vinc=V0/L
      vol=0.0d0

      do i=0, L-1
        vol(L+i, 1) = V0-Vinc*i
      end do

      call relax(vol, L)
    else
      Vinc=V0/Lc
      coarse=interpolate(Lc)
      
      vol = fit(L, Lc, coarse)

      vol(:, [1,L]) = 0.0d0
      vol([1, 2*L], :) = 0.0d0
      Vinc=V0/L
      do i=0, L-1
        vol(L+i, 1) = V0-Vinc*i
      end do

      call relax(vol, L)
  end if

  end function interpolate

  subroutine relax(vol, L)
    real(8), target, intent(inout) :: vol(:, :)
    integer, intent(in) :: L
    real(8), target, allocatable, dimension(:,:) :: next_vol, diff
    real(8), pointer :: old(:,:), new(:,:), tmp(:,:)
    real(8) :: eps=1E-5
    allocate(next_vol(2*L,L), diff(2*L,L))

    next_vol = vol
    old => vol
    new => next_vol
    diff = 1.0d0

    do while (maxval(abs(diff)) > eps) !(n < 5)!
      new(2:2*L-1, 2:L-1) = (old(2:2*L-1, 1:L-2) + old(2:2*L-1, 3:L) + old(1:2*L-2, 2:L-1) + old(3:2*L, 2:L-1)) / 4.0d0 
      new(2:L-1,1) = new(L,L-1:2:-1)
      diff = new-old
      tmp => new
      new => old
      old => tmp
    end do
    vol = new
  end subroutine relax

  function fit(L, Lc, coarse) result(vol)
    integer, intent(in) :: L, Lc
    real(8), intent(in) :: coarse(:,:)
    integer :: ierr
    real(8) :: vol(2*L,L)
    type(fitpack_grid_surface) :: surface

    ierr = surface%new_fit([(real(i, kind=8)/(Lc-1), i=0, Lc-1)], [(real(i, kind=8)/(Lc-1./2.), i=0, 2*Lc-1)], coarse, 0.0d0, 3)
    surface%lwrk = 2 * (2*L*4 + L*4)
    surface%liwrk = 2*(2*L + L)
    deallocate(surface%wrk, surface%iwrk)
    allocate(surface%wrk(surface%lwrk), surface%iwrk(surface%liwrk))
    vol(2:2*L-1,2:L-1) = surface%eval( [(real(i, kind=8)/(real(L-3, kind=8)-1.0), i=0, L-3)], &
                        [(real(i, kind=8)/(real(L-3./2., kind=8)-0.5), i=0, 2*L-3)], ierr)
  end function fit

end program main

