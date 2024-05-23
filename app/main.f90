program main
  use fitpack_grid_surfaces
  use Utils
  implicit none
  integer :: L, i
  real(8) :: V0, T1, T2, Vinc
  character(len=80) :: filename
  real(8), allocatable :: vol(:,:), grid(:,:)
  
  call parse_command_line(L, V0, filename)
  allocate(vol(2*L,L), grid(2*L,2*L))
  
  call cpu_time(T1)
  !vol = interpolate(L, V0)
  Vinc=V0/L
  vol=0.0d0

  do i=0, L-1
    vol(L+i, 1) = V0-Vinc*i
  end do

  call relax(vol, L)
  call cpu_time(T2)

  grid = combine(vol)
  print *, T2-T1
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
  
    if (L <= 20) then
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
    real(8), target, intent(inout) :: vol(:, :)
    integer, intent(in) :: L
    real(8), target, dimension(2*L,L) :: next_vol, diff
    real(8), pointer :: old(:,:), new(:,:), tmp(:,:)
    real(8) :: eps=1E-5

    next_vol = vol
    old => vol
    new => next_vol
    diff = 1.0d0
    !call Print2D(diff)
    !print *, maxval(abs(diff))
    !print *, maxval(abs(diff)) < eps

    do while (maxval(abs(diff)) > eps) !(n < 5)!
      !call Print2D(old)
      new(2:2*L-1, 2:L-1) = (old(2:2*L-1, 1:L-2) + old(2:2*L-1, 3:L) + old(1:2*L-2, 2:L-1) + old(3:2*L, 2:L-1)) / 4.0d0 
      new(2:L-1,1) = new(L,L-1:2:-1)
      diff = new-old
      tmp => new
      new => old
      old => tmp
    end do
    vol = new
  end subroutine relax

end program main

