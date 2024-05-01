module Utils
  use, intrinsic :: iso_fortran_env
  use stdlib_io_npy, only: save_npy
  implicit none

  public

contains

  subroutine Print2D(array)
    implicit none
    real(8), intent(in) :: array(:,:)
    integer :: i, j

    do i = 1, size(array,1)
       do j = 1, size(array,2)
          write(*, '(F10.2)', advance='no') array(i, j)
       end do
       write(*,*) ! Print new line after each row
    end do
    write(*,*)
  end subroutine Print2D

  subroutine write2file(filename, array)
    implicit none
    character(len=*), intent(in) :: filename
    real(8), dimension(:,:), intent(in) :: array

    call save_npy(filename, array)
  end subroutine write2file

  subroutine parse_command_line(L, V0, filename)
    implicit none
    integer, intent(out) :: L
    real(8), intent(out) :: V0
    character(len=*) :: filename
    integer :: arg_count, ierr
    character(len=256) :: arg_value

    ! Get the number of command-line arguments
    arg_count = command_argument_count()

    ! Check if the number of arguments is correct
    if (arg_count /= 3) then
        print *, "Usage: ./relaxation <L> <V0> <filename.npy>"
        stop
    end if

    ! Read the first argument (L)
    call get_command_argument(1, arg_value)
    read(arg_value, *, iostat=ierr) L
    if (ierr /= 0) then
        print *, "Error: Invalid argument for L."
        stop
    end if

    ! Read the second argument (V0)
    call get_command_argument(2, arg_value)
    read(arg_value, *, iostat=ierr) V0
    if (ierr /= 0) then
        print *, "Error: Invalid argument for V0."
        stop
    end if

    ! Read the third argument (filename)
    call get_command_argument(3, arg_value)
    read(arg_value, '(A)', iostat=ierr) filename
    if (ierr /= 0) then
        print *, "Error: Invalid argument for filename."
        stop
    end if
  end subroutine parse_command_line
  
  function combine(vol) result(grid)
    integer :: i,j, L
    real(8), intent(in) :: vol(:,:)
    real(8), dimension(size(vol,1), size(vol,1)) :: grid

    L = size(vol,2)
    do j=1, L
      do i=L+1-j, 2*L
        grid(i, L+j) = vol(i,j)
        grid(L+1-j, 2*L+1-i) = vol(i,j)
      end do
    end do

    grid(L+1:size(vol,1), 1:size(vol,2)) = 0.0d0
  end function combine
end module Utils

