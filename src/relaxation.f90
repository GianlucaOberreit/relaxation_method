module relaxation
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, relaxation!"
  end subroutine say_hello
end module relaxation
