module dgl_utils_c
    use dgl_interface

    integer, parameter :: dp = dgl_real

contains

subroutine C_string_ptr_to_F_string(C_string, F_string)
    use ISO_C_BINDING
    type(C_PTR), intent(in) :: C_string
    character(len=*), intent(out) :: F_string
    character(len=1, kind=C_CHAR), dimension(:), pointer :: p_chars
    integer :: i
    if (.not. C_associated(C_string)) then
      F_string = ' '
    else
      call C_F_pointer(C_string, p_chars, [huge(0)])
      do i = 1, len(F_string)
        if (p_chars(i) == C_NULL_CHAR) exit
        F_string(i:i) = p_chars(i)
      end do
      if (i <= len(F_string)) F_string(i:) = ' '
    end if
end subroutine

end module dgl_utils_c