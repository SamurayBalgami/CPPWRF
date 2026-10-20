! wrf_error_fatal3 (frame/module_wrf_error.F) for the three WRF phys objects linked by run_test.sh
subroutine wrf_error_fatal3(file_str, line, str)
  implicit none
  character(len=*), intent(in) :: file_str, str
  integer, intent(in) :: line
  write(*, *) 'FATAL ERROR: ', trim(file_str), ':', line, ': ', trim(str)
  stop 1
end subroutine wrf_error_fatal3
