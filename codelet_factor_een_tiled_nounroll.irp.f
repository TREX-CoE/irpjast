
program codelet_factor_een_tiled_nounroll
  implicit none
  integer :: i
  double precision :: ticks_0, ticks_1, cpu_0, cpu_1
  integer*8 :: irp_imax 


  !PROVIDE factor_een_blas tmp_c_tiled
  PROVIDE factor_een_tiled_nounroll tmp_c_tiled_nounroll

  call provide_factor_een_tiled_nounroll

  double precision :: irp_rdtsc

  irp_imax = max(1_8,20_8 * 125000000_8 /(int(nelec,8) * int(nelec,8) * int(nnuc,8) * ncord))
  irp_imax = 1
  print *,"irp_imax=",irp_imax

  call cpu_time(cpu_0)
  ticks_0 = irp_rdtsc()
  do i=1,irp_imax
  print *,i
    !call bld_tmp_c_tiled
    call bld_tmp_c_tiled_nounroll
    call bld_factor_een_tiled_nounroll
  enddo
  ticks_1 = irp_rdtsc()
  call cpu_time(cpu_1)
  print *, 'factor_een_tiled_nounroll'
  print *, '-----------'
  print *, 'Cycles:'
  print *,  (ticks_1-ticks_0)/dble(irp_imax)
  print *, 'Seconds:'
  print *,  (cpu_1-cpu_0)/dble(irp_imax)
end

  