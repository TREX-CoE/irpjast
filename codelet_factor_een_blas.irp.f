
program codelet_factor_een_blas
  implicit none
  integer :: i
  double precision :: ticks_0, ticks_1, cpu_0, cpu_1
  integer*8 :: irp_imax 


  PROVIDE factor_een_blas tmp_c 

  call provide_factor_een_blas

  double precision :: irp_rdtsc

  irp_imax = max(1_8,20_8 * 125000000_8 /(int(nelec,8) * int(nelec,8) * int(nnuc,8) * ncord))

  call cpu_time(cpu_0)
  ticks_0 = irp_rdtsc()
  do i=1,irp_imax
    call bld_tmp_c
    call bld_factor_een_blas
  enddo
  ticks_1 = irp_rdtsc()
  call cpu_time(cpu_1)
  print *, 'factor_een_blas'
  print *, '-----------'
  print *, 'Cycles:'
  print *,  (ticks_1-ticks_0)/dble(irp_imax)
  print *, 'Seconds:'
  print *,  (cpu_1-cpu_0)/dble(irp_imax)
end

  
