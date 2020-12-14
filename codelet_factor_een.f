
program codelet_factor_een
  implicit none
  integer :: i
  double precision :: ticks_0, ticks_1, cpu_0, cpu_1
  integer, parameter :: irp_imax = 100000

  

  call provide_factor_een

  double precision :: irp_rdtsc

  call cpu_time(cpu_0)
  ticks_0 = irp_rdtsc()
  do i=1,irp_imax
    call bld_factor_een
  enddo
  ticks_1 = irp_rdtsc()
  call cpu_time(cpu_1)
  print *, 'factor_een'
  print *, '-----------'
  print *, 'Cycles:'
  print *,  (ticks_1-ticks_0)/dble(irp_imax)
  print *, 'Seconds:'
  print *,  (cpu_1-cpu_0)/dble(irp_imax)
end

  