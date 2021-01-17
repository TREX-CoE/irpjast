program jastrow_irp
  implicit none
  TOUCH elec_coord
  print *, "J_{IRP} = ", jastrow_full
  print *, "\nabla J_{IRP} = ", jastrow_derivs(1:3, :)
  print *, "\nabla^2 J_{IRP} = ", jastrow_derivs(4, :)
end program
