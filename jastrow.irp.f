program jastrow
  implicit none
  print *, 'Number of electrons: ', nelec
  print *, 'The total Jastrow factor'
  print *, jastrow_full
  print *, 'REF'
  print *, factor_een_deriv_e_ref
  print *, 'X'
  print *, factor_een_deriv_e
  print *, 'BLAS'
  print *, factor_een_deriv_e_blas(1:nelec,1:4)
  print *, ''
  print *, factor_een_deriv_e_ref(1:nelec,1) - factor_een_deriv_e_blas(1:nelec,1)
  print *, factor_een_deriv_e_ref(1:nelec,2) - factor_een_deriv_e_blas(1:nelec,2)
  print *, factor_een_deriv_e_ref(1:nelec,3) - factor_een_deriv_e_blas(1:nelec,3)
  print *, factor_een_deriv_e_ref(1:nelec,4) - factor_een_deriv_e_blas(1:nelec,4)
  !PROVIDE jastrow_full

end program
