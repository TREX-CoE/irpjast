program jastrow
  implicit none
  print *, 'Nabla J1'
  integer::k
  double precision :: j1, j2, j0, deriv, dt, lapl
  dt = 1.d-4

BEGIN_TEMPLATE 
lapl = 0.d0
j0 = $X $Y 
do k=1,3

  elec_coord(1,k) -= dt
  TOUCH elec_coord
  j1 = $X $Y 

  elec_coord(1,k) += 2.d0*dt
  TOUCH elec_coord
  j2 = $X $Y 

  deriv = (j2 - j1)/(2.d0*dt)
  lapl += (j2 - 2.d0*j0 + j1)/(dt*dt)
  print *, 'deriv $X '
  print *, deriv
  print *, $X_deriv_e(k,$Z)
  print *, ''

  elec_coord(1,k) -= dt
  TOUCH elec_coord

enddo
  print *, 'lapl $X '
  print *, lapl
  print *, $X_deriv_e(4 ,$Z)
  print *, ''

SUBST [X,Y,Z]
factor_een ;     ;  1   ;;
END_TEMPLATE
!rescale_een_e ; (1,3,1) ; 1,3,1 ;;
!rescale_een_n ; (1,1,2) ; 1,1,2 ;;
!rescale_een_e ; (1,2,2) ; 1,2,2 ;;
!elnuc_dist ; (1,1); 1,1 ;;
!elec_dist  ; (1,2); 1,2 ;;


end program
