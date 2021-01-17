      subroutine jast_elec_champ(nelec, elec_coord)
	! This subroutine allows for a correct interfacing between 
	! the Jastrow IRPF90 files and the CHAMP variables 
        implicit real*8(a-h, o-z)
	! This files must be included when compiling in CHAMP
	include 'vmc.h'
	include 'force.h'
        ! This common blocks are defined in CHAMP
        common /config/ xold(3, MELEC), xnew(3, MELEC), vold(3, MELEC)

	! Electron coordinates
        dimension elec_coord(nelec, 3)

        do j = 1, 3
           do i = 1, nelec
              !elec_coord(i, j) = xold(j, i)
              ! Temporarily set some values for testing without champ
              elec_coord(i, j) = dble(i + j) * 0.1d0
           end do
        end do

      end subroutine jast_elec_champ
 
      subroutine jast_nuc_champ(nnuc, xnuc_coord)
	! This subroutine allows for a correct interfacing between 
	! the Jastrow IRPF90 files and the CHAMP variables 
        implicit real*8(a-h, o-z)
	! This files must be included when compiling in CHAMP
	include 'vmc.h'
	include 'force.h'
        ! This common blocks are defined in CHAMP
        common/atom/znuc(MCTYPE), cent(3, MCENT), pecent

	! Nuclear coordinates
        dimension xnuc_coord(nnuc, 3)

        do j = 1, 3
           do i = 1, nnuc
              !xnuc_coord(i, j) = cent(j, i)
              ! Temporarily set some values for testing without champ
              xnuc_coord(i, j) = dble(i + j) * 0.01d0
           end do
        end do                          

      end subroutine jast_nuc_champ

      subroutine jast_pars_champ(naord, ntypenuc, aord_vect, nbord,
     &                           bord_vect, ndim_cord_vect, cord_vect)
	! This subroutine allows for a correct interfacing between 
	! the Jastrow IRPF90 files and the CHAMP variables 
        implicit real*4(a-h, o-z)
	! This files must be included when compiling in CHAMP
	include 'vmc.h'
        include 'force.h'
        ! This common blocks are defined in CHAMP
        common/jaspar3/a(MORDJ1,MWF),b(MORDJ1,2,MWF),c(83,MCTYPE,MWF)
        common/jaspar4/a4(MORDJ1,MCTYPE,MWF),norda,nordb,nordc

	! Jastrow parameters
        double precision aord_vect(naord + 1, ntypenuc)
        double precision bord_vect(nbord + 1)
        double precision cord_vect(ndim_cord_vect, ntypenuc)

        do j = 1, ntypenuc
           do i = 1, naord + 1
              !aord_vect(i, j) = a4(i, j, 1)
              ! Temporarily set some values for testing without champ
              aord_vect(i, j) = 1.0d0
           end do
        end do                          

        do i = 1, nbord + 1
           !bord_vect(i) = b(i, 1, 1)
           ! Temporarily set some values for testing without champ
           bord_vect(i) = 0.5d0
        end do

        do j = 1, ntypenuc
           do i = 1, ndim_cord_vect
              !cord_vect(i, j) = c(i, j, 1)
              ! Temporarily set some values for testing without champ
              cord_vect(i, j) = 1.0d0
           end do
        end do

      end subroutine jast_pars_champ
