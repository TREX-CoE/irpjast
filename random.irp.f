BEGIN_PROVIDER [ integer, seed_size ]
 implicit none
 BEGIN_DOC
 ! Size of the random seed
 END_DOC
  call random_seed(size=seed_size)

END_PROVIDER

BEGIN_PROVIDER [ integer, seed, (seed_size) ]
 implicit none
 BEGIN_DOC
 ! Random seed
 END_DOC
 integer :: i
 do i=1,seed_size
   seed(i) = i
 enddo
 call random_seed(put=seed)

END_PROVIDER
