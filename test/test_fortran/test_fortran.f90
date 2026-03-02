program test_fortran
    use dgl_interface
    use solvers
    implicit none

    call test_davidson(500, 10, 10, verbose=.true.)

    stop "testing new tests"
!  allocate (eig(n_max), evec(n, n_max))
!!
!! test davidson:
!!
!  eig  = dgl_zero
!  evec = dgl_zero
!  do i = 1, n_max
!    evec(i,i) = dgl_one
!  end do
!  ok = .false.
!!
!  mx_p => mx
!  write(6,*) ' testing Davidson:'
!  call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
!                        dgl_verbose = .true.,&
!                        dgl_memory = memory,&
!                        dgl_memory_unit = memory_unit)
!!
!  if (ok) then
!    write(6,*) ' Davidson converged.'
!    write(lutest,1000) 'Davidson'
!    write(lutest,*)
!    write(lutest,1010)
!    do i = 1, n_targ
!      write(lutest,1020) i, eig(i)
!    end do
!!
!!     fix the phase
!!
!    do i = 1, n_targ
!      if (evec(1,i).lt.dgl_zero) evec(:,i) = - evec(:,i)
!    end do
!    write(lutest,1031) (i, i = 1, n_targ)
!    do j = 1, n
!      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
!    end do
!    write(lutest,*)
!  else
!    write(6,*) ' Davidson failed to converge.'
!  end if
!!
!! test generalized davidson:
!!
!  eig  = dgl_zero
!  evec = dgl_zero
!  do i = 1, n_max
!    evec(i,i) = dgl_one
!  end do
!  ok = .false.
!!
!  write(6,*) ' testing Generalized Davidson:'
!  call dgl_davidson_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, metvec=mx_p, &
!                        dgl_verbose = .true.,&
!                        dgl_memory = memory,&
!                        dgl_memory_unit = memory_unit)
!!
!  if (ok) then
!    write(6,*) ' Generalized Davidson converged.'
!    write(lutest,1000) 'Generalized Davidson'
!    write(lutest,*)
!    write(lutest,1010)
!    do i = 1, n_targ
!      write(lutest,1020) i, eig(i)
!    end do
!    do i = 1, n_targ
!!
!! fix the phase
!!
!  if (evec(1,i).lt.dgl_zero) evec(:,i) = - evec(:,i)
!    end do
!    write(lutest,1031) (i, i = 1, n_targ)
!    do j = 1, n
!      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
!    end do
!    write(lutest,*)
!  else
!    write(6,*) 'Generalized Davidson failed to converge.'
!  end if
!!
!! test non-symmetric davidson:
!!
!  allocate (evec_l(n,n_max))
!  evec = dgl_zero
!  do i = 1, n_max
!    evec(i,i) = dgl_one
!  end do
!  ok = .false.
!!
!  write(6,*) ' testing non-symmetric Davidson:'
!  call dgl_davidson_nosym_driver(n, n_targ, n_max, arx, alx, dx, "LR", eig, evec, ok, evec_2 = evec_l, &
!                              dgl_verbose = .true.,&
!                              dgl_memory = memory, &
!                              dgl_memory_unit = memory_unit)
!!
!  if (ok) then
!    write(6,*) ' non-symmetric Davidson converged.'
!    write(lutest,1000) 'Non-Symmetric Davidson'
!    write(lutest,*)
!    write(lutest,1010)
!
!    do i = 1, n_targ
!      write(lutest,1020) i, eig(i)
!    end do
!
!    write(lutest,1032) 'Right ', (i, i = 1, n_targ)
!!
!!   fix the phase
!!
!    if (evec(1,i).lt.dgl_zero) evec(:,i) = - evec(:,i)
!    do j = 1, n
!      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
!    end do
!    write(lutest,*)
!
!    write(lutest,1032) 'Left ', (i, i = 1, n_targ)
!!
!!     fix the phase
!!
!    if (evec_l(1,i).lt.dgl_zero) evec_l(:,i) = - evec_l(:,i)
!    do j = 1, n
!      write(lutest,1022) j, (evec_l(j,i), i = 1, n_targ)
!    end do
!    write(lutest,*)
!
!  else
!    write(6,*) ' non-symmetric Davidson failed to converge.'
!  end if
!!
!  deallocate (evec_l)
!!
!! test lobpcg:
!!
!  eig  = dgl_zero
!  evec = dgl_zero
!!
!  do i = 1, n_max
!    evec(i,i) = dgl_one
!  end do
!  ok = .false.
!!
!  write(6,*) ' testing LOBPCG:'
!
!  call dgl_lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, &
!                        dgl_verbose = .true.,&
!                        dgl_memory = memory,&
!                        dgl_memory_unit = memory_unit)
!!
!  if (ok) then
!    write(6,*) ' LOBPCG converged.'
!    write(lutest,1000) 'LOBPCG'
!    write(lutest,*)
!    write(lutest,1010)
!    do i = 1, n_targ
!      write(lutest,1020) i, eig(i)
!    end do
!    do i = 1, n_targ
!!
!!     fix the phase
!!
!  if (evec(1,i).lt.dgl_zero) evec(:,i) = - evec(:,i)
!    end do
!    write(lutest,1031) (i, i = 1, n_targ)
!    do j = 1, n
!      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
!    end do
!    write(lutest,*)
!  else
!    write(6,*) ' LOBPCG failed to converge.'
!  end if
!!
!! test lobpcg:
!!
!  eig  = dgl_zero
!  evec = dgl_zero
!!
!  do i = 1, n_max
!    evec(i,i) = dgl_one
!  end do
!  ok = .false.
!!
!  write(6,*) ' testing Generalized LOBPCG:'
!
!  call dgl_lobpcg_driver(n, n_targ, n_max, ax, dx, eig, evec, ok, metvec=mx_p, &
!                        dgl_verbose = .true.,&
!                        dgl_memory = memory,&
!                        dgl_memory_unit = memory_unit)
!!
!  if (ok) then
!    write(6,*) ' Generalized LOBPCG converged.'
!    write(lutest,1000) 'Generalized LOBPCG'
!    write(lutest,*)
!    write(lutest,1010)
!    do i = 1, n_targ
!      write(lutest,1020) i, eig(i)
!    end do
!    do i = 1, n_targ
!!
!!     fix the phase
!!
!  if (evec(1,i).lt.dgl_zero) evec(:,i) = - evec(:,i)
!    end do
!    write(lutest,1031) (i, i = 1, n_targ)
!    do j = 1, n
!      write(lutest,1022) j, (evec(j,i), i = 1, n_targ)
!    end do
!    write(lutest,*)
!  else
!    write(6,*) ' Generalized LOBPCG failed to converge.'
!  end if
!!
!! test smogd:
!!
!  deallocate (evec)
!  allocate (evec(2*n, n_max))
!  eig  = dgl_zero
!  evec = dgl_zero
!!
!  do i = 1, n_max
!    evec(i,i)   = dgl_one
!  end do
!  ok = .false.
!!
!  write(6,*) ' testing SMOGD:'
!  call dgl_smogd_driver(2*n, n_targ, n_max, apbx, ambx, spdx, smdx, lrprc, &
!                        eig, evec, ok, &
!                        dgl_verbose = .true., &
!                        dgl_memory = memory, &
!                        dgl_memory_unit = memory_unit)
!  if (ok) then
!    write(6,*) ' SMOGD converged.'
!    write(lutest,1000) 'SMOGD'
!    write(lutest,*)
!    write(lutest,1010)
!    do i = 1, n_targ
!      write(lutest,1020) i, eig(i)
!    end do
!    write(lutest,*)
!  else
!    write(6,*) ' SMOGD failed to converge.'
!  end if
!!
!! close the output file:
!!
!  close (lutest)
!!
!! free the memory:
!!
!  deallocate (evec, eig)
!!
end program test_fortran
