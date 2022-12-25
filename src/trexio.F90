#ifdef HAVE_TREXIO
#include <trexio_f.f90>

subroutine trexio_write_geometry(trexio_file)
  use trexio
  use common_data
  implicit none
  integer(trexio_t), intent(in)  :: trexio_file
  integer                        :: rc

  rc = trexio_write_nucleus_num(trexio_file, number_nuclei)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_nucleus_charge(trexio_file, charge)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_nucleus_coord(trexio_file, centers_nuclei)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_nucleus_label(trexio_file, ATOM, len(ATOM(1))+1)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_nucleus_repulsion(trexio_file, enucl)
  call trexio_assert(rc, TREXIO_SUCCESS)

end

subroutine trexio_write_basis(trexio_file)
  use trexio
  use common_data
  implicit none
  integer(trexio_t), intent(in)  :: trexio_file
  integer                        :: rc
  integer                        :: i, m, shell_num, prim_num

  integer, allocatable :: nucleus_index(:), shell_ang_mom(:), shell_index(:)
  integer, allocatable :: r_power(:), ao_shell(:)
  double precision, allocatable :: exponent(:), coefficient(:), shell_factor(:)
  double precision, allocatable :: prim_factor(:), ao_normalization(:)

  rc = trexio_write_basis_type(trexio_file, 'S', 1)
  call trexio_assert(rc, TREXIO_SUCCESS)

  shell_num = 0
  do i=1,nbasis
    shell_num = shell_num+1
    select case (orb_name_full(i))
      case ('1S')
      case ('2S')
      case ('3S')
      case ('2P_Z')
      case ('3P_Z')
      case ('3D_ZZ')
      case ('4F_ZZZ')
      case ('5G_ZZZZ')
      case default
        shell_num = shell_num-1
    end select
  end do
  rc = trexio_write_basis_shell_num(trexio_file, shell_num)
  call trexio_assert(rc, TREXIO_SUCCESS)

  prim_num = shell_num

  rc = trexio_write_basis_prim_num(trexio_file, prim_num)
  call trexio_assert(rc, TREXIO_SUCCESS)

  allocate(nucleus_index(shell_num), shell_ang_mom(shell_num),       &
      shell_factor(shell_num), r_power(shell_num), ao_shell(nbasis))
  allocate(shell_index(prim_num), exponent(prim_num), coefficient(prim_num),&
      prim_factor(prim_num))

  m = 0
  do i=1,nbasis
    m = m+1
    ao_shell(i) = m
    select case (orb_name_full(i))
      case ('1S')
        shell_ang_mom(m) = 0
        r_power(m) = 0
      case ('2S')
        shell_ang_mom(m) = 0
        r_power(m) = 1
      case ('3S')
        shell_ang_mom(m) = 0
        r_power(m) = 2
      case ('2P_Z')
        shell_ang_mom(m) = 1
        r_power(m) = 1
      case ('3P_Z')
        shell_ang_mom(m) = 1
        r_power(m) = 2
      case ('3D_ZZ')
        shell_ang_mom(m) = 2
        r_power(m) = 2
      case ('4F_ZZZ')
        shell_ang_mom(m) = 3
        r_power(m) = 3
      case ('5G_ZZZZ')
        shell_ang_mom(m) = 4
        r_power(m) = 4
      case default
        m = m-1
        cycle
    end select
    nucleus_index(m) = nucleus_number(i)
    exponent(m) = g_contract(1,i)
  end do

  do m=1,shell_num
     shell_index(m) = m
     shell_factor(m) = 1.d0
     coefficient(m) = 1.d0
     prim_factor(m) = 1.d0
  end do

  rc = trexio_write_basis_nucleus_index(trexio_file, nucleus_index)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_basis_shell_ang_mom(trexio_file, shell_ang_mom)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_basis_shell_factor(trexio_file, shell_factor)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_basis_r_power(trexio_file, r_power)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_basis_shell_index(trexio_file, shell_index)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_basis_exponent(trexio_file, exponent)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_basis_coefficient(trexio_file, coefficient)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_basis_prim_factor(trexio_file, prim_factor)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_ao_cartesian(trexio_file, 1)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_ao_num(trexio_file, nbasis)
  call trexio_assert(rc, TREXIO_SUCCESS)

  rc = trexio_write_ao_shell(trexio_file, ao_shell)
  call trexio_assert(rc, TREXIO_SUCCESS)

  allocate(ao_normalization(nbasis))
  ao_normalization = 1.d0

  rc = trexio_write_ao_normalization(trexio_file, ao_normalization)
  call trexio_assert(rc, TREXIO_SUCCESS)

end

#endif

