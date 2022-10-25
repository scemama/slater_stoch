subroutine write_STO_in_file_info_basis(i_orb)
  use common_data
  integer, intent(in) :: i_orb

  integer :: mm

  write(21,'(a,i4,5a,i4,a)')'Orb #',i_orb,'  ',trim(orb_name_full(i_orb)),'  ',trim(ATOM(nucleus_number(i_orb))),&
        ' STO fitted with',ng0,' gaussians'
  write(21,'(a,i2)')'Number of contractions= ',n_contract(i_orb)
  do mm=1,n_contract(i_orb)
    write(21,'(a5,f10.3,a6,f10.3)')'  c_i= ',c_contract(mm,i_orb),' g_i= ',g_contract(mm,i_orb)
  enddo
  write(21,*)
  write(22,'(i5, x, i5)')i_orb,n_gauss(i_orb)
  write(22,'(4e22.15)')(c_gauss(mm,i_orb,1),mm=1,n_gauss(i_orb))
  write(22,'(4e22.15)')(c_gauss(mm,i_orb,2),mm=1,n_gauss(i_orb))
  write(22,'(4e22.15)')(g_gauss(mm,i_orb)  ,mm=1,n_gauss(i_orb))
  write(22,*)
end

