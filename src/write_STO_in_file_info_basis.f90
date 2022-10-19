subroutine write_STO_in_file_info_basis(i_orb)
  include 'j.inc'

  write(21,'(a,i4,5a,i4,a)')'Orb #',i_orb,'  ',trim(orb_name_full(i_orb)),'  ',trim(ATOM(nucleus_number(i_orb))),&
        ' STO fitted with',ng0,' gaussians'
  write(21,'(a,i2)')'Number of contractions= ',n_contract(i_orb)
  rcmax=-1.d0
  do mm=1,n_contract(i_orb)
    write(21,'(a5,f10.3,a6,f10.3)')'  c_i= ',c_contract(mm,i_orb),' g_i= ',g_contract(mm,i_orb)
  enddo
  write(21,'(a,e22.15)')'chi2 (STO orbital and its gaussian representation)= ',chi2(nbasis)
  write(21,*)
  write(22,'(i4,i3)')i_orb,n_gauss(i_orb)
  write(22,'(4e22.15)')(c_gauss(1,mm,i_orb),mm=1,n_gauss(i_orb))
  write(22,'(4e22.15)')(c_gauss(2,mm,i_orb),mm=1,n_gauss(i_orb))
  write(22,'(4e22.15)')(g_gauss(mm,i_orb)  ,mm=1,n_gauss(i_orb))
  write(22,*)
end

