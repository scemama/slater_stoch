!! derivative of u_orb
double precision function up_orb(i_orb,r)
include 'j.inc'
up_orb=0.d0
 do mm=1,n_contract(i_orb)
  contrib=(n_sto(i_orb)/r-g_contract(mm,i_orb))*c_contract(mm,i_orb)*r**n_sto(i_orb)*dexp(-g_contract(mm,i_orb)*r   )
  up_orb=up_orb+contrib
 enddo  !mm
end

