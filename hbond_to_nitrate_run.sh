gfortran -o achbond_k_pair.x achbond_k_pair.f90
cp hbond_to_nitrate_input achbond_k_pair_input
./achbond_k_pair.x
xmgrace achbond_k_pair_acf_h.dat achbond_k_pair_acf_h_d.dat |
xmgrace achbond_k_pair_ln_acf_h.dat achbond_k_pair_ln_acf_h_d.dat |
xmgrace achbond_k_pair_acf_hh.dat achbond_k_pair_acf_hh_d.dat |
xmgrace achbond_k_pair_ln_acf_hh.dat achbond_k_pair_ln_acf_hh_d.dat
