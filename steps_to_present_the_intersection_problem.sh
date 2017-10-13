#

cp Thermo_functions_with_P_and_G_sorted_with_EL_from_EOS__VOLUME_and_All_energies_divided_by_F_unit_SORTED_in_Ts_blocked_with_header.dat dd.dat

rm -rf to_post_math_problem

mkdir to_post_math_problem

mv dd.dat ./to_post_math_problem

cd to_post_math_problem
####

grep -v "^#" dd.dat | awk '{print $1,$3,$4,$5,$6,$7,$8,$10,$11}'> dd_2.dat # Remove P1 an Cv columns

# Introduce blank line at the beginning:
sed '1i\\' dd_2.dat > dd_3.dat

# substitute 1st (blank line) by the header form Thermo*_SORTED_in_Ts.dat:
var="# VOLUME              P                  EL                E0              ET           ENTROPY           TS               G            Temperatures"
#var="oror"

# Because the $var contains "/", it does not work, e.g., if you try var="oror" We can use an alternate regex delimiter (~) as sed allows you to use any delimiter
sed "1s~.*~$var~" dd_3.dat > dd_4.dat 

mv dd_4.dat  solid_1__blocked_as_T_wise.dat

#####

# Remove blank lines:
sed '/^\s*$/d' solid_1__blocked_as_T_wise.dat >  bridge.dat

grep -v "^#" bridge.dat | sort -k2 -n   > kk.dat # sort as P-wise

#grep -v "^#" kk.dat | awk '{print $1,$3,$4,$5,$6,$7,$8,$10,$11}'> kk_2.dat # Remove P1 an Cv columns

sed '1i\\' kk.dat > kk_2.dat


# substitute 1st (blank line) by the header form Thermo*_SORTED_in_Ts.dat:
var="# VOLUME              P                  EL                E0              ET           ENTROPY           TS               G            Temperatures"
#var="oror"

# Because the $var contains "/", it does not work, e.g., if you try var="oror" We can use an alternate regex delimiter (~) as sed allows you to use any delimiter
sed "1s~.*~$var~" kk_2.dat > kk_3.dat 

mv kk_3.dat solid_1__sorted_as_P_wise.dat
