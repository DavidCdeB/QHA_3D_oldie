#
# Preparing data for surface analysis:

# Remove if they exist:
for i in dd.dat \
         dd_2.dat\
         dd_3.dat\
         dd_4.dat\
         solid_1__xyz_blocked_as_T_wise.dat\
         bridge.dat\
         kk.dat\
         kk_2.dat\
         kk_3.dat\
         solid_1__xyz_sorted_as_P_wise.dat
do 

rm -Rf $i

done

# Copy to be safe:
cp Thermo_functions_with_P_and_G_sorted_with_EL_from_EOS__VOLUME_and_All_energies_divided_by_F_unit_SORTED_in_Ts_blocked_with_header.dat dd.dat

# Only maintain P, G, T columns
grep -v "^#" dd.dat | awk '{print $3,$10,$11}'> dd_2.dat

# Introduce blank line at the beginning:
sed '1i\\' dd_2.dat > dd_3.dat

# Substitute 1st (blank line) by a "P G T " header:
var="# P (Gpa)       G / F. unit  (a.u.)     T (K)"

# Because the $var contains "/", it does not work, e.g., if you try var="oror" We can use an alternate regex delimiter (~) as sed allows you to use any delimiter
sed "1s~.*~$var~" dd_3.dat > dd_4.dat

mv dd_4.dat  solid_1__xyz_blocked_as_T_wise.dat

# Remove blank lines:
sed '/^\s*$/d' solid_1__xyz_blocked_as_T_wise.dat >  bridge.dat

# Sort as P-wise:
grep -v "^#" bridge.dat | sort -k1 -n   > kk.dat

# Introduce blank line at the beginning:
sed '1i\\' kk.dat > kk_2.dat

# Substitute 1st (blank line) by a "P G T " header:
var="# P (Gpa)       G / F. unit  (a.u.)     T (K)"

# Because the $var contains "/", it does not work, e.g., if you try var="oror" We can use an alternate regex delimiter (~) as sed allows you to use any delimiter
sed "1s~.*~$var~" kk_2.dat > kk_3.dat

mv kk_3.dat solid_1__xyz_sorted_as_P_wise.dat


