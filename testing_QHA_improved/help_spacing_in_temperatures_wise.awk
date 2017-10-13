# 
# awk to add blank line every time Temperature (x value) changes. 
# This will end up with a blocking pattern of 11 pressures.

/^[[:blank:]]*#/ {next} # ignore comments (lines starting with #)  # This will not destroy the header from Thermo*.dat
NF < 3 {next} # ignore lines which don’t have at least 3 columns  
#$11 != prev {printf ~\n~; prev=$11} # print blank line whenever the value of column $11 "Temperature (K)" changes. #for some reason this does not work
#{print} # print the line 
$11 != prev {print ""; prev=$11}
{print}
