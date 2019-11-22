#!/bin/bash
#read var_year
#echo "The year is: $var_year"

echo -n "Enter  the summary of activities make before save "
read var_activity
echo "Your name is: $var_activity"
#echo -n "Press [ENTER] to continue,...: "
echo "You can go on!...."
git add *.f90
git add ma*
git add *.sh
git add *.py
git add dados.in
git add en*
git commit -m "$var_activity"
git push origin 
