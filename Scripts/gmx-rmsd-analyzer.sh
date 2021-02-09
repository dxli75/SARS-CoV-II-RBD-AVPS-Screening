#!/bin/bash
echo "================================================================================"
if [ $# != 1 ]; then
  echo "This script is used to do a complete mmpbsa."
  echo "Usage: $0 [job_name]"
  echo "Author: Daixi Li"
  exit 1
else
  echo "The current usage: gmx-ncovs-run.sh $@"
fi
echo "################################################################################"
echo "STARTED AT: " `date`

echo ""
echo "PLEASE NAME THE CURRENT PROJECT!"
model=$1
echo "Make index file:"
#echo -e "ri 1-${rnum1}\nname 19 strimer\nri ${rnum3}-${rnum4}\nname 20 peptide\n19 & 3\nname 21 host-ca\n20 & 3\nname 22 guest-ca\nq" | gmx make_ndx -f ${model}-ion.gro -o index.ndx
echo -e "n18 & 4\nname 21 hostbb\n19 & 4\nname 22 guestbb\nq" | gmx make_ndx -f ${model}-ion.gro -n index.ndx -o index.ndx
host=18 && guest=19 && hostbb=21 && guestbb=22
charge=`grep "qtot" topol.top | tail -1 | awk '{ print $11}'`
echo "net charge:" $charge
#if [ $charge != 0 ]; then
#   #host=18 && guest=19 && fix=20 # For host-protein only
#   host=19 && guest=20 && fix=21  # For host-guest complex
#else
#   #host=15 && guest=16 && fix=17 # If No NaCl or KCl salts
#   #host=18 && guest=19 && fix=20 # For host-protein only and NaCl_Conc=0.1538 mol/L
#   host=19 && guest=20 && fix=21  # For host-guest complex and NaCl_Conc=0.1538 mol/L
#fi

echo "STEP5: Make preparation for RMSD"
dir=mmpbsa
mkdir ${dir}
echo 1 | gmx trjconv -s $model-eqmd01.tpr  -f $model-eqmd01.trr  -o ${dir}/${model}1.trr -dt 10
echo 1 | gmx trjconv -s $model-eqmd02.tpr  -f $model-eqmd02.trr  -o ${dir}/${model}2.trr -dt 10
echo 1 | gmx trjconv -s $model-eqmd03.tpr  -f $model-eqmd03.trr  -o ${dir}/${model}3.trr -dt 10
echo 1 | gmx trjconv -s $model-eqmd04.tpr  -f $model-eqmd04.trr  -o ${dir}/${model}4.trr -dt 10
echo 1 | gmx trjconv -s $model-eqmd05.tpr  -f $model-eqmd05.trr  -o ${dir}/${model}5.trr -dt 10 

cd ${dir}
echo -e "c\nc\nc\nc\nc\n" | gmx trjcat -f ${model}1.trr ${model}2.trr ${model}3.trr ${model}4.trr ${model}5.trr -o ${model}-m1.trr -cat -settime
echo -e "1\n1\n" | gmx trjconv -s ../$model-emsd.tpr -f $model-m1.trr -n ../index.ndx -o $model-m2.trr -pbc nojump -center 
echo -e "1\n1\n" | gmx trjconv -s ../$model-emsd.tpr -f $model-m2.trr -n ../index.ndx -o $model.trr -fit rot+trans

echo ""
echo "Check Results(rms) for Equilibration"
echo 4        4            | gmx rms -f $model.trr -s ../$model-emsd.tpr -n ../index.ndx -o $model-complex-rmsd.xvg
echo ${hostbb}  ${hostbb}  | gmx rms -f $model.trr -s ../$model-emsd.tpr -n ../index.ndx -o $model-host-rmsd.xvg
echo ${guestbb} ${guestbb} | gmx rms -f $model.trr -s ../$model-emsd.tpr -n ../index.ndx -o $model-guest-rmsd.xvg

rm $model-m?.trr  ${model}?.trr

echo "Well done!"