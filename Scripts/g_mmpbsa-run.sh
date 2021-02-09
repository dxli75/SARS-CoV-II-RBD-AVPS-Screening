#!/bin/bash
echo "================================================================================"
if [ $# != 1 ]; then
  echo "This script is used to do a complete mm-pbsa."
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
charge=`grep "qtot" topol.top | tail -1 | awk '{ print $11}'`
echo "net charge:" $charge
if [ $charge != 0 ]; then
   host=18 && guest=19
else
   host=15 && guest=16
fi

echo "STEP5: Make preparation for mmpbsa"
dir=mmpbsa
mkdir ${dir}
gmx grompp  -f emsd.mdp -c ${model}-eqmd05.gro -p topol.top -o ${dir}/${model}.tpr -maxwarn 2
echo 1 | gmx trjconv -s $model-eqmd05.tpr  -f $model-eqmd05.trr  -o ${dir}/${model}-m1.trr -b 15000 -e 20000 -dt 10

cd ${dir}
#echo -e "c\nc\n" | gmx trjcat -f ${model}a.trr ${model}b.trr -o ${model}-c.trr -cat -settime
echo -e "1\n1\n" | gmx trjconv -s ../${model}-emsd.tpr -f ${model}-m1.trr -n index.ndx -o ${model}-m2.trr -pbc nojump -center 
echo -e "1\n1\n" | gmx trjconv -s ../${model}-emsd.tpr -f ${model}-m2.trr -n index.ndx -o ${model}.trr -fit rot+trans

GMXPATH=/home/dxli/simulation/gromacs2020.2/bin/g_mmpbsa
cp ../index.ndx           ./
cp ${GMXPATH}/g_mmpbsa    ./
cp ${GMXPATH}/*.mdp       ./
cp ${GMXPATH}/*.py        ./
cp ${GMXPATH}/energy2bfac ./
chmod 755 ./g_mmpbsa ./energy2bfac

echo "STEP5.1: Calculate MM binding energy in vacuum+coulomb + vdw"
echo $host $guest | ./g_mmpbsa -s $model.tpr -f $model.trr -n index.ndx -mme -pdie 2 -decomp

echo "STEP5.2: Calculate PB apolar binding energy in implicit solution"
echo $host $guest | ./g_mmpbsa -s $model.tpr -f $model.trr -n index.ndx -i polar.mdp -nomme -pbsa -decomp

echo "STEP5.3: Calculate SA apolar binding energy in implicit solution"
echo "Use sav model"
echo $host $guest | ./g_mmpbsa -s $model.tpr -f $model.trr -n index.ndx -i apolar_sav.mdp -nomme -pbsa -decomp -apol sav.xvg -apcon sav_contrib.dat

echo "STEP5.4: Calculate average binding energy"
echo "Use the sav model"
python3 ./MmPbSaStat.py -bs -nbs 1000 -m energy_MM.xvg -p polar.xvg -a sav.xvg -os summary_energy_sav.dat

echo "STEP5.5: Calculate final contribute energy"
# Note: the order of input, -m MM -p polar -a apolar must be kept.
# Or you can write the paths of the three files into metafile.dat.
echo "Use the sav model"
python3 ./MmPbSaDecomp.py -m contrib_MM.dat -p contrib_pol.dat -a sav_contrib.dat -bs -nbs 1000 -ct 999 -o final_contrib_energy_sav.dat -om energymaping_sav.dat

echo "STEP5.6: Convert contribute energy into b-factor and write into pdb file for pymol"
echo $host $guest | ./energy2bfac -s $model.tpr -i energymaping_sav.dat -n index.ndx -c ${model}-sav0.pdb -s1 ${model}-sav1.pdb -s2 ${model}-sav2.pdb

