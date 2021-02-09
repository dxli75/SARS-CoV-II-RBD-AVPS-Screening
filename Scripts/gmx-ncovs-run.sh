#!/bin/bash
echo "================================================================================"
if [ $# != 4 ]; then
  echo "This script is used to do a complete molecular simulation,"
  echo "which includes energy minimization, pre-equilibration and "
  echo "equilibration or production simulation."
  echo "Usage: $0 [job_name] [total_residue_number] [gpu_id] [ntomp]"
  echo "Author: Daixi Li"
  exit 1
else
  echo "The current usage: gmx-ncovs-run.sh $@"
fi
echo "################################################################################"
echo "STARTED AT: " `date`
echo "Check GPU running:"
nvidia-smi

echo "STEP1: Check Files And Make Model"
echo ""
echo "PLEASE NAME THE CURRENT PROJECT!"
model=$1
MDPPATH=/home/dxli/simulation/projects/mdp
echo "Copy mdp file:"
cp ${MDPPATH}/emsd.mdp  ./
cp ${MDPPATH}/prmd.mdp  ./
cp ${MDPPATH}/prnvt.mdp ./
cp ${MDPPATH}/eqmd.mdp  ./
cp ${MDPPATH}/tip3p.gro ./

echo "Convert pdb to gmx:"
echo 6 | gmx pdb2gmx -f $model.pdb -o $model.gro -chainsep id -merge all -water tip3p -ignh > gmx.log

echo "Check topol file:"
if [ -e topol.top ]; then
   echo "topol.top file is available."
   sed -i "/^Protein_chain_A/s/Protein_chain_A/Protein/" topol.top
   sed -i "/^Protein_chain_B/s/Protein_chain_B/Protein/" topol.top
   sed -i "/^Protein_chain_C/s/Protein_chain_C/Protein/" topol.top
   sed -i "/^Protein_chain_D/s/Protein_chain_D/Protein/" topol.top
   sed -i "/^Protein_chain_E/s/Protein_chain_E/Protein/" topol.top
else
   echo "Error: topol.top file not found ."
   exit 1
fi

echo "Get system size:"
echo 1 | gmx editconf -f $model.gro -o $model.gro -princ > gmx.log
maxa=`grep "new system size :" gmx.log | awk '{ print $5+2.4}'`
maxb=`grep "new system size :" gmx.log | awk '{ print $6+2.4}'`
maxc=`grep "new system size :" gmx.log | awk '{ print $7+2.4}'`
echo "system size:" $maxa $maxb $maxc

echo "Set PBC size:"
#maxv=`echo -e "$maxa\n$maxb\n$maxc" | awk 'BEGIN {max = 0} {if ($1+0>max+0) max=$1 fi} END {print max}'`
gmx editconf -f $model.gro -o $model.gro -box $maxa $maxb $maxc -c 

echo "Solvate the system:"
gmx solvate -cp $model.gro -cs tip3p.gro -p topol.top -o $model-sol.gro

echo "Genion the system:"
gmx grompp -f emsd.mdp -c $model-sol.gro -p topol.top -o $model-em4ion.tpr -maxwarn 1

charge=`grep "qtot" topol.top | tail -1 | awk '{ print $11}'`
echo "net charge:" $charge
if [ $charge -ge 0 ]; then
   echo 13 | gmx genion -s $model-em4ion.tpr -p topol.top -o $model-ion.gro -nn $charge
else
   echo 13 | gmx genion -s $model-em4ion.tpr -p topol.top -o $model-ion.gro -np $((0-$charge))
fi

echo "Make index file:"
#echo -e "ri 1-597\nname 18 ace2\nri 598-$2\nname 19 peptide\n18&3\nname 20 fix\nq" | gmx make_ndx -f $model-ion.gro -o index.ndx
echo -e "ri 1-195\nname 18 ncovs\nri 196-$2\nname 19 peptide\n18&3\nname 20 fix\nq" | gmx make_ndx -f $model-ion.gro -o index.ndx
echo "Produce PDB file for check:"
echo 1 | gmx editconf -f $model-ion.gro -o $model-ion.pdb -n index.ndx

echo "Make posre file:"
echo 5 | gmx genrestr -f $model-ion.gro -n index.ndx -o posre-mc.itp
echo 4 | gmx genrestr -f $model-ion.gro -n index.ndx -o posre-bb.itp
echo 3 | gmx genrestr -f $model-ion.gro -n index.ndx -o posre-ca.itp
echo 20| gmx genrestr -f $model-ion.gro -n index.ndx -o posre-fx.itp

echo "Make top file:"
cat topol.top | sed "s/posre.itp/posre-mc.itp/" > topol-mc.top
cat topol.top | sed "s/posre.itp/posre-bb.itp/" > topol-bb.top
cat topol.top | sed "s/posre.itp/posre-ca.itp/" > topol-ca.top
cat topol.top | sed "s/posre.itp/posre-fx.itp/" > topol-fx.top

echo ""
echo "STEP2: Energy Minimization"
gmx grompp -f emsd.mdp -c ${model}-ion.gro -p topol.top -o $model-emsd.tpr -maxwarn 1
gmx mdrun -deffnm $model-emsd -gpu_id $3 -ntmpi 1 -ntomp $4 -v

echo ""
echo "STEP3: Pre-equilibration(each for 1 ns)"
gmx grompp -f prnvt.mdp -c $model-emsd.gro -r $model-emsd.gro -p topol.top -o $model-prmd01.tpr -maxwarn 1
gmx mdrun -deffnm $model-prmd01 -gpu_id $3 -ntmpi 1 -ntomp $4 -v

gmx grompp -f prmd.mdp -c $model-prmd01.gro -r $model-prmd01.gro -t $model-prmd01.trr \
                       -e $model-prmd01.edr -p topol-mc.top     -o $model-prmd02.tpr -maxwarn 1
gmx mdrun -deffnm $model-prmd02 -gpu_id $3 -ntmpi 1 -ntomp $4 -v

gmx grompp -f prmd.mdp -c $model-prmd02.gro -r $model-prmd02.gro -t $model-prmd02.trr \
                       -e $model-prmd02.edr -p topol-bb.top     -o $model-prmd03.tpr -maxwarn 1
gmx mdrun -deffnm $model-prmd03 -gpu_id $3 -ntmpi 1 -ntomp $4 -v

gmx grompp -f prmd.mdp -c $model-prmd03.gro -r $model-prmd03.gro -t $model-prmd03.trr \
                       -e $model-prmd03.edr -p topol-ca.top     -o $model-prmd04.tpr -maxwarn 1
gmx mdrun -deffnm $model-prmd04 -gpu_id $3 -ntmpi 1 -ntomp $4 -v

gmx grompp -f prmd.mdp -c $model-prmd04.gro -r $model-prmd04.gro -t $model-prmd04.trr \
                       -e $model-prmd04.edr -p topol-fx.top     -o $model-prmd05.tpr -maxwarn 1
gmx mdrun -deffnm $model-prmd05 -gpu_id $3 -ntmpi 1 -ntomp $4 -v

echo ""
echo "STEP4: Equilibration(20 ns)"
gmx grompp -f eqmd.mdp  -c $model-prmd05.gro -r $model-prmd05.gro -p topol.top -o $model-eqmd01.tpr -maxwarn 2
gmx mdrun -deffnm $model-eqmd01 -gpu_id $3 -ntmpi 1 -ntomp $4 -v
gmx grompp -f eqmd.mdp  -c $model-eqmd01.gro -r $model-eqmd01.gro -p topol.top -o $model-eqmd02.tpr -maxwarn 2
gmx mdrun -deffnm $model-eqmd02 -gpu_id $3 -ntmpi 1 -ntomp $4 -v
gmx grompp -f eqmd.mdp  -c $model-eqmd02.gro -r $model-eqmd02.gro -p topol.top -o $model-eqmd03.tpr -maxwarn 2
gmx mdrun -deffnm $model-eqmd03 -gpu_id $3 -ntmpi 1 -ntomp $4 -v
gmx grompp -f eqmd.mdp  -c $model-eqmd03.gro -r $model-eqmd03.gro -p topol.top -o $model-eqmd04.tpr -maxwarn 2
gmx mdrun -deffnm $model-eqmd04 -gpu_id $3 -ntmpi 1 -ntomp $4 -v
gmx grompp -f eqmd.mdp  -c $model-eqmd04.gro -r $model-eqmd04.gro -p topol.top -o $model-eqmd05.tpr -maxwarn 2
gmx mdrun -deffnm $model-eqmd05 -gpu_id $3 -ntmpi 1 -ntomp $4 -v

echo ""
echo "Check Results(pdb) for Equilibration"
echo 1 | gmx editconf -f $model-eqmd01.gro -n index.ndx -o $model-eqmd01.pdb
echo 1 | gmx editconf -f $model-eqmd02.gro -n index.ndx -o $model-eqmd02.pdb
echo 1 | gmx editconf -f $model-eqmd03.gro -n index.ndx -o $model-eqmd03.pdb
echo 1 | gmx editconf -f $model-eqmd04.gro -n index.ndx -o $model-eqmd04.pdb
echo 1 | gmx editconf -f $model-eqmd05.gro -n index.ndx -o $model-eqmd05.pdb

echo "Well Done!"