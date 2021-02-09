#!/bin/csh
echo "This script is used to merge 10 PDB files into one."
echo ""
date
set model=$1
echo "REMARK   Merge 10 PDB files into one" > complex-all.pdb

set i = 1
while ( $i <= 10 )
  echo "MODEL   " $i >> complex-all.pdb
  cat ${model}$i.pdb >> complex-all.pdb
  echo "ENDMDL"      >> complex-all.pdb
  @ i++
end

echo "END"   >> complex-all.pdb