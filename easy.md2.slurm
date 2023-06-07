#!/bin/bash

# file: easy_md.slurm

#SBATCH --job-name=gromaces_md
#SBATCH --partition=gpu
#SBATCH --nodelist=c003
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --gres=gpu:1
#SBATCH --mem=400G
#SBATCH --output=gromacs_md_%j.log
#SBATCH --error=gromacs_md_%j.err
#SBATCH --distribution=cyclic:block

export OMP_NUM_THREADS=36

echo TIME subset start $(date)
source /share/apps/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate md_env

#------------------
box="2.0"
ion1="NA"
ion2="CL"
gpu="1"
lig="atp"
workdir="[pwd]/auto_md/TgaRuvbl1_ATP_md"
SRC="[pwd]/auto_md/src"
GMX="/usr/local/gromacs/bin/gmx"

# SET NUMBER OF MPI TASKS 
#mpi="1"
mpi=$(($SLURM_NTASKS_PER_NODE*$SLURM_JOB_NUM_NODES))
# SET NUMBER OF OPENMP THREADS
#omp="36"
omp=$SLURM_CPUS_PER_TASK

mkdir -p $workdir 
cd $workdir
cp $SRC/50ns_mdp/* ./
cp -r $SRC/sort_mol2_bonds.pl ./
cp -r $SRC/complex_topol.py ./

sed /JZ4/s//${lig^^}/ minim.mdp -i
sed /JZ4/s//${lig^^}/ nvt.mdp -i
sed /JZ4/s//${lig^^}/ npt.mdp -i
sed /JZ4/s//${lig^^}/ md.mdp -i 

echo "--------- Prepared Complex.gro Files---------------"

obabel -ismi ${lig}.smi -opdb --conformer --gen3d --weighted -h --ff GAFF --partialcharge eem -O ${lig}_h.pdb
grep -E 'HETATM |ATOM ' ${lig}_h.pdb > ${lig}_new.pdb
sed -i "s/UNL/${lig^^}/g" ${lig}_new.pdb
echo TER >> ${lig}_new.pdb

antechamber -i ${lig}_new.pdb -fi pdb -o ${lig}.gjf -fo gcrt -pf y -gn "%nproc=30" -gm "%mem=240MB" -ch ${lig}.chk -gk "#HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) iop(6/50=1) opt" -ge ${lig}.gesp -gv 1
g16 < ${lig}.gjf > ${lig}.out
antechamber -i ${lig}.out -fi gout -c resp -o ${lig}_antechamber.mol2 -fo mol2 -pf y  -rn ${lig^^}

#antechamber -i ${lig}_new.pdb -fi pdb -o ${lig}_antechamber.mol2 -fo mol2 -c bcc -rn ${lig^^} -s 2 -at gaff2 -pf y -m 1

parmchk2 -i ${lig}_antechamber.mol2 -f mol2 -o ${lig}_antechamber.frcmod

cat <<EOF > tl.in
source oldff/leaprc.ff99SB
source leaprc.gaff2
loadamberparams ${lig}_antechamber.frcmod
lig=loadmol2 ${lig}_antechamber.mol2
check lig
saveamberparm lig ${lig}.prmtop ${lig}.inpcrd
quit
EOF

tleap -f tl.in

acpype -p ${lig}.prmtop -x ${lig}.inpcrd -d
cp ./${lig^^}.amb2gmx/${lig^^}_GMX.gro ./${lig}.gro
cp ./${lig^^}.amb2gmx/${lig^^}_GMX.top ./${lig}.top
#mv ./${lig^^}.amb2gmx/posre_${lig^^}.itp ./posre_${lig}.itp
sed -n '/atomtypes/,/moleculetype/{/moleculetype/b;p}' ${lig}.top > ${lig}.prm
sed -n '/moleculetype/,/system/{/system/b;p}' ${lig}.top > ${lig}.itp
rm -rf ${lig^^}.amb2gmx
rm -rf QOUT punch tl.in
rm -rf ${lig}.out ${lig}.gesp ${lig}.chk ${lig}.gjf
rm ${lig}.prmtop ${lig}.inpcrd ${lig}_antechamber.mol2 ${lig}_antechamber.frcmod leap.log
 
pdbfixer receptor.pdb --keep-heterogens=none --add-atoms=heavy --replace-nonstandard --add-residues --output protein.pdb
$GMX pdb2gmx -f protein.pdb -o protein_processed.gro -ff amber99sb -water spc -ignh -nochargegrp -p topol.top -i posre_protein.itp
mv topol.top protein.top
sed -n '/moleculetype/,/restraint/{/restraint/b;p}' protein.top > protein.itp
sed -i 's/Protein_chain_A/Protein /' protein.itp
cp -r $SRC/topol.top ./
python3 complex_topol.py protein_processed.gro ${lig}.gro
sed -i "s/1  ${lig^^}/1${lig^^}  /g" complex.gro

echo "---------- Molecular Dynamic Processing ---------------"

(echo "0 & ! a H*"; echo "name 3 ${lig^^}-H"; echo "q") | $GMX make_ndx -f ${lig}.gro -o index_${lig}.ndx
echo "3" | $GMX genrestr -f ${lig}.gro -n index_${lig}.ndx -o posre_${lig}.itp -fc 1000 1000 1000
sed -i '/; Include water topology/i\; Include Ligand position restraints\n#ifdef POSRES_LIGAND\n#include "'"./posre_${lig}.itp"'"\n#endif\n' topol.top
 
$GMX editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0
$GMX solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
$GMX grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo 15 | $GMX genion -s ions.tpr -o solv_ions.gro -p topol.top -pname ${ion1} -nname ${ion2} -neutral -conc 0.15

(echo "1 | 13"; echo "q") | $GMX make_ndx -f solv_ions.gro -o index.ndx

$GMX grompp -f minim.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -o em.tpr
$GMX mdrun -s em.tpr -v -deffnm em -ntmpi $mpi -ntomp $omp

$GMX grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2
$GMX mdrun -s nvt.tpr -v -deffnm nvt -ntmpi $mpi -ntomp $omp

$GMX grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 2
$GMX mdrun -s npt.tpr -v -deffnm npt -ntmpi $mpi -ntomp $omp

$GMX grompp -f md.mdp -c npt.gro -t npt.cpt -r npt.gro -p topol.top -n index.ndx -o md_0_50.tpr
$GMX mdrun -deffnm md_0_50 -ntmpi $mpi -ntomp $omp -nb gpu -bonded gpu -pme gpu -s md_0_50.tpr
mkdir -p ok

echo TIME subset Finish $(date)
