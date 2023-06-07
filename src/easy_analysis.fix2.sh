#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate md_env

#------------------
lig="atp"
workdir="./"
gmx="/usr/local/gromacs/bin/gmx"
#------------------
cd $workdir
mkdir -p analysis_result

#(echo " 1 | 13"; echo 'q') | $gmx make_ndx -f md_0_50.gro -n index.ndx -o index_analysis.ndx
echo "0 0" | $gmx trjconv -s md_0_50.tpr -f md_0_50.xtc -o md_0_50_nojump.xtc -center -pbc nojump -ur compact -n index.ndx
echo "21 0" | $gmx trjconv -s md_0_50.tpr -f md_0_50_nojump.xtc -o md_0_50_fit.xtc -fit rot+trans -dt 100 -n index.ndx

echo "11 0" | $gmx energy -f em.edr -o ./analysis_result/potential.xvg
echo "16 0" | $gmx energy -f nvt.edr -o ./analysis_result/temperature.xvg
echo "17 0" | $gmx energy -f npt.edr -o ./analysis_result/pressure.xvg
echo "23 0" | $gmx energy -f npt.edr -o ./analysis_result/density.xvg
echo "9 0"   | $gmx energy -f md_0_50.edr -o ./analysis_result/Coul_SR.xvg
echo "8 0"   | $gmx energy -f md_0_50.edr -o ./analysis_result/LJ_SR.xvg

# ----------------- Backbone/蛋白的分子结构变化的程度-----------------------------------------------------

echo "4 4" | $gmx rms -s em.tpr -f md_0_50_fit.xtc -o ./analysis_result/rmsd_xtal.xvg -tu ns -n index.ndx 
echo "4 4" | $gmx rms -s md_0_50.tpr -f md_0_50_fit.xtc -o ./analysis_result/rmsd_protein.xvg -tu ns -n index.ndx
echo "3 3" | $gmx rms -s md_0_50.tpr -f md_0_50_fit.xtc -o ./analysis_result/rmsd_complex.xvg -tu ns -n index.ndx 
echo "13 13" | $gmx rms -s md_0_50.tpr -f md_0_50_fit.xtc -o ./analysis_result/rmsd_ligand.xvg -tu ns -n index.ndx

## ---------------每个原子的涨落，计算均方根涨落(RMSF)------------------------------

echo "4" | $gmx rmsf -s em.tpr -f md_0_50_fit.xtc -o ./analysis_result/rmsf_xtal.xvg -n index.ndx
echo "1" | $gmx rmsf -s md_0_50.tpr -f md_0_50_fit.xtc -o ./analysis_result/rmsf_protein.xvg -res -n index.ndx
echo "3" | $gmx rmsf -s md_0_50.tpr -f md_0_50_fit.xtc -o ./analysis_result/rmsf_complex.xvg -n index.ndx

## ---------------- 溶剂可及性表面积，并计算平均值和方差值--------------------------

echo "1" | $gmx sasa -s md_0_50.tpr -f md_0_50_fit.xtc -o ./analysis_result/sasa.xvg -oa ./analysis_result/oa.xvg -or ./analysis_result/resasa.xvg -n index.ndx
$gmx sasa -s md_0_50.tpr -f md_0_50_fit.xtc -o ./analysis_result/sasa_Hydrophobe.xvg -surface 'group protein' -output '"Hydrophobic" group protein and charge {-0.2 to 0.2};"Hydrophilic" group protein and not charge {-0.2 to 0.2}'
$gmx analyze -f ./analysis_result/resasa.xvg

## ---------------- 回旋半径，并计算平均值和方差值---------------

echo "1" | $gmx gyrate -s md_0_50.tpr -f md_0_50_fit.xtc -o ./analysis_result/gyrate.xvg -n index.ndx
$gmx analyze -f ./analysis_result/gyrate.xvg

## ----------------- 氢键 ，并计算平均值和方差值---------------------

echo "1 13"| $gmx hbond -s md_0_50.tpr -f md_0_50_fit.xtc -num ./analysis_result/hbnum.xvg -n index.ndx
$gmx analyze -f ./analysis_result/hbnum.xvg

## ----------------- 蛋白质二级结构分析 -----------------------

conda activate dssp_env
echo "1" | $gmx do_dssp -f md_0_50_fit.xtc -s md_0_50.tpr -n index.ndx -o ./analysis_result/secondary-structure.xpm -sc secondary-structure.xvg -ssdump dump.dat -map ss.map -ver 4
conda activate md_env
$gmx xpm2ps -f ./analysis_result/secondary-structure.xpm -di md.m2p -o ./analysis_result/ss.eps2
ps2pdf -sPAPERSIZE=ledger ./analysis_result/ss.eps2 ./analysis_result/ss.pdf

## ------------------ 对比前后结构变化情况 -----------------------

echo "1 1" | $gmx confrms -f1 solv_ions.gro -f2 md_0_50.gro -o ./analysis_result/fit.pdb

## ------------------ 电影形成 --------------------------

echo "1 21" | $gmx trjconv -s md_0_50.tpr -f md_0_50_fit.xtc -center -n index.ndx -o ./analysis_result/complex.pdb -tu ps -dt 100

## ------------------ 基于协方差矩阵分析特征向量即Free energy landscape (FEL)的PCA分析(20~50ns之间的结果) ------------------- 

echo "1 1" | $gmx trjconv -f md_0_50_nojump.xtc -b 20000 -s md_0_50.tpr -fit rot+trans -o pcamdfit.xtc
echo "3 3" | $gmx covar -s md_0_50.gro -f pcamdfit -o eigenvalues.xvg -v eigenvectors.trr -xpma covapic.xpm
echo "3 3" | $gmx anaeig -f pcamdfit -s md_0_50.gro -v eigenvectors.trr -first 1 -last 2 -2d 2d.xvg
echo "3 3" | $gmx anaeig -f pcamdfit -s md_0_50.gro -v eigenvectors.trr -last 1 -proj pc1.xvg
echo "3 3" | $gmx anaeig -f pcamdfit -s md_0_50.gro -v eigenvectors.trr -first 2 -last 2 -proj pc2.xvg
perl sham.pl -i1 pc1.xvg -i2 pc2.xvg -data 1 -data2 1 -o gsham_pcainput.xvg
$gmx sham -f gsham_pcainput.xvg -ls FEL_shampca.xpm -tsham 300
python2 xpm2txt.py -f FEL_shampca.xpm -o ./analysis_result/free-energy-landscapePCA.txt
#rm -f gsham_pcainput.xvg

## ------------------ 分析最后20~50ns的Free energy landscape (FEL) RMSD和Rg ------------------- 

echo "1" | $gmx gyrate -s md_0_50.tpr -f pcamdfit.xtc -o FEL_gyrate.xvg
echo "1 1" | $gmx rms -s  md_0_50.tpr -f pcamdfit.xtc -o FEL_rmsd.xvg
perl sham.pl -i1 FEL_rmsd.xvg -i2 FEL_gyrate.xvg -data 1 -data2 1 -o gsham_RMSDRginput.xvg
$gmx sham -f gsham_RMSDRginput.xvg -ls FEL_shamRMSDRg.xpm -tsham 300
python2 xpm2txt.py -f FEL_shamRMSDRg.xpm -o ./analysis_result/free-energy-landscapeRMSDRg.txt
#rm -f gsham_RMSDRginput.xvg

## ------------------ 从最后20ns获取典型结构用于药效团模拟 -------------------

echo "21 0" | $gmx trjconv -s md_0_50.tpr -f md_0_50_nojump.xtc -b 30000 -o md_30_50_fit.xtc -fit rot+trans -dt 100 -n index.ndx
echo "21 21" | $gmx cluster -f md_30_50_fit.xtc -s md_0_50.tpr -method gromos -n index.ndx -o rmsd-clust.xpm -g cluster.log -dist rmsd-dist.xvg -cutoff 0.3 -ev -sz -tr -ntr -clid -clndx -cl -method gromos 
sed -n '/MODEL        1/,/TER/p' clusters.pdb | head -n -1 | tail -n+2 > ref-complex.pdb

## ------------------ 计算结合能Method1[gmx_mmpbsa] -----------------------

# 提取最后稳定构象进行计算，不然时间无限长。提取最后1ns。总长10ns

mkdir -p binding_analysis_1
cp ./md_0_50.xtc ./binding_analysis_1/md_0_50.xtc 
cp ./md_0_50.tpr ./binding_analysis_1/ie.tpr
cp ./index.ndx ./binding_analysis_1/index.ndx
sed /BEC/s//${lig^^}/ gmx_mmpbsa.bsh -i 
cp ./gmx_mmpbsa.bsh ./binding_analysis_1/gmx_mmpbsa.bsh
cp ./plotmmpb.py ./binding_analysis_1/plotmmpb.py
cd binding_analysis_1
./gmx_mmpbsa.bsh
cp _pid~MMPBSA.dat _pid~MMPBSA.dat.backup
sed "s/([^)]*)//g" -i _pid~MMPBSA.dat
sed -n -e :a -e '1,5!{P;N;D;};N;ba' -i _pid~MMPBSA.dat
python3 plotmmpb.py
cd ..

mkdir analysis
cp *.xvg ./analysis
cp md_0_10.log ./analysis
cp *.pdb ./analysis
#cp resp.log ./analysis
#cp analysis.log ./analysis
#cp rmsd-clust.xpm cluster.pdb ./analysis

## ------------------ 计算结合能，关键在于polar.mdp和apolar_sasa.mdp文件设置 -----------------------

# echo 1 0 | gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_9_10_noPBC.xtc -n index.ndx -pbc mol -ur compact -b 9000 -e 10000 -dt 50
# mkdir 1EBZ2
#cp md_0_10_fit.xtc ./1EBZ2/npt.xtc
#cp md_0_10.tpr ./1EBZ2/npt.tpr
#cp protein.top ./1EBZ2/protein.itp
#cp jz4.top ./1EBZ2/ligand.itp
#cp topol.top ./1EBZ2/
#(echo "del 20-22"; echo 'q') | gmx make_ndx -f md_0_10.gro -n index.ndx -o index_binding.ndx
#cp index_binding.ndx ./1EBZ2/index.ndx
#cp INPUT.dat ./1EBZ2/INPUT.dat
#cp -r charmm36-feb2021.ff ./1EBZ2/
# $GMXPBSAHOME/gmxpbsa0.sh
# $GMXPBSAHOME/gmxpbsa1.sh
# $GMXPBSAHOME/gmxpbsa2.sh

#   ###Three steps calculation
#   ####(a) Calculation of potential energy in Vacuum
#   echo "1 13" | g_mmpbsa -f md_0_10_fit.xtc -s md_0_10.tpr -n index.ndx -pdie 2 -decomp -mm PAvM_energy_mm.xvg -mmcon PAvM_contrib_mm.dat -dt 10
#   ####(b) Calculation of polar solvation energy
#   echo "1 13" | g_mmpbsa -f md_0_10_fit.xtc -s md_0_10.tpr -n index.ndx -i polar.mdp -nomme -pbsa -decomp -pol PAvM_polar.xvg -pcon PAvM_contrib_pol.dat -dt 10
#   ####(c) Calculation of non-polar solvation energy --> sasa model
#   echo "1 13" | g_mmpbsa -f md_0_10_fit.xtc -s md_0_10.tpr -n index.ndx -i apolar_sasa.mdp -nomme -pbsa -decomp -apol PAvM_wca.xvg -apcon  PAvM_contrib_wca.dat -dt 10
#   python3 MmPbSaStat.py -m PAvM_energy_mm.xvg  -p PAvM_polar.xvg -a PAvM_wca.xvg -bs -nbs 2000 -of PAvM_full_energy.dat -os PAvM_summary_energy.dat
#   python3 MmPbSaDecomp.py -bs -nbs 2000 -m PAvM_contrib_mm.dat -p PAvM_contrib_pol.dat -a PAvM_contrib_wca.dat -o PAvM_final_contrib_energy.dat -om PAvM_energyMapIn.dat

###One step calculation
#echo "1 13" | g_mmpbsa -f md_0_1.xtc -s md_0_10.tpr -n index.ndx -i pbsa.mdp -pdie 2 -pbsa -decomp
###Average Binding Energy Calculation
#echo " " | sudo -S python MmPbSaStat.py -m energy_MM.xvg -p polar.xvg -a apolar.xvg

## ------------------ 计算结合能Method2[gmx_MMPBSA] -----------------------

mkdir -p binding_analysis_2
cp /share/home/dw_user0001/auto_md/scr/*.in binding_analysis_2/
cp npt.tpr index.ndx md_30_50_fit.xtc md_0_50.tpr atp.itp topol.top binding_analysis_2/
cp -r charmm36-feb2021.ff binding_analysis_2/
cp atp* binding_analysis_2/

cd binding_analysis_2
#gmx_MMPBSA -O -i mmgbsa.in -cs npt.tpr -ci index.ndx -cg 1 13 -ct mdfit80_100ns.xtc  -cp topol.top -o FINAL_RESULTS_MMGBSA.dat -ls atp_ini.pdb -eo FINAL_RESULTS_MMGBSA.csv -nogui
#mv FINAL_RESULTS* ../$project_name/mm_gbsa
gmx_MMPBSA -O -i mmpbsa.in -cs md_0_50.tpr -ci index.ndx -cg 1 13 -ct md_30_50_fit.xtc -cp topol.top -rg 1 -lg 13 -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -nogui
mv FINAL_RESULTS* ../$project_name/mm_pbsa
