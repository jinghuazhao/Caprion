#!/usr/bin/bash

export TMPDIR=${HPC_WORK}/work
export dir=~/rda/projects/Caprion_proteomics
export caprion=${dir}/pilot

. /etc/profile.d/modules.sh 
module load ceuadmin/stata

function pilot_setup()
# RCN3_HUMAN	Q96D15	RCN3	Reticulocalbin-3 (EF-hand calcium-binding protein RLP49)	endoplasmic reticulum [GO:0005783]; endoplasmic reticulum lumen [GO:0005788]	calcium ion binding [GO:0005509]	collagen biosynthetic process [GO:0032964]; ERAD pathway [GO:0036503]; lung epithelium development [GO:0060428]; phospholipid homeostasis [GO:0055091]; positive regulation of peptidase activity [GO:0010952]; protein secretion [GO:0009306]; protein transport [GO:0015031]; regulation of protein kinase B signaling [GO:0051896]; surfactant homeostasis [GO:0043129]
# FCGRN_HUMAN	P55899	FCGRT	IgG receptor FcRn large subunit p51 (FcRn) (IgG Fc fragment receptor transporter alpha chain) (Neonatal Fc receptor)	external side of plasma membrane [GO:0009897]; extracellular space [GO:0005615]; integral component of membrane [GO:0016021]; plasma membrane [GO:0005886]	beta-2-microglobulin binding [GO:0030881]; IgG binding [GO:0019864]	IgG immunoglobulin transcytosis in epithelial cells mediated by FcRn immunoglobulin receptor [GO:0002416]; immune response [GO:0006955]
{
  sed '2,3d' ${caprion}/data/caprion.sample > ${caprion}/data/caprion.dat
  cd data
  for chr in {1..22}
  do
    ln -sf caprion-${chr}.bgen chr${chr}.bgen
#   cp caprion-${chr}.bgen.bgi chr${chr}.bgen.bgi
    ln -sf ${HPC_WORK}/data/interval/chr${chr}.bgen.bgi
  done
  ln -sf ${HPC_WORK}/data/interval/SNPinfo.dta.gz
  ln -sf ${HPC_WORK}/data/interval/Chunks_15.dta
  cut -d' ' -f1-3 caprion.sample > interval.sample
  sed '1,2d' interval.sample  > sample_info.txt
  stata <<\ \ END
    insheet id idorder missing using sample_info.txt, delim(" ")
    format id %15.0g
    format idorder %15.0g
    gzsave sample_info, replace
  END
  rm sample_info.txt
  cd -
}

function phase2_setup()
{
  sed '2,3d' ${caprion}/data2/phase2.sample > ${caprion}/data2/phase2.dat
  sed '2,3d' ${caprion}/data2/phase2_invn.sample > ${caprion}/data2/phase2_invn.dat
  cd data2
  for chr in {1..22}
  do
    ln -sf caprion-${chr}.bgen chr${chr}.bgen
#   cp caprion-${chr}.bgen.bgi chr${chr}.bgen.bgi
    ln -sf ${HPC_WORK}/data/interval/chr${chr}.bgen.bgi
  done
  ln -sf ${HPC_WORK}/data/interval/SNPinfo.dta.gz
  ln -sf ${HPC_WORK}/data/interval/Chunks_15.dta
  cut -d' ' -f1-3 phase2.sample > interval.sample
  sed '1,2d' interval.sample  > sample_info.txt
  stata <<\ \ END
    insheet id idorder missing using sample_info.txt, delim(" ")
    format id %15.0g
    format idorder %15.0g
    gzsave sample_info, replace
  END
  rm sample_info.txt
  cd -
}

function pilot_run()
{
  export uniprot=(Q96D15 P55899)
  export prot=(RCN3 FCGRT)
  for i in {0..1}
  do
    export y=${uniprot[$i]}
    export trait=${prot[$i]}
    if [ ! -d ${caprion}/${trait} ]; then mkdir ${caprion}/${trait}; fi
    echo ${y} -- ${trait}
    stata-mp -b do utils/gwas.do
  done
  for i in {0..1}
  do
    export y=${uniprot[$i]}
    export trait=${prot[$i]}
    echo ${y} -- ${trait}
    head -1 ${caprion}/${trait}/interval.*.All.txt
    grep -w rs113886122 ${caprion}/${trait}/interval.*.All.txt
  done
  for i in {0..1}
  do
    export y=${uniprot[$i]}_invn
    export trait=${prot[$i]}_invn
    if [ ! -d ${caprion}/${trait} ]; then mkdir ${caprion}/${trait}; fi
    echo ${y} -- ${trait}
    stata-mp -b do utils/gwas.do
  done
  for i in {0..1}
  do
    export y=${uniprot[$i]}_invn
    export trait=${prot[$i]}_invn
    echo ${y} -- ${trait}
    head -1 ${caprion}/${trait}/interval.*.All.txt
    grep -w rs113886122 ${caprion}/${trait}/interval.*.All.txt
  done
}

function phase2_run()
{
  export full=(RCN3_442625488_VADQDGDSMATR RCN3_442666668_EVAKEFDQLTPEESQAR RCN3_All RCN3_DR)
  export abbrev=(RCN3_44262548~R RCN3_44266666~R RCN3_All RCN3_DR)
  for i in {0..3}
  do
    export y=${full[$i]}
    export trait=${abbrev[$i]}
    if [ ! -d ${caprion}/${y} ]; then mkdir ${caprion}/${y}; fi
    echo ${y} -- ${trait}
    stata-mp -b do ${caprion}/utils/gwas2.do
  done
  for i in {0..3}
  do
    export y=${full[$i]}
    export trait=${abbrev[$i]}
    echo ${y} -- ${trait}
    head -1 ${caprion}/${y}/interval.*.All.txt
    grep -w rs113886122 ${caprion}/${y}/interval.*.All.txt
  done
# _invn
  export abbrev=(RCN3_44262548~n RCN3_44266666~R RCN3_All_invn RCN3_DR_invn)
  for i in {0..3} 
  do
    export y=${full[$i]}_invn
    export trait=${abbrev[$i]}
    if [ ! -d ${caprion}/${y} ]; then mkdir ${caprion}/${y}; fi
    echo ${y} -- ${trait}
    stata-mp -b do ${caprion}/utils/gwas2_invn.do
  done
  for i in {0..3}
  do
    export y=${full[$i]}_invn
    export trait=${abbrev[$i]}_invn
    echo ${y} -- ${trait}
    head -1 ${caprion}/${y}/interval.*.All.txt
    grep -w rs113886122 ${caprion}/${y}/interval.*.All.txt
  done
}

function RCN3_FCGRT_plink()
{
  gunzip -c bgen/Q96D15-plink2.gz | head -1
  zgrep rs113886122 bgen/Q96D15-plink2.gz
  zgrep rs113886122 bgen/P55899-plink2.gz
  gunzip -c bgen2/RCN3_All_invn-plink2.gz | head -1
  zgrep rs113886122 bgen2/RCN3_All_invn-plink2.gz
}

function qq_manhattan_plot()
{
  parallel -j1 --env caprion --env dir -C' ' '
    echo {}
    cut -d" " -f1,4,12  --output-delimiter=" " ${dir}/{}/interval.{}.All.txt | \
    grep -v NA > ${caprion}/${dir}/plot/{}.txt
    (
      echo chromosome position
      awk "\$3<5e-8{print \$1,\$2}" ${dir}/plot/{}.txt
    ) > ${caprion}/${dir}/plot/{}.annotate
    R --slave --vanilla --args \
      input_data_path=${caprion}/${dir}/plot/{}.txt \
      output_data_rootname=${dir}/plot/{}_turboqq \
      plot_title="{}" < ~/cambridge-ceu/turboqq/turboqq.r
    R --slave --vanilla --args \
      input_data_path=${caprion}/${dir}/plot/{}.txt \
      output_data_rootname=${dir}/plot/{}_turboman \
      custom_peak_annotation_file_path=${caprion}/${dir}/plot/{}.annotate \
      reference_file_path=~/cambridge-ceu/turboman/turboman_hg19_reference_data.rda \
      pvalue_sign=5e-8 \
      plot_title="{}" < ~/cambridge-ceu/turboman/turboman.r
    # rm ${caprion}/${dir}/plot/{}.txt
  '
}

export dir=RCN3-FCGRT-genomewide
export dirplot=${dir}/plot
if [ ! -d ${dirplot} ]; then mkdir ${dirplot}/plot; fi
ls ${dir} | grep "/" | sed 's|/||' | grep -v -e plot | qq_manhattan_plot

R --no-save -q <<END
  eda <- function(d,v)
  {
    plot(d[[v]],main=v,xlab="",ylab="")
    hist(d[[v]],main="",xlab="")
    boxplot(d[[v]],horizontal=TRUE)
  }
  setwd("/home/jhz22/Caprion/pilot")
  phase2 <- read.table("data2/phase2.dat",header=TRUE)
  phase2_invn <- read.table("data2/phase2_invn.dat",header=TRUE)
  pilot <- read.table("data/caprion.dat",header=TRUE)
  pdf("RCN3-FCGRT-genomewide/plot/eda.pdf")
  par(mfrow=c(3,1))
  eda(pilot,"Q96D15")
  eda(pilot,"Q96D15_invn")
  eda(pilot,"P55899")
  eda(pilot,"P55899_invn")
  eda(phase2,"RCN3_442625488_VADQDGDSMATR")
  eda(phase2_invn,"RCN3_442625488_VADQDGDSMATR_invn")
  eda(phase2,"RCN3_442666668_EVAKEFDQLTPEESQAR")
  eda(phase2_invn,"RCN3_442666668_EVAKEFDQLTPEESQAR_invn")
  eda(phase2,"RCN3_All")
  eda(phase2_invn,"RCN3_All_invn")
  eda(phase2,"RCN3_DR")
  eda(phase2_invn,"RCN3_DR_invn")
  dev.off()
END
