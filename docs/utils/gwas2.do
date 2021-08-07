//
// export TMPDIR=$HPC_WORK/work
// 

local dir : env HPC_WORK
set maxvar 21000
insheet using data2/phase2.pheno, delim(" ") case clear
rename IID id
format id %15.0g
merge 1:1 id using "`dir'/data/interval/interval_data"
gwas2 FIBB_44258003~n, studyname(interval) dirgenefiles("`dir'/data/interval")
