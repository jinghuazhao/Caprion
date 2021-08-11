local dir : env HPC_WORK
local caprion : env caprion
local y : env y
local trait : env trait
set maxvar 21000
insheet using "`caprion'/data2/phase2_invn.dat", delim(" ") case clear
// To fix the glitch of SNPTEST on covariates
foreach x of varlist PC1-PC20 {
   replace `x'=`x'*10000
}
rename ID_1 id
format id %15.0g
merge 1:1 id using "`dir'/data/interval/interval_data"
gwas2 `trait', studyname(interval) dirgenefiles("`caprion'/data2") chr(19) /*
*/    covariates(agePulse bmi PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 PC11 PC12 PC13 PC14 PC15 PC16 PC17 PC18 PC19 PC20) /*
*/    class(sexPulse ethnicity classification) dirwork("`caprion'/`y'") outfmt(txt)
