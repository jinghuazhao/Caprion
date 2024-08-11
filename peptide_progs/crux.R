#!/usr/local/Cluster-Apps/ceuadmin/R/4.4.1-icelake/bin/Rscript

# NOTE: The usage and distribution of this software is governed by the Apache License Version 2.0  (https://www.apache.org/licenses/LICENSE-2.0)

library(multicomp)

fasta_fi1e=Sys.getenv("uniprot")
spectra_file=Sys.getenv("mgf")
scoring="xcorr"
n_decoys=5
fdr_threshold=0.01
output_file=Sys.getenv("out")
keep_psms=TRUE
name="ZWK"
extra_crux="--precursor-window 10 --precursor-window-type ppm --mz-bin-width 0.02 --pm-min-peak-pairs 100 --pm-charges 2"
seed=123

fasta = paste0(fasta_fi1e,".fasta")
ms2 = paste0(spectra_file,".mzML")

CRUX=Sys.which("crux")
val = 0

dataSet = name

# Preprocessing

system(paste("mkdir -p ", dataSet, sep = ""))

for(decoy in 0:(n_decoys - 1))
{
  index = paste(dataSet, "/dcy", decoy, ".index", sep = "")
  if(!file.exists(index))
  {
    system(paste(CRUX, " tide-index --allow-dups T --seed ", decoy," --overwrite T --output-dir ", index, " ", fasta, " ", index, sep = ""))
  }
  
  if(scoring %in% c("xcorr", "pvalue"))
  {
    root = paste(dataSet, ".", scoring, sep = "")
    if(scoring == "xcorr")
    {
      params = "" 
    }
    else
    {
      params="--exact-p-value T --mz-bin-width 1.0005079"
    }
    tideResults = paste(root, "/dcy", decoy, ".tide-search.decoy.txt", sep = "")
    if(!file.exists(tideResults))
    {
      system(paste(CRUX, " tide-search --output-dir ", root, " --fileroot dcy", decoy, " --overwrite T --concat F --top-match 1 --auto-precursor-window warn --auto-mz-bin-width warn ", extra_crux, " ", params, " ", ms2, " ", index, sep = ""))
    }
    
    if(decoy == 0)
    {
      tideResults = paste(root, "/dcy", decoy, ".tide-search.target.txt", sep = "")
      for(method in c("tdc", "mix-max"))
      {
        confidenceResults = paste(root, "/", method, ".dcy", decoy, ".assign-confidence.target.txt", sep = "")
        if(!file.exists(confidenceResults))
        {
        # extract-rows has been removed from 4.0
        # system(paste(CRUX, " extract-rows --comparison eq --column-type int ", tideResults, " \"xcorr rank\" 1 > target.tmp", sep = ""))
        # system(paste(CRUX, " extract-rows --comparison eq --column-type int \`echo ", tideResults, " | sed \'s/target/decoy/\'\` \"xcorr rank\" 1 > decoy.tmp", sep = ""))
          system(paste(CRUX, " assign-confidence --estimation-method ", method, " --output-dir ", root, " --fileroot ", method, ".dcy", decoy, " target.tmp", sep = ""))
        # system("rm target.tmp decoy.tmp")
        }
      }
    }
  }else{
	  print("ERROR: Invalid scoring method, please use \"xcorr\" or \"pvalue\"")
	}
}
	
# Remove duplicated targets
if(n_decoys > 1)
{
	for(decoy in 1:(n_decoys - 1))
	{
  		system(paste("rm -f */dcy", decoy, ".tide-search.target.txt", sep = ""))
	}

	for(decoy in 0:(n_decoys - 1))
	{
  		system(paste("rm -f */dcy", decoy, ".tide-search.log.txt", sep = ""))
 		system(paste("rm -f */dcy", decoy, ".tide-search.params.txt", sep = ""))
	}
}

if(scoring == "xcorr")
{
  is_xcorr = T
} else if(scoring == "pvalue")
{
  is_xcorr = F
} else
{
  print("Invalid scoring type. Defaulting to xcorr score")
  is_xcorr = T
}


get_peptide_scores = function(dat_list)
{
  tar = dat_list[[1]]
  total_peps = tar$sequence
  for(i in 2:length(dat_list))
  {
    total_peps = c(total_peps, dat_list[[i]]$original.target.sequence)
  }
  uniq_total_peps = unique(total_peps)
  
  tar_scores_sequence = uniq_total_peps
  tar_scores_scores = c()
  dec_score_matrix = matrix(-100, nrow = length(uniq_total_peps), ncol = length(dat_list) - 1)
  
  for(i in 1:length(uniq_total_peps))
  {
    if(uniq_total_peps[i] %in% tar$sequence)
    {
      tar_scores_scores[i] = max(tar[which(tar$sequence == uniq_total_peps[i]),]$xcorr.score) 
    }
    else{
      tar_scores_scores[i] = -100
    }
    
    for(j in 2:length(dat_list))
    {
      if(uniq_total_peps[i] %in% dat_list[[j]]$original.target.sequence)
      {
        dec_score_matrix[i, j-1] = max(dat_list[[j]][which(dat_list[[j]]$original.target.sequence == uniq_total_peps[i]),]$xcorr.score) 
      }
      else
      {
        dec_score_matrix[i, j-1] = -100
      }
    }
  }
  
  tar_scores = data.frame(sequence = tar_scores_sequence, xcorr.score = tar_scores_scores)
  return(list(tar_scores, dec_score_matrix))
}


read_all = function(min, max, data_name, use_xcorr)
{
  if(use_xcorr)
  {
    df.list = list()
    df.list[[1]] = read.csv(paste("./", data_name, ".xcorr/dcy0.tide-search.target.txt", sep = ""), sep = "\t", stringsAsFactors = F)
    for(i in min:max)
    {
      df.list[[i+2]] = read.csv(paste("./", data_name, ".xcorr/dcy", toString(i),".tide-search.decoy.txt", sep = ""), sep = "\t", stringsAsFactors = F)
    }
    return(df.list) 
  }
  else{
    df.list = list()
    df.list[[1]] = read.csv(paste("./", data_name, ".pvalue/dcy0.tide-search.target.txt", sep = ""), sep = "\t", stringsAsFactors = F)
    for(i in min:max)
    {
      df.list[[i+2]] = read.csv(paste("./", data_name, ".pvalue/dcy", toString(i),".tide-search.decoy.txt", sep = ""), sep = "\t", stringsAsFactors = F)
    }
    return(df.list) 
  }
}

dat_list = read_all(0, n_decoys - 1 , dataSet, is_xcorr) 

if(!is_xcorr)
{
    for(i in 1:length(dat_list))
    {
	old_name = names(dat_list[[i]])
	old_name[which(old_name == "exact.p.value")] = "xcorr.score"
	dat_list[[i]] = setNames(dat_list[[i]], old_name)
	dat_list[[i]]$xcorr.score = -1* dat_list[[i]]$xcorr.score
    }
}

if(seed != "NULL")
{
  set.seed(as.integer(seed))
}

pep_scores = get_peptide_scores(dat_list)

tar_scores = matrix(pep_scores[[1]]$xcorr.score, nrow = 1)
dec_scores = t(pep_scores[[2]])

alpha_range = as.double(fdr_threshold)
res_index = multidecoy_comp(tar_scores, dec_scores, alpha_range)[[1]]
res_peps = pep_scores[[1]][res_index,]$sequence

write(as.character(res_peps), output_file)

if(keep_psms == "F")
{
   system(paste("rm -r ", dataSet, ".xcorr", sep = ""))
   system(paste("rm -r ", dataSet, ".pvalue", sep = ""))
   system(paste("rm -r ", dataSet, sep = ""))
}
