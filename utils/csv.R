load("2020.rda")
write.csv(Comp_Neq1,file="ZYQ_Comp_Neq1_Norm_Int_20200812.csv", quote=FALSE, row.names=FALSE)
write.csv(Normalized_All, file="ZYQ_Protein_Norm_All_20200813_v1.csv", quote=FALSE, row.names=FALSE)
write.csv(Protein_DR_Filt, file="ZYQ_Protein_Norm_DR_filt_20200813_v1.csv", quote=FALSE, row.names=FALSE)

