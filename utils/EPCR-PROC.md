# EPCR-PROC in Caprion pilot study

Working directory,

~/rds/projects/olink_proteomics/scallop/Caprion/EPCR-PROC

## Annotations
```
#    Protein Accession  Gene            Protein.Description
# EPCR_HUMAN    Q9UNN8 PROCR Endothelial protein C receptor
# PROC_HUMAN    P04070  PROC Vitamin K-dependent protein C

# chrom chromStart  chromEnd strand    acc uniprotName geneName geneSynonyms hgncSym         ensGene
# chr20   33759957  33764613      + Q9UNN8  EPCR_HUMAN    PROCR         EPCR   PROCR ENSG00000101000
#  chr2  128177518 128186519      + P04070  PROC_HUMAN     PROC                 PROC ENSG00000115718
```

# Phase I data

## Correlation of piptides

* EPCR-PROC.tsv. The raw data.
* EPCR-PROC-corr.tsv. The Pearson correlation.
* EPCR-PROC-corr.png. The heatmap.

![Heatmap](EPCR-PROC/EPCR-PROC-corr.png)

```
                          names.EPCR.
1 EPCR_442603139_LHM.147.0354.LQISYFR
2 EPCR_442605396_LHM.147.0354.LQISYFR
3           EPCR_442582461_LHMLQISYFR
4            EPCR_442581804_TLAFPLTIR
5                            EPCR_All
6                             EPCR_DR

               EPCR_442603139 EPCR_442605396 EPCR_442582461 EPCR_442581804 EPCR_All  EPCR_DR
EPCR_442603139 
EPCR_442605396 "-0.0546"
EPCR_442582461 " 0.2541"      " 0.0815"
EPCR_442581804 " 0.3482"      " 0.0190"      " 0.3011"
EPCR_All       " 0.2660"      " 0.9359"      " 0.2662"      " 0.2138"
EPCR_DR        " 0.2660"      " 0.9359"      " 0.2662"      " 0.2138"      " 1.0000"

Call:
lm(formula = EPCR_All ~ EPCR_442603139 + EPCR_442605396 + EPCR_442582461 +
    EPCR_442581804, data = EPCR)

Residuals:
      Min        1Q    Median        3Q       Max
-0.060649 -0.010128  0.000392  0.011153  0.056049

Coefficients:
                Estimate Std. Error t value Pr(>|t|)
(Intercept)    8.7522184  0.0853057  102.60   <2e-16 ***
EPCR_442603139 0.0111745  0.0002604   42.92   <2e-16 ***
EPCR_442605396 0.4346305  0.0026603  163.38   <2e-16 ***
EPCR_442582461 0.0470646  0.0028821   16.33   <2e-16 ***
EPCR_442581804 0.0323885  0.0027846   11.63   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.01921 on 195 degrees of freedom
Multiple R-squared:  0.9936,    Adjusted R-squared:  0.9935
F-statistic:  7588 on 4 and 195 DF,  p-value: < 2.2e-16
```

## Scatter, histogram and boxplots

![EPCR](EPCR-PROC/EPCR-PROC-desc-0.png)
![PROC](EPCR-PROC/EPCR-PROC-desc-1.png)

## Q-Q/Manhattan plots

Associaiton statistics were based on invn(ormal) transformation, invnormal(prot)~genotypes.

![EPCR](EPCR-PROC/Q9UNN8_invn_turboqq.png)
![EPCR](EPCR-PROC/Q9UNN8_invn_turboman.png)
![PROC](EPCR-PROC/P04070_invn_turboqq.png)
![PROC](EPCR-PROC/P04070_invn_turboman.png)

## Sentinels

Sentinels are Q9UNN8_invn.sentinels and P04070_invn.sentinels.

### EPCR

MarkerName|P|prot|rsid
----------|-|----|----
chr1:198680470|7.7746e-06|Q9UNN8_invn|rs150905683
chr1:242492519|4.52035e-06|Q9UNN8_invn|rs4658474
chr2:40640161|8.47293e-06|Q9UNN8_invn|rs6742589
chr2:108904323|1.79526e-06|Q9UNN8_invn|rs71421322
chr2:129343910|3.38733e-06|Q9UNN8_invn|rs538566482
chr3:3589208|6.97483e-06|Q9UNN8_invn|rs73010413
chr3:63622010|4.92919e-06|Q9UNN8_invn|rs73126516
chr3:115857670|8.18486e-06|Q9UNN8_invn|rs36060209
chr4:26467624|3.08996e-06|Q9UNN8_invn|rs78545240
chr5:117693310|2.17759e-06|Q9UNN8_invn|rs73791119
chr5:178789756|6.07192e-06|Q9UNN8_invn|rs11749206
chr6:10437912|8.87339e-06|Q9UNN8_invn|rs76522721
chr6:31385204|5.65772e-07|Q9UNN8_invn|rs113425623
chr6:150708928|6.90209e-06|Q9UNN8_invn|rs9478711
chr7:46878807|2.58103e-06|Q9UNN8_invn|rs73119156
chr8:113975571|7.68998e-06|Q9UNN8_invn|rs10505197
chr8:141604230|1.29017e-06|Q9UNN8_invn|rs11997990
chr11:8066317|2.82897e-06|Q9UNN8_invn|rs72848407
chr11:82121298|8.22729e-06|Q9UNN8_invn|rs7114139
chr12:119574662|3.5405e-06|Q9UNN8_invn|rs12229183
chr14:45522256|5.26276e-07|Q9UNN8_invn|rs374947851
chr16:54065346|8.6531e-06|Q9UNN8_invn|rs12051239
chr16:86878088|1.07322e-06|Q9UNN8_invn|rs11645358
chr18:65690741|3.02894e-06|Q9UNN8_invn|rs17078446

### PROC

MarkerName|P|prot|rsid
----------|-|----|----
chr1:14700486|4.41035e-06|P04070_invn|rs1539170
chr1:226790403|4.91707e-06|P04070_invn|rs878149
chr2:85379656|1.95094e-06|P04070_invn|rs2568224
chr3:39025856|2.77768e-06|P04070_invn|rs78768764
chr3:70686181|5.40327e-06|P04070_invn|rs9826254
chr3:106597199|8.58278e-07|P04070_invn|rs56400218
chr6:3251943|7.99169e-06|P04070_invn|rs11963805
chr6:7937069|3.63501e-06|P04070_invn|rs776682351
chr6:145369170|3.93743e-06|P04070_invn|rs117840209
chr7:136707968|5.2286e-06|P04070_invn|rs10246819
chr8:7128386|9.31995e-07|P04070_invn|rs370324853
chr8:99183900|5.71606e-06|P04070_invn|rs2443574
chr9:138730664|2.63476e-06|P04070_invn|rs141376886
chr10:54410153|3.65163e-06|P04070_invn|rs1348351
chr10:70244181|4.25257e-06|P04070_invn|rs35608830
chr11:13498427|8.78973e-07|P04070_invn|rs11022845
chr11:13498427|8.78973e-07|P04070_invn|rs556750517
chr11:23539628|9.03995e-06|P04070_invn|rs7943543
chr12:133801548|9.85288e-06|P04070_invn|rs556045604
chr13:67082089|5.07758e-06|P04070_invn|rs148386003
chr13:81938729|4.35324e-06|P04070_invn|rs117099502
chr13:104582407|6.21026e-06|P04070_invn|rs11842347
chr14:93804153|4.33826e-07|P04070_invn|rs117720114
chr15:65003284|1.67562e-06|P04070_invn|rs144408155
chr18:12739631|9.43079e-06|P04070_invn|rs62096050
chr20:36889971|8.93867e-06|P04070_invn|rs6098808
chr20:36889971|8.93867e-06|P04070_invn|rs767174705
chr21:34516707|4.38852e-06|P04070_invn|rs963950

# Phase II data

```
                                     EPCR_All   EPCR_DR
EPCR_442581804_TLAFPLTIR            0.5113919 0.6174716
EPCR_442582461_LHMLQISYFR           0.2986443 0.4357448
EPCR_442603139_LHM[147.0354]LQISYFR 0.7980766 0.7396526
EPCR_442605396_LHM[147.0354]LQISYFR 0.6105076 0.5631775
EPCR_All                            1.0000000 0.9822826
EPCR_DR                             0.9822826 1.0000000
```

![Desc1](EPCR-PROC/EPCR-PROC-phase2-0.png)
![Desc2](EPCR-PROC/EPCR-PROC-phase2-1.png)

![Pattens of correlation](EPCR-PROC/EPCR-PROC-phase2-all.png)
![Correlations](EPCR-PROC/EPCR-PROC-phase2-corr.png)
