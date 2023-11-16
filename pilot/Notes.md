# Notes

## m/z (ChatGPT)

In the context of mass spectrometry, the term "charge" refers to the electrical charge of an ion. Mass spectrometry is a technique used to analyze the mass-to-charge ratio (m/z) of ions. Monoisotopic m/z refers to the mass-to-charge ratio of the monoisotopic peak, which is the peak corresponding to the ion containing the most abundant isotope of each element.

The Max Isotope Time Centroid is the time at which the intensity-weighted average of the isotope distribution reaches its maximum value. This parameter is important for accurately determining the mass and charge state of an ion.

The charge of an ion is a crucial piece of information in mass spectrometry. It is used in the calculation of the mass-to-charge ratio (m/z) of an ion, which is a key parameter in identifying and characterizing molecules. The m/z value is used to determine the mass of an ion, and the charge is needed to convert the raw mass spectrometry data (measured in daltons) into the mass-to-charge ratio.

The formula for calculating m/z is:

$$
m/z=\frac{\mbox{Mass of the ion}}/{\mbox{Mass of the ion​}}
$$

The charge of an ion is essential for interpreting mass spectrometry data, particularly in determining the mass-to-charge ratio (m/z) of ions, which is fundamental for identifying and characterizing molecules in a sample.

When the charge is 0 for monoisotopic m/z, it implies a neutral species, and the mass spectrometry analysis might be focusing on the measurement of neutral molecules or particles.

In proteomics, researchers often analyze proteins and peptides in their neutral forms. Techniques like matrix-assisted laser desorption/ionization (MALDI) are commonly used for the analysis of intact proteins in their neutral state

## Peptide and Protein ID using OpenMS tools

<https://training.galaxyproject.org/training-material/topics/proteomics/tutorials/protein-id-oms/tutorial.html#peptide-identification>

## PoGo

Fast Mapping of Peptides to Genomic Coordinates for Proteogenomic Analyses, <https://www.sanger.ac.uk/tool/pogo/>, GitHub, <https://github.com/cschlaffner/PoGo>.

It uses transcript translations and reference gene annotations to identify the genomic loci of peptides and post-translational modifications. Multiple occurrences of peptides in the input data resulting in the same genomic loci will be collapsed as a single occurrence in the output.

The input format is a tab delimited file with four columns with file extensions such as *.pogo, *.txt, and *.tsv.

Column | Column header | Description
-------|---------------|------------
1      | Sample        | Name of sample or experiment
2      | Peptide       | Peptide sequence with PSI-MS nodification names in round brackets following the mpdified amino acid, e.g. PEPT(Phopsho)IDE for a phosphorylated threonine
3      | PSMs          | Number of peptide-spectrum matches (PSMs) for the given peptide, including those redundantly identified (peptides can be “seen” more than once in a run)
4      | Quant         | Quantitative value for the given peptide in the given sample

An example is established as follows,

```bash
wget -S ftp://ftp.sanger.ac.uk/pub/teams/17/software/PoGo/PoGo_Testprocedures.zip
unzip PoGo_Testprocedures.zip
cd PoGo_Testprocedures/Testfiles
module load ceuadmin/PoGo
for Peptides in Testfile_experimental Testfile_small
do
  PoGo -fasta input/gencode.v25.pc_translations.fa -gtf input/gencode.v25.annotation.gtf -in input/${Peptides}.txt
done
```

Output files are also contained in the `input/` directory.

GENCODE annotation data are available from <https://www.gencodegenes.org/human/> and <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/>.

The Java GUI, <https://github.com/cschlaffner/PoGoGUI>, is run as follows,

```bash
java -jar PoGoGUI-v1.0.0.jar
```
which requires PoGo executable as well.

The source is compiled with maven, <https://maven.apache.org/>, e.g.,

```bash
module load maven-3.5.0-gcc-5.4.0-3sgaeze
mvn install
```
assuming that `pom.xml` is available, e.g., `/usr/local/Cluster-Apps/ceuadmin/PoGo/1.0.0/PoGoGUI/PoGoGUI`.

## Proteoform Analysis

### ETH Zurich, U Toronto Team Develops Tool for Bottom-Up Proteomics Proteoform Analysis

Jul 28, 2021 | Adam Bonislawski

NEW YORK – A team led by researchers at ETH Zurich and the University of Toronto has developed a tool that allows for the detection of protein proteoforms in bottom-up proteomics data.

Described in a paper published in June in Nature Communications, the tool, called COPF (COrrelation-based functional ProteoForm) uses peptide correlation analysis to detect differences in proteoform populations across different samples or conditions and could aid researchers as they seek to better understand the role different protein forms play in biology and disease.

The human genome is thought to have around 20,000 protein coding genes, but many of these 20,000 proteins exist in the body in various forms, differentiated by, for instance, post-translational modifications or amino acid substitutions. These different forms are called proteoforms, and it is widely believed that biological processes are guided not just by the proteins present but by what proteoforms are present and in what proportions.

Traditional bottom-up proteomics workflows have provided only limited insight into proteoform populations however, due to the fact that the presence of a particular protein is typically inferred by the detection of just a few of its peptides and that digesting proteins into peptides for mass spec analysis makes it near impossible to link a modified peptide back to a particular proteoform.

Some proteomics researchers have addressed this issue by moving to top-down proteomics, which looks at intact proteins, allowing them to better distinguish between different proteoforms. Top-down proteomics is very technically challenging, however, and is not yet able to analyze proteins with the breadth and depth of bottom-up workflows.

Recently, the development of more reproducible and higher-throughput bottom-up workflows, and particular workflows using data independent-acquisition (DIA) mass spectrometry, have allowed researchers like the Nature Communications authors to apply peptide correlation analysis to the study of proteoforms.

Peptide correlation analysis looks at differences in peptide behavior within and across proteins in bottom-up data.

Researchers have developed a number of approaches for turning peptide measurements into protein data, with most working under the assumption that peptides from the same protein will behave the same way.

In practice, though, that isn't the case. On one hand, there are a number of technical reasons why two peptides from the same protein may not behave the same way. For instance, different digestion efficiencies could lead to some peptides being more abundant than others. Different ionization efficiencies could similarly make one peptide more likely than another to be detected by the mass spec.

The presence of different proteoforms could also play a role. For instance, if a protein is present in both a full-length and truncated form, expression changes observed in the full-length form wouldn't be observable if the peptide being measured wasn't present in the truncated form. Not only would this throw off protein-level quantitation, but it would also mask relative changes in the two protein forms that could be biologically important.

A major challenge to applying this insight has been determining which differences in peptide behavior reflect real technical or biological variation and which are just noise, noted Hannes Röst, research chair in mass spectrometry-based personalized medicine at the University of Toronto and an author on the Nature Communications study.

"In many cases [such variation] was noise," he said. "When you look at traditional shotgun proteomics workflows and data analyses, really the power is not at the peptide-level quantification but at the protein level from the aggregation of multiple peptides. On the peptide level you see a lot of noise, and I think that has prevented us from using this observation that individual peptides could yield a lot of interested information because people really only looked at the protein-level data, because that is what they trusted."

Röst said that the development of targeted protein quantitation approaches like multiple-reaction monitoring (MRM) has demonstrated that individual peptides can be measured with high accuracy, and the development of DIA mass spec approaches has enabled MRM-style peptide quantitation at the proteome scale.

At the same time, improvements in mass spec technology have allowed researchers to collect the kind of large and reproducible datasets required for peptide correlation analysis, he said.

"These are types of experiments we wouldn't have imagined 10 years ago, because for correlation-based approaches to work, you need a relatively large number of samples, and you need low variance," he said. "We are not detecting [proteoforms] that are not changing between different [conditions], we are only detecting those that change. And for this to work we need to have multiple replicates and we need to have different conditions and to be able to measure these peptides with high quantitative accuracy across these conditions."

The COPF tool looks at the intensities of peptides coming from a particular protein across all the samples measured in an experiment and then calculates peptide correlations for all the pairs of peptides coming from that protein and uses hierarchical clustering to divide the peptides into two clusters. It then scores the likelihood that multiple proteoforms of a protein are present by comparing the level of peptide correlation between the clusters to the level of in-cluster variation.

The tool does not identify the specific modifications or variations that distinguish the different proteoforms but rather the peptides that appear to differentiate between the forms of the protein in the different biological contexts investigated. Analyzing a DIA dataset that looked at five different tissue types across eight different mice, COPF identified 63 proteins that exhibited different proteoform groups, including proteins with known tissue-specific splice variants. The researchers also identified proteoforms created by proteolytic and autocatalytic cleavage and phosphorylation, indicating, they wrote, that the tool is "agnostic to the different mechanisms by which proteoforms can be generated inside the cell."

The development of COPF follows the publication last year of a study by researchers at Barts Cancer Institute and the University of Wisconsin-Madison detailing another peptide correlation analysis tool for identifying proteoforms in bottom-up data called PeCorA.

Unlike COPF, which requires proteoforms to differ by two or more peptides, PeCorA can detect proteoforms based on single peptide differences. This makes it a potentially more sensitive tool but also less specific than COPF, Röst said.

More generally, he said that he expected ongoing improvements in mass spec technology would further improve peptide correlation-based approaches like COPF and PeCorA by boosting peptide coverage.

"To kind of cover every possible protein isoform we would need to have complete coverage of every protein, and unfortunately we are currently quite far away from having peptide-level coverage of every protein," he said. "I think that is currently one of the limitations where we are kind of hitting a wall."

Röst added that his lab has begun acquiring data on Bruker's timsTOF Pro platform, "and there we definitely see both an increase in protein coverage and also in the number of peptides we can measure."

"That's why I'm very optimistic that while this is just the first implementation of the method, the data we are producing at this moment is much more complete, and therefore I think it would be even more suitable to our approach than the data we used in the paper," he said.
