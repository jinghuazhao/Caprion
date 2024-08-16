# Notes

(**Sections I -- III are due to Claude**)

## I. Meta-data

1. **Isotope.Group.ID** is a unique identifier for a group of isotopes that belong to the same peptide or molecule. In mass spectrometry, isotopes are atoms of the same element that have the same number of protons but differ in the number of neutrons. This ID helps to group together isotopes that arise from the same peptide, allowing for easier identification and analysis.
2. **Protein** contains the name or identifier of the protein that the peptide (or molecule) is derived from. This information is typically obtained by searching the MS data against a protein database.
3. **Modified.Peptide.Sequence** is the amino acid sequence of the peptide, including any post-translational modifications (PTMs) that have been identified. PTMs are chemical modifications that occur after protein synthesis, such as phosphorylation, ubiquitination, or methylation. The sequence is usually represented in a standard format, such as using lowercase letters for modified residues.
4. **Monoisotopic.m/z** is the monoisotopic mass-to-charge ratio (m/z) of the peptide or molecule. The monoisotopic mass is the mass of the most abundant isotope of each element in the molecule, which is typically the lightest isotope (e.g., 12C, 1H, 14N, 16O, etc.). This value is used as a reference point for identifying the peptide or molecule.
5. **Max.Isotope.Time.Centroid** is the time centroid (or apex) of the most intense isotope in the isotope group. In liquid chromatography-mass spectrometry (LC-MS), peptides are separated based on their retention time (the time it takes for the peptide to elute from the column). The time centroid is the time point at which the peptide signal is most intense, which can be used to quantify the peptide abundance.
6. **Charge** is the charge state of the peptide or molecule. In mass spectrometry, peptides can be ionized to different charge states (e.g., +1, +2, +3, etc.), which affects their mass-to-charge ratio (m/z). The charge state is an important parameter for identifying peptides and molecules.

They are invaluable for analyzing and interpreting MS data, including peptide identification, quantification, and characterization of post-translational modifications.

## II. MS1/MS2

In mass spectrometry-based proteomics, the typical workflow for identifying peptides and proteins involves using tandem mass spectrometry (MS/MS or MS2). In this process, precursor ions (peptides) are selected in the first stage of mass spectrometry (MS1) and then fragmented to produce a series of smaller ions in the second stage (MS2). The resulting fragment ions (product ions) are analyzed to infer the sequence of the peptide and, by extension, identify the proteins from which they were derived.

However, it is possible to infer peptides and proteins using only MS1 data through a process known as "MS1-only" or "untargeted" analysis. This approach can be particularly useful in the following scenarios:

1. **Label-based quantification**: Techniques like SILAC (Stable Isotope Labeling by Amino acids in Cell culture) or chemical labeling (e.g., TMT, iTRAQ) rely on MS1 data for quantification. The mass shift introduced by labels allows for the direct comparison of peptide abundances based on their MS1 ion intensities.

2  **Label-free quantification**: Proteins can be quantified by comparing the intensities of their corresponding peptide ions in MS1 across different samples. This requires accurate mass and retention time alignment and often uses algorithms to detect and quantify features (peptide ions) consistently across multiple runs.

3. **Accurate Mass and Time tags (AMT)**: This approach relies on a previously established library of peptide identifications, where each peptide is characterized by its accurate mass and normalized retention time. In subsequent analyses, peptides can be inferred by matching the observed accurate mass and retention time to the library without the need for MS2 fragmentation.

4. **Data-independent acquisition (DIA)**: In some DIA workflows, proteins can be inferred from MS1 data when coupled with complex data analysis strategies and spectral libraries. It is important to note that while DIA collects MS1 spectra, it also involves the simultaneous fragmentation of all ions in a given mass range, and thus MS2-level data is typically available and used for identification.

It is important to note that MS1-only approaches may have limitations in terms of identification specificity and sensitivity compared to traditional MS2-based methods. MS1-based protein inference is generally less confident because it lacks sequence-specific information that can only be obtained from fragment ions in MS2. For this reason, MS1-based methods are often complemented by MS2 data or rely on extensive peptide libraries and sophisticated computational algorithms to increase the confidence of peptide and protein identification.

## III. OpenMS/crux/maxQuant/FragPipe

OpenMS, Crux, MaxQuant, and FragPipe are all prominent software platforms for analyzing proteomics data, each with its own strengths and weaknesses. Here's a comparison:

* OpenMS:

    * Focus: Provides a flexible and open-source framework for developing and executing various mass spectrometry data analysis workflows.
    * Strengths:
          *  Highly modular and customizable: Offers a vast collection of algorithms and tools that can be combined and customized to create tailored workflows.
          *  Open-source and extensible: Encourages community contributions and allows for the development of new tools and algorithms.
          *  Supports various data formats and instruments: Compatible with a wide range of data formats and mass spectrometry platforms.
          *  Strong support for metabolomics data: While primarily used for proteomics, it also offers tools for analyzing metabolomics data.
    * Limitations:
          *  Steeper learning curve: Requires programming knowledge and familiarity with command-line interfaces.
          *  Less user-friendly: Lacks a comprehensive graphical user interface (GUI), making it less accessible for beginners.
          *  Limited pre-built workflows: While highly customizable, it requires more manual effort to set up standard workflows compared to MaxQuant or FragPipe.

* Crux:

    * Focus: A command-line toolkit designed for peptide identification, protein quantification, and statistical validation of proteomics data.
    * Strengths:
          *  Open-source and well-documented: Provides clear documentation and allows for community contributions.
          *  Fast and efficient: Known for its computational efficiency and speed.
          *  Strong statistical validation: Offers rigorous statistical methods for validating peptide and protein identifications.
          *  Supports various search engines: Compatible with multiple search engines, including Comet and Tide.
    * Limitations:
          *  Command-line interface only: Requires familiarity with command-line operations.
          *  Less user-friendly: Lacks a GUI, making it less accessible for beginners.
          *  Limited pre-built workflows: Requires more manual effort to set up complete analysis pipelines.

* MaxQuant:

    * Focus: Primarily known for its robust and sensitive peptide and protein identification and quantification using its proprietary Andromeda search engine.
    * Strengths:
          *  User-friendly interface: Provides a GUI for easier data processing and analysis.
          *  Robust and sensitive identification and quantification: Offers high-quality results for standard DDA-based proteomics experiments.
          *  Strong support for label-free quantification (LFQ) and match between runs (MBR).
          *  Extensive post-translational modification (PTM) analysis: Offers comprehensive support for identifying and quantifying various PTMs.
    * Limitations:
          *  Less flexible for specialized workflows: Primarily designed for standard bottom-up proteomics experiments.
          *  Limited support for DIA data: While it can handle DIA data, it's not its primary strength.
          *  Closed-source: The core algorithms are not open-source, limiting community contributions and customization.

* FragPipe:

    * Focus: Offers a more modular and flexible platform with various tools for different proteomics workflows, including both DDA and DIA.
    * Strengths:
          *  Versatile and modular: Includes a suite of tools for various tasks, including peptide identification, quantification, and statistical analysis.
          *  Extensive support for DIA data: Features DIA-Umpire, a dedicated tool for analyzing DIA data using various algorithms.
          *  Highly accurate and sensitive quantification: Employs IonQuant for precise quantification using extracted ion chromatograms.
          *  Open-source and actively developed: Encourages community contributions and continuous improvement.
    * Limitations:
          *  Steeper learning curve: Primarily operates through a command-line interface, requiring more technical expertise.
          *  Less user-friendly interface: Lacks a comprehensive GUI, making it less intuitive for beginners.

Here's a table summarizing the key differences:

Feature |OpenMS |Crux |MaxQuant |FragPipe
--------|-------|-----|---------|--------
Primary Focus |Flexible Framework |Peptide ID & Quantification |Peptide & Protein ID/Quant |Modular Platform
Open Source |Yes |Yes |No (core algorithms) |Yes
User Interface |Primarily CLI |CLI |GUI |Primarily CLI
Learning Curve |Steep |Moderate |Easier |Steep
Flexibility |Highly Flexible |Moderate |Less Flexible |More Flexible
DIA Support |Limited |Limited |Limited |Extensive (DIA-Umpire)
Quantification Methods |Various |Various |LFQ, iBAQ |IonQuant
PTM Analysis |Supported |Supported |Extensive |PTM-Shepherd
Community Support |Strong |Moderate |Moderate |Strong

Choosing the Right Tool:

* OpenMS: Ideal for researchers with programming skills who need a highly customizable and extensible platform for developing specialized workflows.
* Crux: Suitable for researchers comfortable with command-line interfaces and seeking fast and efficient tools for peptide identification, protein quantification, and statistical validation.
* MaxQuant: Best for researchers looking for a user-friendly platform with robust performance for standard DDA-based proteomics experiments, especially those focusing on label-free quantification.
* FragPipe: Ideal for researchers seeking a highly flexible and customizable platform for various workflows, including DIA analysis, and who are comfortable with command-line operations.

Remember to consider your specific research goals, data type, and bioinformatics expertise when choosing the best tool for your needs. You might even explore combining different tools to leverage their unique strengths for different aspects of your analysis.

## IV. Galaxy tutorials

Web: <https://usegalaxy.org/> (https://training.galaxyproject.org/training-material/)

```bash
docker run -p 8080:80 quay.io/galaxy/introduction-training
```

Visit <http://localhost:8080>. Login as **admin** with password **password** to access.

See <https://training.galaxyproject.org/training-material/topics/proteomics/tutorials/protein-id-oms/tutorial.html>

## V. PoGo

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

## VI. Proteoform Analysis

#### ETH Zurich, U Toronto Team Develops Tool for Bottom-Up Proteomics Proteoform Analysis

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
