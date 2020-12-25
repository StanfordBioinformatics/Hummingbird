## Section 6: Using Different Input File Formats and Tools for Format Conversions

In case users want to leverage the downsampling step in Hummingbird but have input files in formats different than BAM or fastq/fastq.gz, please follow the examples below:

1. Input file(s) are in CRAM format

   The user can convert CRAM to SAM using `samtools` (http://www.htslib.org/doc/samtools-view.html): ```samtools view -C -T ref.fa aln.bam > aln.cram``` 
   Please note that the original reference fasta file is required for this conversion. Generating the index file after conversion may be necessary for subsequent    analysis using software tools or pipelines. For more information on using CRAM files with samtools, please see http://www.htslib.org/workflow/.

2. Input file(s) are in SAM format

   `Samtools` has a functionality that does this conversion: ```samtools view -bS file.sam | samtools sort - file_sorted```
   Generating the index file after conversion may be necessary for subsequent analysis using software tools or pipelines.

3. Input file(s) are in FASTQ but need uBAM

   In some cases, the bioinformatics pipeline to evaluate accepts unaligned BAM files so conversion of the FASTQ files to uBAM is needed.
   FastqToSam tool within the Picard suite of tools (https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam) can be used: 
   ```java -jar picard.jar FastqToSam F1=file_1.fastq O=fastq_to_bam.bam SM=for_tool_testing```
      
4. Aligned BAM to unmapped BAM

   The RevertSam tool from the Picard tools suite can be used: ```java -jar picard.jar RevertSam I=input.bam O=reverted.bam``` where the `input.bam` is the          aligned BAM and `reverted.bam` is the output unmapped BAM.

5. BED to BAM 

   `Bedtools` can be used to convert a file from BED format to BAM: ```bedtools bedToBam -i input.bed -g genome_file > input_converted.bam``` where `genome_file`    is not a fasta file but a two column file with list of chromosomes and the corresponding chromosome sizes in basepairs. 

   The `genome_file` can be fasta index file in the `.fai` format where the first two columns are extracted (indexing can be done using ```samtools faidx            reference.fasta```) and the "chr" prefix is added to the chromosome names. For more details on `genome_file` and pre-defined genome files available with          bedtools distribution, please see https://bedtools.readthedocs.io/en/latest/content/general-usage.html#genome-file-format.

   In case the input file is in BED12 format and spliced BAM entries are to be generated, use: ```bedToBam -i input_bed12format.bed -g genome_file -bed12 >          input_converted_spliced.bam```

   If required, BED12 file (has blocked features) can be converted to BED6 (each feature listed in a separate line) format using: ```bedtools bed12ToBed6 -i          input_bed12.bed```. For details on options, please refer to bedtools documentation (https://bedtools.readthedocs.io/en/latest/content/overview.html).

6. GFF to BAM

   `Bedtools` can help convert a feature file such as ones in GFF format to BAM: ```bedtools bedToBam -i input.gff -g genome_file > input_converted.bam``` where      `genome_file` is not a fasta file but a two column file with list of chromosomes and the corresponding chromosome sizes in basepairs.

   Please refer to BED to BAM for more details on `genome_file` format.

7. VCF to BAM

   `Bedtools` can help convert a feature file such as ones in GFF format to BAM: ```bedtools bedToBam -i input.vcf -g genome_file > input_converted.bam``` where      `genome_file` is not a fasta file but a two column file with list of chromosomes and the corresponding chromosome sizes in basepairs.


8. BAM to BED/BEDPE formats

   Please note that currently Hummingbird does not natively support BEDPE format. However, users can skip the downsampling step (please check Downsample option in    Hummingbird for more details) and continue using other features of Hummingbird.

    a) An input BAM file can be converted to a BED file (BED6 format by default) using `bedtools`: ```bedtools bamtobed -i input.bam > output.bed```

      For further details, please see https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html.

    b) An input BAM file can be converted to a BEDPE file using `bedtools`: ```bedtools bamtobed -i input.bam -bedpe > output.bedpe```

    c) For non-lossy conversions of BAM to BED can be performed using `bam2bed`functionality in the `BEDOPS` suite of tools: ```bam2bed --keep-header < input.bam        > output.bed```
  
      The `--keep-header` option is needed for the header information to be included in the output file. 
  
      In some cases, sorting and indexing of the input.bam file may be required for subsequent analyses.

9. VCF to BED/BEDPE formats

   A number of structural variant based tools help with conversion from VCF to BED or BEDPE formats such as 

    a) `lumpysv` (https://github.com/arq5x/lumpy-sv): Please refer to the script `vcfToBedpe.py`.
    
    b) `SURVIVOR`(https://github.com/fritzsedlazeck/SURVIVOR): Please refer to the `bedpetovcf` functionality.
    
    c) `svtools`(https://github.com/hall-lab/svtools): Please refer to the `bedpetovcf` and `vcftobedpe` subcommands. The benchmarking details on these                  subcommands can be found in Table 3 of their publication (https://academic.oup.com/bioinformatics/article/35/22/4782/5520944) which gives an idea of the          computational resources required and execution times.

10. BEDPE to BED12 format

    The subcommand `bedpetobed12` within the `svtools` (https://github.com/hall-lab/svtools) can convert a BEDPE file to a BED12 format.

11. BAM to FASTQ

     a) Using the `bamtofastq` functionality in `bedtools`, for paired-end data: ```bedtools bamtofastq -i aln.qsort.bam -fq aln.R1.fq -fq2 aln.R2.fq```
  
       The input bam file has to be sorted by query name that can be done using ```samtools sort -n aln.bam aln.qsort```
       In case of single-end reads, conversion is done using: ```bedtools bamtofastq  -i aln.bam -fq aln.fq```
       Further information on the various options that can be used in the `bamtofastq` command line, please see                                                          https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html

     b) Using `samtools`, for paired-end data: ```samtools fastq -1 paired1.fq -2 paired2.fq -0 /dev/null -s /dev/null -n in_sorted.bam```
        The input bam file has to be sorted before providing it for the conversion similar to the explanation in 11 (a) above. 
        For more details on the options, please refer to http://www.htslib.org/doc/samtools-fasta.html.

    Alternative tools for conversion of BAM to FASTQ can be found here: https://sites.google.com/site/wiki4metagenomics/tools/samtools/converting-bam-to-fastq.


For other file format conversions not listed here refer to,

  a) GALAXY suite of tools can be used. Please refer to "Convert Formats" in https://vclv99-241.hpc.ncsu.edu/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fvipints%2Ffml_gff3togtf%2Ffml_bed2gff%2F2.1.0&version=2.1.0&__identifer=ibvartqtce.
 
  b) Jvarkit: Java utilities for Bioinfomatics can be used. Please see http://lindenb.github.io/jvarkit/. 


Any of the above functionalities from various format conversion tools can be incoporated to the Hummingbird code by building a docker image of the tool and adding the required command lines for an automated execution.


Please note that sorting and generating the index file after the conversion may be necessary for subsequent analysis using software tools or bioinformatics pipelines.
