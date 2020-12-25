## Section 7: Alternative Downsampling Methods

For users interested in downsampling techniques other than the ones supported by Hummingbird, please refer to the examples below:

1. Downsampling SAM/BAM files using `DownsampleSam` tool from Picard 

   NOTE: Implemented in Hummingbird but is fixed in terms of options used within the Picard tool as in the following command line ```java -jar picard.jar            DownsampleSam I=input.bam O=downsampled.bam STRATEGY=Chained P=0.02 ACCURACY=0.0001``` which retains only 2% of the reads in the input file and this percentage    comes from the downsampling fraction(s) provided by the Hummingbird user.

   For a better accuracy when dealing with smaller fractions such as retaining 0.001% of the reads, one can use ```java -jar picard.jar DownsampleSam I=input.bam    O=downsampled.bam STRATEGY=HighAccuracy P=0.00001 ACCURACY=0.0000001```

   This tool offers a number of strategies for downsampling as well as levels of accuracy (combinations of which are not offered by Hummingbird currently) which      can be dependent on memory availability.
   For more options that can be used with DownsampleSam, please see https://gatk.broadinstitute.org/hc/en-us/articles/360036431292-DownsampleSam-Picard-.

2. Downsampling SAM/BAM files at the chromosome level using `samtools`

   For example, to extract chromosome 22 from a bam file and obtain a bam file with only chr 22: ```samtools view -b -o <output.bam> -@<INT_threads> <input.bam>      chr22``` where `threads` is an optional parameter to speed up the process of extraction. Regions within a specific chromosome can also be specified for            extraction.

   The user needs to ensure that the `input.bam` file is sorted (```samtools sort example.bam -o example_sorted.bam```), indexed (```samtools index                  example_sorted.bam```) and chromosome specified matches the chromosome name in `input.bam`. Further details on samtools view usage can be found here:              http://www.htslib.org/doc/samtools-view.html.

3. Downsampling BAM files at the chromosomal region level using ENSEMBL tool `Data Slicer` 

   Another way downsampling BAM files is via the `Data Slicer` tool (http://grch37.ensembl.org/Homo_sapiens/Tools/DataSlicer) available as part of the ENSEMBLE      suite of tools that provides a GUI for the users. The downsampling of the BAM is done based on chromosome and coordinates provided by the user. Details on the    usage can be found here grch37.ensembl.org/Help/View?id=575.

4. Downsampling VCF files

    a) Using the `DownSampleVcf` function in the `jvarkit` tool (http://lindenb.github.io/jvarkit/DownSampleVcf.html), a VCF file can downsampled by specifying          the number of random variants to be extracted: ```curl -skL "ftp://ftp-                                                                                            trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz" | gunzip -c |java -jar                  dist/downsamplevcf.jar -n 100 > output.vcf```

    b) Using the `Data Slicer` tool (http://grch37.ensembl.org/Homo_sapiens/Tools/DataSlicer) from the ENSEMBLE project which provides a GUI for the users and            subsamples based on chromosome and coordinates provided by the user. For further help on usage, please refer to grch37.ensembl.org/Help/View?id=575.


NOTE: Some of the file formats other than BAM or fastq/fastq.gz if provided in the gunzip compressed format, can be downsampled by the `zless` functionality in Hummingbird as long as the total number of lines in the original input file is provided using the `target` flag in the `Downsample` option.

The above downsampled files can be provided to Hummingbird to run the Memory Profiler step and then receive the recommended instance types from the Recommendation Engine. 

The examples provided in this documentation are by no means an exhaustive list but just a guide for users to consider when running different bioinformatics pipelines that accept varied input file formats. We will be adding support to different input file formats in the downsampling step of Hummingbird in the future.
