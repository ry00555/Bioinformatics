#!/bin/bash
#SBATCH --job-name=j_GATK
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ry00555@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=GATK.%j.out
#SBATCH --error=GATK.%j.err

cd $SLURM_SUBMIT_DIR
Working Directory = /home/ry00555/Bioinformatics/CrassaGenome

#Load these modules that are compatible with GATK version 4.3
ml GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8
ml picard/2.27.4-Java-13.0.2
ml BWA/0.7.17-GCC-8.3.0
module load BEDTools/2.29.2-GCC-8.3.0
module load SAMtools/1.10-GCC-8.3.0
module load BEDOPS/2.4.39-foss-2019b
ml Bowtie2/2.4.1-GCC-8.3.0
ml R/3.6.2-foss-2019b

#This section will create the reference files needed for most GATK tool commands.
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
  -R GCF_000182925.2.fasta \
  -O GCF_000182925.2.dict

java -jar $EBROOTPICARD/picard.jar BedToIntervalList \
-I crassa.bed \
-R GCF_000182925.2.fasta \
-SD GCF_000182925.2.dict \
-O Crassa.interval_list


gatk PreprocessIntervals \
          -R GCF_000182925.2.fasta \
          -L Crassa.interval_list \
          --interval-merging-rule OVERLAPPING_ONLY \
          --bin-length 1000 \
          --padding 0 \
          -O Crassa.preprocessed_intervals.interval_list

gatk AnnotateIntervals \
  -R GCF_000182925.2.fasta \
  -L Crassa.preprocessed_intervals.interval_list \
  --interval-merging-rule OVERLAPPING_ONLY \
  -O Crassa_preprocessed_annotated_intervals.tsv

#109_58 is a wildtype sample. 109_59 is a ∆mus30, ∆mei3 sample. The bam files were aligned using MapCutandRun.sh using the reference fasta file, GCF_000182925.2.fasta, to keep consistent. The sorted bam files are in  /scratch/ry00555/OutputRun109/Run109Bam
#There is a way to create bam files from pair end reads fasta.gz files, but when mapped to IGV, they are empty.
#In MapCutandRun.sh, read group ID, sample name, library number are not annotated in the bam file. I ran GATK AddOrReplaceReadGroups after samtools index 109_nn_Genomic.bam files.
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
                -I /scratch/ry00555/OutputRun109/Run109Bam/109_59_Genomic.bam \
                -O Run109ReadGroups/109_59_Genomic.bamoutput.bam \
                -RGID 1 \
                -RGLB lib1 \
                -RGPL illumina \
                -RGPU S34 \
                -RGSM 109_59

                java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
                                -I /scratch/ry00555/OutputRun109/Run109Bam/109_58_Genomic.bam \
                                -O Run109ReadGroups/109_58_Genomic.bamoutput.bam \
                                -RGID 2 \
                                -RGLB lib1 \
                                -RGPL illumina \
                                -RGPU S34 \
                                -RGSM 109_58


                                ry00555@ss-sub2 CrassaGenome$ gatk CollectReadCounts \
                                >        -I Run109ReadGroups/109_58_Genomic.bamoutput.bam \
                                >        -R GCF_000182925.2.fasta \
                                >        -L Crassa.interval_list \
                                >        --interval-merging-rule OVERLAPPING_ONLY \
                                >        --disable-sequence-dictionary-validation false \
                                >        --sequence-dictionary GCF_000182925.2.dict \
                                >        -O 109tsv/109_58.counts.tsv
                                Using GATK jar /apps/eb/GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8/gatk-package-4.3.0.0-local.jar
                                Running:
                                    java -Dsamjdk.use_async_io_read_samtools=false -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Dsamjdk.compression_level=2 -jar /apps/eb/GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8/gatk-package-4.3.0.0-local.jar CollectReadCounts -I Run109ReadGroups/109_58_Genomic.bamoutput.bam -R GCF_000182925.2.fasta -L Crassa.interval_list --interval-merging-rule OVERLAPPING_ONLY --disable-sequence-dictionary-validation false --sequence-dictionary GCF_000182925.2.dict -O 109tsv/109_58.counts.tsv
                                16:04:37.605 INFO  NativeLibraryLoader - Loading libgkl_compression.so from jar:file:/apps/eb/GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8/gatk-package-4.3.0.0-local.jar!/com/intel/gkl/native/libgkl_compression.so
                                16:04:37.719 INFO  CollectReadCounts - ------------------------------------------------------------
                                16:04:37.720 INFO  CollectReadCounts - The Genome Analysis Toolkit (GATK) v4.3.0.0
                                16:04:37.720 INFO  CollectReadCounts - For support and documentation go to https://software.broadinstitute.org/gatk/
                                16:04:37.720 INFO  CollectReadCounts - Executing as ry00555@ss-sub2.gacrc.uga.edu on Linux v3.10.0-1160.36.2.el7.x86_64 amd64
                                16:04:37.720 INFO  CollectReadCounts - Java runtime: OpenJDK 64-Bit Server VM v13.0.2+8
                                16:04:37.720 INFO  CollectReadCounts - Start Date/Time: November 30, 2022 at 4:04:37 PM EST
                                16:04:37.720 INFO  CollectReadCounts - ------------------------------------------------------------
                                16:04:37.720 INFO  CollectReadCounts - ------------------------------------------------------------
                                16:04:37.721 INFO  CollectReadCounts - HTSJDK Version: 3.0.1
                                16:04:37.721 INFO  CollectReadCounts - Picard Version: 2.27.5
                                16:04:37.721 INFO  CollectReadCounts - Built for Spark Version: 2.4.5
                                16:04:37.721 INFO  CollectReadCounts - HTSJDK Defaults.COMPRESSION_LEVEL : 2
                                16:04:37.721 INFO  CollectReadCounts - HTSJDK Defaults.USE_ASYNC_IO_READ_FOR_SAMTOOLS : false
                                16:04:37.721 INFO  CollectReadCounts - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_SAMTOOLS : true
                                16:04:37.721 INFO  CollectReadCounts - HTSJDK Defaults.USE_ASYNC_IO_WRITE_FOR_TRIBBLE : false
                                16:04:37.722 INFO  CollectReadCounts - Deflater: IntelDeflater
                                16:04:37.722 INFO  CollectReadCounts - Inflater: IntelInflater
                                16:04:37.722 INFO  CollectReadCounts - GCS max retries/reopens: 20
                                16:04:37.722 INFO  CollectReadCounts - Requester pays: disabled
                                16:04:37.722 INFO  CollectReadCounts - Initializing engine
                                WARNING: BAM index file /home/ry00555/Bioinformatics/CrassaGenome/Run109ReadGroups/109_58_Genomic.bamoutput.bam.bai is older than BAM /home/ry00555/Bioinformatics/CrassaGenome/Run109ReadGroups/109_58_Genomic.bamoutput.bam
                                16:04:37.905 INFO  FeatureManager - Using codec IntervalListCodec to read file file:///home/ry00555/Bioinformatics/CrassaGenome/Crassa.interval_list
                                16:04:37.915 INFO  IntervalArgumentCollection - Processing 40463072 bp from intervals
                                16:04:37.919 INFO  CollectReadCounts - Done initializing engine
                                16:04:37.921 INFO  CollectReadCounts - Collecting read counts...
                                16:04:37.921 INFO  ProgressMeter - Starting traversal
                                16:04:37.921 INFO  ProgressMeter -        Current Locus  Elapsed Minutes       Reads Processed     Reads/Minute
                                16:04:47.922 INFO  ProgressMeter -         chr5:5758647              0.2               6050000       36296370.4
                                16:04:50.643 INFO  CollectReadCounts - 0 read(s) filtered by: WellformedReadFilter
                                0 read(s) filtered by: MappedReadFilter
                                0 read(s) filtered by: NonZeroReferenceLengthAlignmentReadFilter
                                0 read(s) filtered by: NotDuplicateReadFilter
                                117377 read(s) filtered by: MappingQualityReadFilter
                                117377 total reads filtered
                                16:04:50.644 INFO  ProgressMeter -         chr7:4249388              0.2               7749541       36548692.0
                                16:04:50.644 INFO  ProgressMeter - Traversal complete. Processed 7749541 total reads in 0.2 minutes.
                                16:04:50.644 INFO  CollectReadCounts - Writing read counts to /home/ry00555/Bioinformatics/CrassaGenome/109tsv/109_58.counts.tsv...
                                log4j:WARN No appenders could be found for logger (org.broadinstitute.hdf5.HDF5Library).
                                log4j:WARN Please initialize the log4j system properly.
                                log4j:WARN See http://logging.apache.org/log4j/1.2/faq.html#noconfig for more info.
                                16:04:50.719 INFO  CollectReadCounts - CollectReadCounts complete.
                                16:04:50.719 INFO  CollectReadCounts - Shutting down engine
                                [November 30, 2022 at 4:04:50 PM EST] org.broadinstitute.hellbender.tools.copynumber.CollectReadCounts done. Elapsed time: 0.22 minutes.
                                Runtime.totalMemory()=236978176
                                ry00555@ss-sub2 CrassaGenome$

  gatk CollectReadCounts \
       -I Run109ReadGroups/109_58_Genomic.bamoutput.bam \
       -R GCF_000182925.2.fasta \
       -L Crassa.preprocessed_intervals.interval_list  \
       --interval-merging-rule OVERLAPPING_ONLY \
       -O 109tsv/109_58preprocessedcnv.counts.tsv


       gatk CollectReadCounts \
            -I /scratch/ry00555/OutputRun109/Run109Bam/109_60_Genomic.bam \
            -R GCF_000182925.2.fasta \
            -L Crassa.interval_list \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O 109tsv/109_60.counts.tsv

            gatk CollectReadCounts \
                 -I Run109ReadGroups/109_58_Genomic.bamoutput.bam \
                 -R GCF_000182925.2.fasta \
                 -L Crassa.interval_list \
                 --interval-merging-rule OVERLAPPING_ONLY \
                 -O 109tsv/109_59.counts.tsv


                 gatk CollectReadCounts \
                      -I /scratch/ry00555/OutputRun109/Run109Bam/109_60_Genomic.bam \
                      -R GCF_000182925.2.fasta \
                      -L Crassa.interval_list \
                      --interval-merging-rule OVERLAPPING_ONLY \
                      -O 109tsv/109_60.counts.tsv

            java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
                            -I /scratch/ry00555/OutputRun109/Run109Bam/109_60_Genomic.bam \
                            -O Run109ReadGroups/109_60_Genomic.bamoutput.bam \
                            -RGID 3 \
                            -RGLB lib1 \
                            -RGPL illumina \
                            -RGPU S34 \
                            -RGSM 109_60
                            gatk CollectReadCounts \
                                 -I Run109ReadGroups/109_60_Genomic.bamoutput.bam  \
                                 -R GCF_000182925.2.fasta \
                                 -L Crassa.interval_list \
                                 --interval-merging-rule OVERLAPPING_ONLY \
                                 -O 109tsv/109_60.counts.tsv

                                 java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
                                                 -I /scratch/ry00555/OutputRun109/Run109Bam/109_63_Genomic.bam \
                                                 -O Run109ReadGroups/109_63_Genomic.bamoutput.bam \
                                                 -RGID 6 \
                                                 -RGLB lib1 \
                                                 -RGPL illumina \
                                                 -RGPU S34 \
                                                 -RGSM 109_63
                                                 gatk CollectReadCounts \
                                                      -I Run109ReadGroups/109_63_Genomic.bamoutput.bam  \
                                                      -R GCF_000182925.2.fasta \
                                                      -L Crassa.preprocessed_intervals.interval_list \
                                                      --interval-merging-rule OVERLAPPING_ONLY \
                                                      -O 109tsv/109_63preprocessed.counts.tsv

#ERROR occurs with the above command: java.lang.IllegalArgumentException: Reference name for '150' not found in sequence dictionary.
#	at htsjdk.samtools.SAMRecord.resolveNameFromIndex(SAMRecord.java:579)
#	at htsjdk.samtools.SAMRecord.setReferenceIndex(SAMRecord.java:432)
#	at htsjdk.samtools.BAMRecord.<init>(BAMRecord.java:110)
#	at htsjdk.samtools.DefaultSAMRecordFactory.createBAMRecord(DefaultSAMRecordFactory.java:42)
#	at htsjdk.samtools.BAMRecordCodec.decode(BAMRecordCodec.java:283)
#	at htsjdk.samtools.BAMFileReader$BAMFileIterator.getNextRecord(BAMFileReader.java:880)
#	at htsjdk.samtools.BAMFileReader$BAMFileIndexIterator.getNextRecord(BAMFileReader.java:1019)
#	at htsjdk.samtools.BAMFileReader$BAMFileIterator.advance(BAMFileReader.java:854)
#	at htsjdk.samtools.BAMFileReader$BAMFileIndexIterator.<init>(BAMFileReader.java:1001)
#	at htsjdk.samtools.BAMFileReader.createIndexIterator(BAMFileReader.java:948)
#	at htsjdk.samtools.BAMFileReader.query(BAMFileReader.java:626)
#	at htsjdk.samtools.SamReader$PrimitiveSamReaderToSamReaderAdapter.query(SamReader.java:550)
#	at htsjdk.samtools.SamReader$PrimitiveSamReaderToSamReaderAdapter.queryOverlapping(SamReader.java:417)
#	at org.broadinstitute.hellbender.utils.iterators.SamReaderQueryingIterator.loadNextIterator(SamReaderQueryingIterator.java:130)
#	at org.broadinstitute.hellbender.utils.iterators.SamReaderQueryingIterator.<init>(SamReaderQueryingIterator.java:69)
#	at org.broadinstitute.hellbender.engine.ReadsPathDataSource.prepareIteratorsForTraversal(ReadsPathDataSource.java:412)
#	at org.broadinstitute.hellbender.engine.ReadsPathDataSource.iterator(ReadsPathDataSource.java:336)
#	at java.base/java.lang.Iterable.spliterator(Iterable.java:101)
#	at org.broadinstitute.hellbender.utils.Utils.stream(Utils.java:1176)
#	at org.broadinstitute.hellbender.engine.GATKTool.getTransformedReadStream(GATKTool.java:384)
#	at org.broadinstitute.hellbender.engine.ReadWalker.traverse(ReadWalker.java:97)
#	at org.broadinstitute.hellbender.engine.GATKTool.doWork(GATKTool.java:1095)
#	at org.broadinstitute.hellbender.cmdline.CommandLineProgram.instanceMainPostParseArgs(CommandLineProgram.java:192)
#	at org.broadinstitute.hellbender.cmdline.CommandLineProgram.instanceMain(CommandLineProgram.java:211)
#	at org.broadinstitute.hellbender.Main.runCommandLineProgram(Main.java:160)
#	at org.broadinstitute.hellbender.Main.mainEntry(Main.java:203)
#	at org.broadinstitute.hellbender.Main.main(Main.java:289)


       gatk CreateReadCountPanelOfNormals \
         -I 109tsv/109_58preprocessedcnv.counts.tsv \
         --annotated-intervals Crassa_preprocessed_annotated_intervals.tsv \
         -O PanelofNormals/109_58preprocessedcnv.pon.hdf5

         gatk CollectReadCounts \
              -I Run109ReadGroups/109_63_Genomic.bamoutput.bam  \
              -R GCF_000182925.2.fasta \
              -L Crassa.preprocessed_intervals.interval_list \
              --interval-merging-rule OVERLAPPING_ONLY \
              -O 109tsv/109_63preprocessed.counts.tsv

         gatk DenoiseReadCounts \
                   -I 109tsv/109_63preprocessed.counts.tsv \
                   --annotated-intervals Crassa_annotated_intervals.tsv \
                   --count-panel-of-normals PanelofNormals/109_58preprocessedcnv.pon.hdf5 \
                   --standardized-copy-ratios CopyRatios/109_63preprocessed.standardizedCR.tsv \
                   --denoised-copy-ratios CopyRatios/109_63preprocessed.denoisedCR.tsv

         gatk PlotDenoisedCopyRatios \
                         --standardized-copy-ratios CopyRatios/109_63preprocessed.standardizedCR.tsv \
                         --denoised-copy-ratios CopyRatios/109_63preprocessed.denoisedCR.tsv \
                         --sequence-dictionary GCF_000182925.2.dict \
                         --output-prefix Run109CNV_63preprocessed \
                         --output PlotDenoisedCopyRatios

scp -r ry00555@xfer.gacrc.uga.edu:/home/ry00555/Bioinformatics/CrassaGenome/PlotDenoisedCopyRatios/Run109CNV_63preprocessed.denoised.png /Users/ry00555/Desktop/NeurosporaGenome


 gatk ModelSegments \
           --denoised-copy-ratios CopyRatios/109_59.denoisedCR.tsv \
           --output-prefix 109_59 \
           -O ModelSegments


 gatk PlotModeledSegments \
                --denoised-copy-ratios CopyRatios/109_59.denoisedCR.tsv \
          --segments ModelSegments/109_59.modelFinal.seg \
         --sequence-dictionary GCF_000182925.2.dict \
         --output-prefix 109_59 \
         -O PlotModelSegments
