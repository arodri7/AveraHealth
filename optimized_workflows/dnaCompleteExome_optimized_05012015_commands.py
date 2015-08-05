# -*- coding: utf-8 -*-
#!/usr/bin/env python

"""
Runs the Meissnert SNPiR RNA-variant discovery pipeline in optimized mode in one Node
See below for options
"""

def create_wf_meta(options, output_dir, tmp_dir):

    # reads
    fastq = options.fastq
    if options.rfastq:
        rfastq = options.rfastq
    read_length = get_readLength(fastq)

    # get sample name
    sample_name = get_sampleName(fastq)

    # get ncores
    ncpu = get_ncores()
    nthreads = ncpu/2

    # get memory
    jvm_heap = get_linuxRAM()

    # get the reference datasets
    reference_db = options.gatk_ref
    dbsnp="/mnt/galaxyIndices/genomes/Hsapiens/human_g1k_v37/annotation/dbsnp_138.b37.vcf"
    mills="/mnt/galaxyIndices/genomes/Hsapiens/human_g1k_v37/annotation/Mills_and_1000G_gold_standard.indels.b37.vcf"
    mask="/mnt/galaxyIndices/genomes/Hsapiens/human_g1k_v37/annotation/b37_cosmic_v54_120711.vcf"
    #GATK3_PATH = os.environ['GATK3_PATH']
    GATK3_PATH = "/mnt/galaxyTools/tools/gatk3/GenomeAnalysisTK-3.4-0"
    #JAVA_JAR_PATH = os.environ['JAVA_JAR_PATH']
    JAVA_JAR_PATH = "/mnt/galaxyTools/tools/picard/1.134"

    bwa_index = options.bwa_ref
    ## Setup workflow
    command_meta = {}

    read_group = "--rgid=\"%s\" --rglb=\"%s\" --rgpl=\"%s\" --rgsm=\"%s\"" % (options.rgid, options.rglb, options.rgpl, options.rgsm)
    if options.rgcn:
        read_group += " --rgcn=\"%s\"" % options.rgcn
    if options.rgds:
        read_group += " --rgds=\"%s\"" % options.rgds
    if options.rgdt:
        read_group += " --rgdt=\"%s\"" % options.rgdt
    if options.rgfo:
        read_group += " --rgfo=\"%s\"" % options.rgfo
    if options.rgks:
        read_group += " --rgks=\"%s\"" % options.rgks
    if options.rgpg:
        read_group += " --rgpg=\"%s\"" % options.rgpg
    if options.rgpi:
        read_group += " --rgpi=\"%s\"" % options.rgpi
    if options.rgpu:
        read_group += " --rgpu=\"%s\"" % options.rgpu

    #
    step = 0
    ncpu_step = int(ncpu) - 2
    command_meta[step] = []
    input1 = fastq
    input2 = rfastq
    input_files = [input1, input2]
    outputF = "%s/%s-bwa-out.bam" % (output_dir, sample_name)
    output_files = [outputF]
    cmd = " python /opt/galaxy/tools/sr_mapping/bwa_mem.py --threads=%s --fileSource=indexed --ref=%s --fastq=%s --rfastq=%s --output=%s --genAlignType=paired --params=full --minSeedLength 19 --bandWidth 100 --offDiagonal 100 --internalSeeds 1.5 --seedsOccurrence 10000 --seqMatch 1 --mismatch 4 --gapOpen 6 --gapExtension 1 --clipping 5 --unpairedReadpair 17 --minScore 30 %s -M --bam" % (ncpu_step, bwa_index, input1, input2, outputF, read_group)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files": output_files})

    inputF = rfastq
    fastqc_dir = "%s/%s-FastQCr2_dir" % (output_dir, sample_name)
    name = "FastQCr2"
    input_files = [inputF]
    outputF = "%s/%s-fastqc_out_2.html" % (output_dir, sample_name)
    output_files = [outputF, fastqc_dir]
    cmd = "python /opt/galaxy/tools/rgenetics/rgFastQC.py -i %s -d %s -o %s -n %s -f fastqsanger -j %s -e /nfs/software/galaxy/tool-data/shared/jars/FastQC/fastqc" % (inputF, fastqc_dir, outputF, name, name)
    command_meta[step].append({"cl":cmd, "input_files": input_files, "output_file":outputF, "output_files":output_files})

    inputF = fastq
    fastqc_dir = "%s/%s-FastQCr1_dir" % (output_dir, sample_name)
    name = "FastQCr1"
    input_files = [inputF]
    outputF = "%s/%s-fastqc_out_1.html" % (output_dir, sample_name)
    output_files = [outputF, fastqc_dir]
    cmd = "python /opt/galaxy/tools/rgenetics/rgFastQC.py -i %s -d %s -o %s -n %s -f fastqsanger -j %s -e /nfs/software/galaxy/tool-data/shared/jars/FastQC/fastqc" % (inputF, fastqc_dir, outputF, name, name)
    command_meta[step].append({"cl":cmd, "input_files": input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 1
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.bam" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "python /opt/galaxy/tools/sambamba/sambamba_sort.py --memory %sM --input=%s --order=coordinate --output=%s" % (jvm_heap, inputF, outputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 2
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.bam" % (output_dir, sample_name)
    outputMetrics = "%s/%s-SM1.dups" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputMetrics]
    cmd = "python /opt/galaxy/tools/picard/picard_wrapper.py --maxjheap %sm -i %s -n Dupes_Marked --tmpdir %s -o %s --remdups true --assumesorted true --readregex \"[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\" --optdupdist 100 -j \"%s/picard.jar MarkDuplicates\" -d %s -t %s -e bam" % (jvm_heap, inputF, tmp_dir, outputF, JAVA_JAR_PATH, tmp_dir, outputMetrics)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 3
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.bam.bai" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools index %s" % (inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 4
    command_meta[step] = []
    inputF = command_meta[step-2][0]["output_file"]
    inputIndex = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-output.intervals" % (output_dir, sample_name)
    outputLog = "%s/%s-output.intervals.log" % (output_dir, sample_name)
    input_files = [inputF, inputIndex]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T RealignerTargetCreator -ip %s -o %s --num_threads %s -R %s\"  -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -d \"-known:dbsnp,%%(file_type)s\" \"%s\" \"vcf\" \"input_dbsnp_0\" -d \"-known:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_0\"" % (outputLog, inputF, inputIndex, GATK3_PATH, options.ip, outputF, ncpu, reference_db, options.bed, dbsnp, mills)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 5
    command_meta[step] = []
    inputBam = command_meta[step-3][0]["output_file"]
    inputF = command_meta[step-1][0]["output_file"]
    inputIndex = command_meta[step-2][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.realigned.bam" % (output_dir, sample_name)
    outputLog = "%s/%s-bwa-out.sort.rd.realigned.log" % (output_dir, sample_name)
    input_files = [inputF, inputIndex, inputBam]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T IndelRealigner -ip %s -o %s -R %s -LOD 5.0\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -d \"-known:dbsnp,%%(file_type)s\" \"%s\" \"vcf\" \"input_dbsnp_0\" -d \"-known:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_0\"  -d \"-targetIntervals\" \"%s\" \"gatk_interval\" \"gatk_target_intervals\"" % (outputLog, inputBam, inputIndex, GATK3_PATH, options.ip, outputF, reference_db, options.bed, dbsnp, mills, inputF )
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 6
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.realigned.bam.bai" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools index %s" % (inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 7
    command_meta[step] = []
    inputBam = command_meta[step-2][0]["output_file"]
    inputIndex = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-recal_data.table" % (output_dir, sample_name)
    outputLog = "%s/%s-bwa-out.sort.rd.realigned.log" % (output_dir, sample_name)
    input_files = [inputBam, inputIndex]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\"  \"%s\" \"bam_index\" \"gatk_input\" -p \"java -jar %s/GenomeAnalysisTK.jar -T BaseRecalibrator -ip %s -nct %s -R %s --out %s\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -d \"--knownSites:dbsnp,%%(file_type)s\" \"%s\" \"vcf\" \"input_dbsnp_0\" -d \"--knownSites:indels,%%(file_type)s\" \"%s\" \"vcf\" \"input_indels_0\"" % (outputLog, inputBam, inputIndex, GATK3_PATH, options.ip, ncpu, reference_db, outputF, options.bed, dbsnp, mills )
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 8
    command_meta[step] = []
    inputBam = command_meta[step-3][0]["output_file"]
    inputTable = command_meta[step-1][0]["output_file"]
    inputIndex = command_meta[step-2][0]["output_file"]
    #outputF = "%s/%s-bwa-out.sort.rd.realigned.recal.bam" % (output_dir, sample_name)
    outputF = options.output_bam
    outputLog = "%s/%s-bwa-out.sort.rd.realigned.recal.log" % (output_dir, sample_name)
    input_files = [inputTable, inputBam, inputIndex]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap_fraction 1 --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input\" -d \"\" \"%s\" \"bam_index\" \"gatk_input\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T PrintReads -ip %s -o %s -nct %s -R %s --BQSR %s\"" % (outputLog, inputBam, inputIndex, options.bed, GATK3_PATH, options.ip, outputF, ncpu, reference_db, inputTable)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 9
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-bwa-out.sort.rd.realigned.recal.bam.bai" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF]
    cmd = "samtools index %s" % (inputF)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 10
    command_meta[step] = []
    inputBam = command_meta[step-2][0]["output_file"]
    inputIndex = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-raw_variants.vcf" % (output_dir, sample_name)
    outputLog = "%s/%s-raw_variants.log" % (output_dir, sample_name)
    input_files = [inputBam, inputIndex]
    output_files = [outputF, outputLog, inputIndex]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap %sm --stdout %s -d \"-I\" \"%s\" \"bam\" \"gatk_input_0\"  -d \"\" \"%s\" \"bam_index\" \"gatk_input_0\" -d \"--dbsnp:dbsnp,%%(file_type)s\" \"%s\" \"vcf\" \"input_dbsnp_0\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T HaplotypeCaller -ip %s -nct %s --out %s -R %s --output_mode EMIT_VARIANTS_ONLY\"" % (jvm_heap, outputLog, inputBam, inputIndex, dbsnp, options.bed, GATK3_PATH, options.ip, ncpu, outputF, reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 11
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    outputF = "%s/%s-raw_variants.filter.vcf" % (output_dir, sample_name)
    outputLog = "%s/%s-raw_variants.filter.log" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap %sm --stdout %s -d \"--variant:variant,%%(file_type)s\" \"%s\" \"vcf\" \"input_variant\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T VariantFiltration -ip %s -o %s -R %s --maskExtension 0 --maskName Mask\" -d \"--mask:Mask,%%(file_type)s\" \"%s\" \"vcf\" \"input_mask_Mask\"" % (jvm_heap, outputLog, inputF, options.bed, GATK3_PATH, options.ip,  outputF, reference_db, mask)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    #
    step = 12
    command_meta[step] = []
    inputF = command_meta[step-1][0]["output_file"]
    #outputF = "%s/%s-raw_variants.filter.select.vcf" % (output_dir, sample_name)
    outputF = options.output_vcf
    outputLog = "%s/%s-raw_variants.filter.select.log" % (output_dir, sample_name)
    input_files = [inputF]
    output_files = [outputF, outputLog]
    cmd = "python /opt/galaxy/tools/gatk3/gatk3_wrapper.py --max_jvm_heap %sm --stdout %s -d \"--variant:variant,%%(file_type)s\" \"%s\" \"vcf\" \"input_variant\" -d \"--intervals\" \"%s\" \"bed\" \"input_intervals_0\" -p \"java -jar %s/GenomeAnalysisTK.jar -T SelectVariants -ip %s --num_threads %s -o %s -R %s\"" % (jvm_heap, outputLog, inputF, options.bed, GATK3_PATH, options.ip, ncpu, outputF, reference_db)
    command_meta[step].append({"cl":cmd, "input_files":input_files, "output_file":outputF, "output_files":output_files})

    return command_meta
