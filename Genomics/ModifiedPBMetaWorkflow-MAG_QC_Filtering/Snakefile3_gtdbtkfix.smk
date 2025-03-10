import os


CWD = os.getcwd()

configfile: "config.yaml"
#SAMPLES = "Nitrifying"
SAMPLES = config['samplenames']
CHECKM_DB = config["checkm_database_path"]

rule all:
    input:
        expand(os.path.join(CWD,"6-checkm2","{sample}","checkm2","quality_report.tsv"), sample = SAMPLES),
        expand(os.path.join(CWD,"6-checkm2","{sample}","{sample}.BinCount.txt"), sample = SAMPLES),
        expand(os.path.join(CWD,"8-summary","{sample}","{sample}.Complete.txt"), sample = SAMPLES)

rule Checkm2BinAnalysis:
    input:
        db = CHECKM_DB,
    output:
        qv = os.path.join(CWD, "6-checkm2", "{sample}", "checkm2", "quality_report.tsv")
    conda:
        "envs/checkm2.yml"
    threads:
        config["checkm2"]["threads"]
    params:
        indir = os.path.join(CWD, "5-dereplicated-bins", "{sample}", ""),
        outdir = os.path.join(CWD, "6-checkm2", "{sample}", "checkm2", ""),
        tmp = config["tmpdir"]
    log:
        os.path.join(CWD, "logs", "{sample}.Checkm2BinAnalysis.log")
    benchmark:
        os.path.join(CWD, "benchmarks", "{sample}.Checkm2BinAnalysis.tsv")
    shell:
        "checkm2 predict -i {params.indir} -o {params.outdir} -x fa -t {threads} --force "
        "--remove_intermediates --database_path {input.db} --tmpdir {params.tmp} &> {log}"

# Checkpoint 2 - Ensure there are bins after CheckM2, before running GTDB-Tk and the summary
checkpoint AssessCheckm2Bins:
    input:
        qv = os.path.join(CWD, "6-checkm2", "{sample}", "checkm2", "quality_report.tsv"),
    output:
        gtdb = os.path.join(CWD, "6-checkm2", "{sample}", "{sample}.GTDBTk_batch_file.txt"),
        target = os.path.join(CWD,"6-checkm2","{sample}","{sample}.BinCount.txt"),
        output_tsv = os.path.join(CWD, "6-checkm2", "{sample}", "{sample}.quality_report.tsv")
    conda:
        "envs/python.yml"
    threads:
        1
    params:
        bin_dir = os.path.join(CWD, "5-dereplicated-bins", "{sample}", ""),
        completeness = config['filters']['min_completeness'],
        contamination = config['filters']['max_contamination'],
        contigs = config['filters']['max_contigs']
    log:
        os.path.join(CWD, "logs", "{sample}.AssessCheckm2Bins.log")
    shell:
        "python scripts/Filter-Checkm2-Bins.py -i {input.qv} -b {params.bin_dir} "
        "-c1 {params.completeness} -c2 {params.contamination} -c3 {params.contigs} -o {output.gtdb} "
        "-t {output.target} -u {output.output_tsv} &> {log}"

# Function to get targets of forked workflow
def get_post_checkm2_inputs(wildcards):
    # decision is based on whether bin count file is > 0 (bins passed checkm2 filters) or not (no bins passed filtering)
    with checkpoints.AssessCheckm2Bins.get(sample=wildcards.sample).output[1].open() as fh:
        if int(fh.read().strip()) > 0:
            return os.path.join(CWD, "8-summary", "{sample}", "{sample}.Completeness-Contamination-Contigs.pdf")
        else:
            return os.path.join(CWD,"8-summary","{sample}","{sample}.No_MAGs.summary.txt")

# Checkpoint 2 aggregator; close the fork; outputs '/8-summary/SAMPLE/SAMPLE.Complete.txt'
rule CloseCheckm2Fork:
    input:
        get_post_checkm2_inputs
    output:
        os.path.join(CWD,"8-summary","{sample}","{sample}.Complete.txt")
    shell:
        "touch {output}"

##############################
# Checkpoint 2: Fork 1 - No bins passed the filtering criteria.
rule SkipGTDBAnalysis:
    input:
        gtdb = os.path.join(CWD, "6-checkm2", "{sample}", "{sample}.GTDBTk_batch_file.txt"),
        target = os.path.join(CWD,"6-checkm2","{sample}","{sample}.BinCount.txt"),
        output_tsv = os.path.join(CWD,"6-checkm2","{sample}","{sample}.quality_report.tsv")
    output:
        os.path.join(CWD,"8-summary","{sample}","{sample}.No_MAGs.summary.txt")
    threads:
        1
    params:
        outdir = os.path.join(CWD, "8-summary", "{sample}", "")
    shell:
        "mkdir -p {params.outdir} && echo No bins passed filtering in CheckM2, "
        "see {input.output_tsv} for more information > {output}"

##############################
# Checkpoint 2: Fork 2 - Bins passed filters; GTDBTkAnalysis -> GTDBTkCleanup -> MAGSummary -> MAGCopy -> MAGPlots
rule GTDBTkAnalysis:
    input:
        gtdb = os.path.join(CWD, "6-checkm2", "{sample}", "{sample}.GTDBTk_batch_file.txt"),
        target = os.path.join(CWD,"6-checkm2","{sample}","{sample}.BinCount.txt")
    output:
        dir_classify = directory(os.path.join(CWD, "7-gtdbtk", "{sample}", "classify", "")),
        complete = os.path.join(CWD,"7-gtdbtk","{sample}","{sample}.Complete.txt")
    conda:
        "envs/gtdbtk.yml"
    threads:
        config['gtdbtk']['threads']
    params:
        gtdbtk_data = config['gtdbtk']['gtdbtk_data'],
        outdir = os.path.join(CWD, "7-gtdbtk", "{sample}", "")
    log:
        os.path.join(CWD, "logs", "{sample}.GTDBTkAnalysis.log")
    benchmark:
        os.path.join(CWD, "benchmarks", "{sample}.GTDBTkAnalysis.tsv")
    shell:
        "GTDBTK_DATA_PATH={params.gtdbtk_data:q} gtdbtk classify_wf --batchfile {input.gtdb} "
        "--out_dir {params.outdir} -x fa --prefix {wildcards.sample} --mash_db gtdb-tk-r220.msh --cpus {threads} --scratch_dir /tmp "
        " &> {log} && touch {output.complete}"

# Checkpoint 2: Fork 2 - Bins passed filters; GTDBTkAnalysis -> GTDBTkCleanup -> MAGSummary -> MAGCopy -> MAGPlots
rule GTDBTkCleanup:
    input:
        dir_classify = os.path.join(CWD, "7-gtdbtk", "{sample}", "classify", ""),
        complete = os.path.join(CWD,"7-gtdbtk","{sample}","{sample}.Complete.txt")
    output:
        os.path.join(CWD,"7-gtdbtk","{sample}","{sample}.GTDBTk_Summary.txt")
    conda:
        "envs/python.yml"
    threads:
        1
    log:
        os.path.join(CWD, "logs", "{sample}.GTDBTkCleanup.log")
    shell:
        "python scripts/GTDBTk-Organize.py -i {input.dir_classify} -o {output} &> {log}"

# Checkpoint 2: Fork 2 - Bins passed filters; GTDBTkAnalysis -> GTDBTkCleanup -> MAGSummary -> MAGCopy -> MAGPlots
rule MAGSummary:
    input:
        gtdbtk = os.path.join(CWD,"7-gtdbtk","{sample}","{sample}.GTDBTk_Summary.txt"),
        checkm2 = os.path.join(CWD, "6-checkm2", "{sample}", "{sample}.quality_report.tsv")
    output:
        os.path.join(CWD,"8-summary","{sample}","{sample}.HiFi_MAG.summary.txt")
    conda:
        "envs/python.yml"
    threads:
        1
    log:
        os.path.join(CWD, "logs", "{sample}.MAGSummary.log")
    shell:
        "python scripts/MAG-Summary.py -g {input.gtdbtk} -c {input.checkm2} -o {output} &> {log}"

# Checkpoint 2: Fork 2 - Bins passed filters; GTDBTkAnalysis -> GTDBTkCleanup -> MAGSummary -> MAGCopy -> MAGPlots
rule MAGCopy:
    input:
        mag_sum = os.path.join(CWD,"8-summary","{sample}","{sample}.HiFi_MAG.summary.txt"),
        mag_dir = os.path.join(CWD, "5-dereplicated-bins", "{sample}", "")
    output:
        directory(os.path.join(CWD,"8-summary","{sample}", "MAGs", ""))
    conda:
        "envs/python.yml"
    threads:
        1
    log:
        os.path.join(CWD, "logs", "{sample}.MAGCopy.log")
    shell:
        "python scripts/Copy-Final-MAGs.py -i {input.mag_sum} -m {input.mag_dir} -o {output} &> {log}"

# Checkpoint 2: Fork 2 - Bins passed filters; GTDBTkAnalysis -> GTDBTkCleanup -> MAGSummary -> MAGCopy -> MAGPlots
rule MAGPlots:
    input:
        checkm2_tsv = os.path.join(CWD, "6-checkm2", "{sample}", "{sample}.quality_report.tsv"),
        mag_dir = os.path.join(CWD,"8-summary","{sample}", "MAGs", ""),
        mag_eval = os.path.join(CWD,"8-summary","{sample}","{sample}.HiFi_MAG.summary.txt")
    output:
        o1 = os.path.join(CWD, "8-summary", "{sample}", "{sample}.All-DASTool-Bins.pdf"),
        o2 = os.path.join(CWD, "8-summary", "{sample}", "{sample}.Completeness-Contamination-Contigs.pdf")
    conda:
        "envs/python.yml"
    threads:
        1
    params:
        completeness = config['filters']['min_completeness'],
        contamination = config['filters']['max_contamination']
    log:
        os.path.join(CWD, "logs", "{sample}.MAGPlots.log")
    shell:
        "python scripts/Plot-Figures.py -i1 {input.checkm2_tsv} -i2 {input.mag_eval} -l {wildcards.sample} "
        "-c1 {params.completeness} -c2 {params.contamination} -o1 {output.o1} -o2 {output.o2}"
        "&> {log}"
