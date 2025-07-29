import pandas as pd

configfile: "config.yaml"

# Load sample info
samples = pd.read_csv("samples.tsv", sep="\t")
SAMPLE_NAMES = samples["sample"].tolist()
SAMPLE_DICT = dict(zip(samples["sample"], zip(samples["R1"], samples["R2"])))

rule all:
    input:
        expand("{outdir}/Trimmed_Reads/{sample}/{sample}_1_val_1.fq", outdir=config["dir"]["outdir"], sample=SAMPLE_NAMES) +
        expand("{outdir}/Trimmed_Reads/{sample}/{sample}_2_val_2.fq", outdir=config["dir"]["outdir"], sample=SAMPLE_NAMES) +
        expand("{outdir}/QC/{sample}_1_val_1_fastqc.html", outdir=config["dir"]["outdir"], sample=SAMPLE_NAMES) +
        expand("{outdir}/QC/{sample}_2_val_2_fastqc.html", outdir=config["dir"]["outdir"], sample=SAMPLE_NAMES) +
        expand("{outdir}/GTF_nc_transcript_filtered/{sample}.gtf",outdir=config["dir"]["outdir"], sample=SAMPLE_NAMES)


rule trimgalore:
    input:
        fq1 = lambda wc: SAMPLE_DICT[wc.sample][0],
        fq2 = lambda wc: SAMPLE_DICT[wc.sample][1]
    output:
        fq1_trimmed = config["dir"]["outdir"] + "/Trimmed_Reads/{sample}/{sample}_1_val_1.fq",
        fq2_trimmed = config["dir"]["outdir"] + "/Trimmed_Reads/{sample}/{sample}_2_val_2.fq"
    params:
        outdir = lambda wc: f"{config['dir']['outdir']}/Trimmed_Reads/{wc.sample}"
    conda:
    	"conda/trimgalore.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        trim_galore --illumina --paired -o {params.outdir} {input.fq1} {input.fq2}
        """

rule fastqc:
    input:
        fq1 = rules.trimgalore.output.fq1_trimmed,
        fq2 = rules.trimgalore.output.fq2_trimmed
    output:
        html1 = "{outdir}/QC/{sample}_1_val_1_fastqc.html",
        html2 = "{outdir}/QC/{sample}_2_val_2_fastqc.html",
        zip1 = temp("{outdir}/QC/{sample}_1_val_1_fastqc.zip"),
        zip2 = temp("{outdir}/QC/{sample}_2_val_2_fastqc.zip")
    params:
        outdir = lambda wc: f"{config['dir']['outdir']}/QC"
    conda:
        "conda/fastqc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} {input.fq1} {input.fq2}
        """

rule star_align:
    input:
        fq1 = config["dir"]["outdir"] + "/Trimmed_Reads/{sample}/{sample}_1_val_1.fq",
        fq2 = config["dir"]["outdir"] + "/Trimmed_Reads/{sample}/{sample}_2_val_2.fq"
    output:
        bam = config["dir"]["outdir"] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam"
    params:
        outdir = config["dir"]["outdir"] + "/STAR/{sample}",
        index = config["star_index"]  # Put this in config.yaml
    threads: 8
    conda:
        "conda/star.yaml"
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix {params.outdir}/ \
            --outSAMtype BAM SortedByCoordinate
        """

rule samtools_index:
    input:
        bam = config["dir"]["outdir"] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        bai = config["dir"]["outdir"] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "conda/samtools.yaml"
    shell:
        """
        samtools index {input.bam}
        """

rule stringtie:
    input:
        bam = config["dir"]["outdir"] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam",
        bai = config["dir"]["outdir"] + "/STAR/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        gtf = config["dir"]["outdir"] + "/StringTie/{sample}/transcripts.gtf"
    params:
        outdir = config["dir"]["outdir"] + "/StringTie/{sample}",
        reference_gtf = config["reference_gtf"],
        ballgown = config["dir"]["outdir"]+"/StringTie/{sample}_ballgown_input/",
        label = "{sample}-Stringtie"
    threads: 4
    conda:
        "conda/stringtie.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        stringtie {input.bam} \
            -G {params.reference_gtf} \
            -o {output.gtf} \
            -p {threads} \
            -b {params.ballgown} \
            -l {params.label}
        """

rule gffcompare_for_novel_transcripts:
	input:
		ref_gff = config['reference_gtf'],
		query_gtf = config["dir"]["outdir"]+"/StringTie/{sample}/transcripts.gtf",
	output:
		novel_transcripts_list = config["dir"]["outdir"]+"/StringTie/{sample}_novel_transcripts.txt",
		novel_transcripts_list2 = config["dir"]["outdir"]+"/StringTie/{sample}_novel_transcripts_with_quotes.txt",
		novel_transcripts = config["dir"]["outdir"]+"/StringTie/{sample}_novel_transcripts.gtf",
	params:
		gtf_dir = config["dir"]["outdir"]+"/StringTie/{sample}",
		label = "{sample}",
		gffcompare_file = "{sample}.transcripts.gtf.tmap",
		column = "-f5" # -f4 for genes
	wildcard_constraints:
	    sample = '\w+'
	conda:
	    "conda/gffcompare.yaml"
	shell:
		"""
			cd {params.gtf_dir} && \
			gffcompare -r {input.ref_gff} {input.query_gtf} -o {params.label} && \
			cat {params.gffcompare_file} | awk '$3=="u"{{print $0}}' | cut {params.column} | sort | uniq > {output.novel_transcripts_list} && \
			sed 's/^/"/; s/$/"/' {output.novel_transcripts_list} > {output.novel_transcripts_list2} && \
			grep -F -f {output.novel_transcripts_list2} {input.query_gtf} > {output.novel_transcripts}	
		"""

## novel_transcript_fasta :                     GFFread is used to fetch FASTA for novel transcript from reference genome
rule novel_transcripts_fasta:
	input:
		novel_transcripts_gtf = config["dir"]["outdir"]+"/StringTie/{sample}_novel_transcripts.gtf",
		reference_fasta = config['ref_genome'],
		#reference_fasta = config["dir"]["ref_genome"]
	output:
		novel_transcripts_fasta =  config["dir"]["outdir"]+"/FASTA/{sample}_novel_transcripts.fa"
	wildcard_constraints:
	    sample = '\w+'
	conda:
	    "conda/gffread.yaml"
	shell:
		"""
			gffread -w {output.novel_transcripts_fasta} -g {input.reference_fasta} {input.novel_transcripts_gtf}
		"""

rule cpat:
	input:
		fasta =  config["dir"]["outdir"]+"/FASTA/{sample}_novel_transcripts.fa",
		logitmodel = config["dir"]["cpat"]["training-data"]['logitmodel'],
		hexamer = config["dir"]["cpat"]["training-data"]['hexamer']
	output:
		tsv = config["dir"]["outdir"]+"/cpat/output/{sample}.ORF_prob.tsv",
		tsv_lowercase = config["dir"]["outdir"]+"/cpat/output/{sample}.lowercase.tsv",
	params:
		outdir = config["dir"]["outdir"]+"/cpat/output"
	conda:
	    "conda/cpat.yaml"
	shell:
		"""
			mkdir -p {params.outdir} && cd {params.outdir}
			cpat.py -g {input.fasta}  -d {input.logitmodel} -x {input.hexamer} -o {wildcards.sample} && \
			
			# changing the case of Stringtie, CPAT changes Stringtie to STRINGTIE in previous rule 
			perl -p -e 's/STRINGTIE/Stringtie/g' {output.tsv} > {output.tsv_lowercase}
		"""
		
		
## filter_cpat :                                Separates non-coding transcripts from CPAT analysis using coding potential cutoff
rule filter_cpat:
	input:
		tsv = config["dir"]["outdir"]+"/cpat/output/{sample}.lowercase.tsv",
	output:
		nctsv = config["dir"]["outdir"]+"/cpat/output/{sample}.nc.tsv",
	params:
		cutoff = config['cpat_filter_cutoff']
	shell:
		"""
			awk -F "\t" '{{ if($10 <= {params.cutoff} ) {{print $1}} }}' {input.tsv} > {output.nctsv}
		"""

## sort_cpat :                                  Sorts non-coding transcripts for comparison with CPC prediction and remove ORF id
rule sort_cpat:
	input:
		nctsv = config["dir"]["outdir"]+"/cpat/output/{sample}.nc.tsv",
	output:
	    ncsortedtsv = config["dir"]["outdir"]+"/cpat/output/{sample}.nc.sorted.tsv",
	    nccorrectedtsv = config["dir"]["outdir"]+"/cpat/output/{sample}.nc.corrected.tsv"
	shell:
		"""
			sort {input.nctsv} > {output.ncsortedtsv}
			awk -F '_ORF_' '{{print $1}}' {output.ncsortedtsv} > {output.nccorrectedtsv}
		"""
## cpc :                                        Predict coding potential of novel transcripts using CPC
rule cpc:
	input:
		fasta =  config["dir"]["outdir"]+"/FASTA/{sample}_novel_transcripts.fa",
	output:
		tsv = config["dir"]["outdir"]+"/cpc/output/{sample}.tsv",
	params:
		outdir = config["dir"]["outdir"]+"/cpc/output"
	conda:
	    "conda/cpc.yaml"
	shell:
		"""
			mkdir -p {params.outdir}
			CPC2.py -i {input.fasta}  -o {output.tsv}
		"""


## filter_cpc :                                 Separates non-coding transcripts from CPC analysis
rule filter_cpc:
	input:
		tsv = config["dir"]["outdir"]+"/cpc/output/{sample}.tsv",
	output:
		nctsv = config["dir"]["outdir"]+"/cpc/output/{sample}.nc.tsv",
	params:
		filter = "noncoding"	
	shell:
		"""
			awk '{{ if($8 == "{params.filter}" ) {{print $1}} }}' {input.tsv} > {output.nctsv}
		"""


## sort_cpc :                                   Sorts non-coding transcripts for comparison with CPAT prediction
rule sort_cpc:
	input:
		nctsv = config["dir"]["outdir"]+"/cpc/output/{sample}.nc.tsv",
	output:
		ncsortedtsv = config["dir"]["outdir"]+"/cpc/output/{sample}.nc.sorted.tsv",
	shell:
		"""
			sort {input.nctsv} > {output.ncsortedtsv}
		"""

## cpat_cpc_intersect :                         Finds common non-coding transcripts predicted by both CPAT and CPC
rule cpat_cpc_intersect:
    input:
        cpat = config["dir"]["outdir"]+"/cpat/output/{sample}.nc.corrected.tsv",
        cpc = config["dir"]["outdir"]+"/cpc/output/{sample}.nc.sorted.tsv",
    output:
        comm = config["dir"]["outdir"]+"/nc/{sample}.common_nc.tsv",
    params:
    	outdir = config["dir"]["outdir"]+"/nc"
    shell:
        """
        	mkdir -p {params.outdir}
            comm -12 {input.cpat} {input.cpc} > {output.comm}
		"""


## non_coding_transcript :                      Fetches non-coding transcripts using output of 'cpat_cpc_intersect' from novel transcripts generated by StringTie
rule non_coding_transcript_from_cpat:
	input:
		tsv =  config["dir"]["outdir"]+"/nc/{sample}.common_nc.tsv",
		query_gtf = config["dir"]["outdir"]+"/StringTie/{sample}_novel_transcripts.gtf",
	output:
		tsv =  config["dir"]["outdir"]+"/nc/{sample}.quoted_common_nc.tsv",
		novel_transcripts = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.gtf",
	params:
		outdir = config["dir"]["outdir"]+"/GTF_nc_transcript"
	shell:
		"""
			mkdir -p {params.outdir}
			sed 's/^/"/; s/$/"/' {input.tsv} > {output.tsv} && \
			grep -F -f {output.tsv} {input.query_gtf} > {output.novel_transcripts}
		"""

## genelist_for_non_coding_transcripts :        Fetches gene list from GTF of 'non_coding_transcript'
rule genelist_for_non_coding_transcripts:
	input:
		nc_gtf =  config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.gtf",
	output:
		gene_list = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.tsv"

	shell:
		"""
			#keeping " to avoid matching with wrong gene ids
			perl -lne 'print @m if @m=(/((?:gene_id)\s+\S+)/g);' {input.nc_gtf} | perl -p -i -e 's/gene_id|;| //g' > {output.gene_list} && \
			sort {output.gene_list} | uniq
		"""


## all_genes_for_non_coding_transcript :        Fetches all transcript of gene list to look for coding isoforms
rule all_genes_for_non_coding_transcript:
	input:
		query_gtf = config["dir"]["outdir"]+"/StringTie/{sample}/transcripts.gtf",
		gene_list = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.tsv"
	output:
		transcripts = config["dir"]["outdir"]+"/GTF_iso_nc_transcript/{sample}.gtf"
	shell:
		"""
			# gene list is already quoted so not adding quotes again, could use Mikado
			grep -F -f {input.gene_list} {input.query_gtf} > {output.transcripts}
		"""

## gene_transcript_pair :                       Prepares gene id, transcript id list from GTF of 'non_coding_transcript' and 'all_genes_for_non_coding_transcript'
##                                              First GTF contains non-coding transcripts and later contains all isoforms of non-codong transcripts
rule gene_transcript_pair:
	input:
		nc_gtf = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.gtf",
		combine_gtf = config["dir"]["outdir"]+"/GTF_iso_nc_transcript/{sample}.gtf",
	output:
		nc_pair = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.pair.tsv",
		combine_pair = config["dir"]["outdir"]+"/GTF_iso_nc_transcript/{sample}.pair.tsv",
	params: 
			scripts = config["scripts"]
	shell:
		"""
			python {params.scripts}/tra.py < {input.nc_gtf} > {output.nc_pair} && \
			python {params.scripts}/tra.py < {input.combine_gtf} > {output.combine_pair}
		"""

## coding_isoform_of_non_coding_transcripts :   Prepares a list of genes which has coding isoforms of non-codong transcripts
rule coding_isoform_of_non_coding_transcripts:
	input:
		nc_pair = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.pair.tsv",
		combine_pair = config["dir"]["outdir"]+"/GTF_iso_nc_transcript/{sample}.pair.tsv",
	output:
		diff = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.ncgenes.tsv",
	params: 
			scripts = config["scripts"]
	shell:
		"""
			python {params.scripts}/filter.py {input.nc_pair} {input.combine_pair} > {output.diff}
		"""

## non_coding_genes :                           Prepares GTF of the genes with only non-codong transcripts
##                                              Adds quotes around id using sed 
##                                              uses reverse grep to fetch genes which are not present in input list 
rule non_coding_genes:
    input:
        gtf = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.gtf",
        genes = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.ncgenes.tsv"
    output:
        quoted_genes = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}.ncgenes.quoted.tsv",
        gtf = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}_filtered.gtf"
    shell:
        """
            sed 's/^/"/; s/$/"/' {input.genes} > {output.quoted_genes} && \
            grep -v -F -f {output.quoted_genes} {input.gtf} > {output.gtf}
        """

## filter_non_coding_genes :                    Filters non-coding genes based on transcript length 
rule filter_non_coding_genes:
	input: 
		gtf = config["dir"]["outdir"]+"/GTF_nc_transcript/{sample}_filtered.gtf",
	output: 
		gtf = config["dir"]["outdir"]+"/GTF_nc_transcript_filtered/{sample}.gtf",
	conda:
	    "conda/gffread.yaml"
	shell:
		"""
			gffread {input.gtf} -T -l 200 -o {output.gtf}
		"""