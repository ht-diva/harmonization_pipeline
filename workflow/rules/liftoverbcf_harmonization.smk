rule liftover_bcftools:
    input:
        sumstats=get_sumstats,
    output:
        vcf=temp(ws_path("temp/{sumstat_id}/{sumstat_id}.liftover.vcf.gz")),
    conda:
        "../envs/liftover_bcftools.yaml"
    params:
        hg37=config.get("hg37_fasta_file_path"),
        hg38=config.get("hg38_fasta_file_path"),
        chain_file=config.get("chain_file_path"),
    shell:
        "workflow/scripts/liftover_bcftools.sh {input.sumstats} {params.hg37} {params.hg38} {params.chain_file} {output.vcf} {threads}"


rule harmonize_sumstats:
    input:
        sumstats=rules.liftover_bcftools.output.vcf,
    output:
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz"),
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.log"),
    conda:
        "../envs/gwaspipe.yaml"
    params:
        format=config.get("params").get("harmonize_sumstats").get("input_format"),
        config_file=config.get("params").get("harmonize_sumstats").get("config_file"),
        output_path=config.get("workspace_path"),
        sumstats_sep=config.get("sumstats_sep"),
    shell:
        "gwaspipe "
        "-f {params.format} "
        "-c {params.config_file} "
        "-s '{params.sumstats_sep}' "
        "-i {input.sumstats} "
        "-o {params.output_path}"


rule bgzip_tabix:
    input:
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz"),
    output:
        ws_path("outputs/{sumstat_id}/{sumstat_id}.gwaslab.tsv.gz.tbi"),
    conda:
        "../envs/bgzip_tabix.yaml"
    shell:
        "workflow/scripts/bgzip_tabix.sh {input} {threads}"


rule create_snp_mapping_table:
    input:
        sumstats=expand(ws_path("temp/{sumstat_id}/{sumstat_id}.liftover.vcf.gz"), sumstat_id=analytes.sumstat_id.iloc[0]),
        #sumstats=get_sumstat(),
    output:
        ws_path("snp_mapping/table.snp_mapping.tsv.gz"),
    conda:
        "../envs/gwaspipe.yaml"
    params:
        format=config.get("params").get("snp_mapping").get("input_format"),
        config_file=config.get("params").get("snp_mapping").get("config_file"),
        output_path=config.get("workspace_path"),
        sumstats_sep=config.get("sumstats_sep"),
    shell:
        "gwaspipe "
        "-f {params.format} "
        "-c {params.config_file} "
        "-s '{params.sumstats_sep}' "
        "-i {input.sumstats} "
        "-o {params.output_path}"
