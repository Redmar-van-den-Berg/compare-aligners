
from itertools import combinations 
include: "common.smk"

rule all:
    input:
        images = expand("images/{version}.sif", version=config["aligners"]),
        samfiles = expand("sam/{version}.sam", version=config["aligners"]),
        bamfiles = expand("bam/{version}.sorted.bam", version=config["aligners"]),
        cadbure = [f"{aligner1}--{aligner2}" for aligner1, aligner2 in combinations(config["aligners"],2) ]

rule get_container:
    output:
        "images/{version}.sif"
    params:
        lambda wc: config["aligners"][wc.version]
    log:
        "log/get_container.{version}.txt"
    shell: """
        singularity build {output} {params}
    """
    
rule build_gmap_database:
    input:
        config["reference"]
    output:
        directory("gmap_database")
    log:
        "log/build_gmap_database.err"
    container:
        containers["gsnap-2014.12.23"]
    shell: """
    gmap_build -D `pwd` \
            -d gmap_database \
            {input} 2> {log}
    """

rule align:
    input:
        f1 = config["forward"],
        f2 = config["reverse"],
        index = rules.build_gmap_database.output,
        image = "images/{version}.sif"
    output: 
        "sam/{version}.sam"
    log:
        "log/align.{version}.txt"
    shell: """
        singularity run {input.image} gsnap \
            --dir `dirname {input.index}` \
            --db `basename {input.index}` \
            --batch 4 \
            --novelsplicing 1 \
            --npaths 1 \
            --quiet-if-excessive \
            --format sam --gunzip {input.f1} {input.f2} > {output} \
            2> {log}
    """

rule sort:
    input:
        samfile = "sam/{version}.sam"
    output:
        bam = "bam/{version}.sorted.bam",
        #bai = "bam/{version}.sorted.bai"
    log:
        "log/sort.{version}.txt"
    container:
        containers["picard"]
    shell: """
        mkdir -p bam;
        picard -Xmx4G SortSam \
            I={input.samfile} \
            O={output.bam} \
            SORT_ORDER=queryname \
            CREATE_INDEX=true
    """

rule cadbure:
    input:
        aligner1 = "bam/{version1}.sorted.bam",
        aligner2 = "bam/{version2}.sorted.bam",
    output:
        folder = directory("{version1}--{version2}"),
        html = "{version1}--{version2}/cadbure.html"
    log:
        stderr = "log/cadbure-{version1}--{version2}.err",
        stdout = "log/cadbure-{version1}--{version2}.out"
    container:
        containers["cadbure"]
    shell: """
        set -e
        mkdir -p {wildcards.version1}--{wildcards.version2};

        cd {wildcards.version1}--{wildcards.version2}

        CADBURE -f {wildcards.version1} \
                -fb ../{input.aligner1} \
                -s {wildcards.version2} \
                -sb ../{input.aligner2} -o cadbure \
                2> ../{log.stderr} > ../{log.stdout}
    """
