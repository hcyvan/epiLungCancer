from pathlib import Path

WORK_PLACE = config['local']['methylationLevel']['home']
GROUPS = config['local']['groups']
methyMatrix = config['local']['methyMatrix']

data_dir = Path(WORK_PLACE)


rule all:
    input:
        expand(data_dir / "one2rest80.{group}.dmr.bed",group=GROUPS),
        expand(data_dir / "one2rest80.{group}.top500.hyper.dmr.bed",group=GROUPS),
        expand(data_dir / "one2rest80.{group}.top500.hypo.dmr.bed",group=GROUPS),
        expand(data_dir / "one2rest80.{group}.top500.hyper.dmr.sample.beta.bed", group=GROUPS),
        expand(data_dir / "one2rest80.{group}.top500.hypo.dmr.sample.beta.bed", group=GROUPS)

rule dmc2dmr:
    input:
        data_dir / "one2rest80.dmc.bed"
    output:
        data_dir / "one2rest80.dmr.bed",
    shell:
        """
        mcomppost dmc2dmr -i {input} --value-column=5 -o {output}        
        """

rule top500:
    input:
        data_dir / "one2rest80.dmr.bed"
    output:
        dmr=data_dir / "one2rest80.{sample}.dmr.bed",
        hypo=data_dir / "one2rest80.{sample}.top500.hypo.dmr.bed",
        hyper=data_dir / "one2rest80.{sample}.top500.hyper.dmr.bed"
    shell:
        """
        set +euo pipefail
        cat {input}|grep {wildcards.sample} > {output.dmr}
        cat {output.dmr}|grep hypo|grep -v chrX|grep -v chrY|sort -k5 -k6 -k7 -rn | head -n 500 > {output.hypo}
        cat {output.dmr}|grep hyper|grep -v chrX|grep -v chrY|sort -k5 -k6 -k7 -rn | head -n 500 > {output.hyper}
        """

rule methyExtract:
    input:
        hypo = data_dir / "one2rest80.{sample}.top500.hypo.dmr.bed",
        hyper = data_dir / "one2rest80.{sample}.top500.hyper.dmr.bed"
    output:
        hypo = data_dir / "one2rest80.{sample}.top500.hypo.dmr.sample.beta.bed",
        hyper = data_dir / "one2rest80.{sample}.top500.hyper.dmr.sample.beta.bed"

    shell:
        """
        methytools extract -i {input.hypo} -m {methyMatrix} -o {output.hypo}
        methytools extract -i {input.hyper} -m {methyMatrix} -o {output.hyper}
        """
