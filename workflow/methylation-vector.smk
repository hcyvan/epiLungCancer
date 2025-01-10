from pathlib import Path

MVC_DATA = config['local']['mv']['mvc']

data_dir = Path(MVC_DATA).parent

groups = ['LUAD', 'LUSC', 'LCC', 'SCLC']


def myfunc(wildcards):
    group_dict = {
        "LUAD": "1",
        "LUSC": "2",
        "LCC": "3",
        "SCLC": "4",
    }
    return [group_dict[x] for x in wildcards]


rule all:
    input:
        #expand(data_dir / "{sample}_0.95_0.2.smvc",sample=groups),
        #expand(data_dir / "{sample}_0.95_0.2.smvr.bed",sample=groups),
        expand(data_dir / "{sample}_0.95_0.2.top500.hypo.smvr.bed",sample=groups),
        expand(data_dir / "{sample}_0.95_0.2.top500.hyper.smvr.bed",sample=groups),

rule separating:
    params:
        group=myfunc
    output:
        data_dir / "{sample}_0.95_0.2.smvc"
    shell:
        """
        pattools mv-separating -i {MVC_DATA} -g {params.group} --frac-mvs 0.95 --frac-samples 0.2 --with-meta -o {output}
        """

rule merge:
    input:
        data_dir / "{sample}_0.95_0.2.smvc"
    output:
        data_dir / "{sample}_0.95_0.2.smvr.bed"
    shell:
        """
        set +euo pipefail
        pattools smvc-merge -i {input} -o {output} -e chrX,chrY
        """

rule top500:
    input:
        data_dir / "{sample}_0.95_0.2.smvr.bed"
    output:
        hypo=data_dir / "{sample}_0.95_0.2.top500.hypo.smvr.bed",
        hyper=data_dir / "{sample}_0.95_0.2.top500.hyper.smvr.bed"
    shell:
        """
        set +euo pipefail
        cat {input}|grep hypo|sort -k4 -k5 -n -r |head -n500 > {output.hypo}
        cat {input}|grep hyper|sort -k4 -k5 -n -r |head -n500 > {output.hyper}
        """
