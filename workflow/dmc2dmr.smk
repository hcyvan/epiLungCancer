from pathlib import Path

DATA_DIR = config['local']['dataDir']

data_dir = Path(DATA_DIR) / 'intermediate' / 'dmc' / 'p80'

class_groups = ['hypo.LUAD', 'hyper.LUAD', 'hypo.LUSC', 'hyper.LUSC', 'hypo.LCC', 'hyper.LCC', 'hypo.SCLC',
                'hyper.SCLC']

rule all:
    input:
        expand(data_dir / "one2rest80.{group}.top500.dmr.bed",group=class_groups)

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
        data_dir / "one2rest80.{sample}.top500.dmr.bed"
    shell:
        """
        set +euo pipefail
        cat {input}|grep {wildcards.sample} | sort -rk 7|head -n500 > {output}
        """
