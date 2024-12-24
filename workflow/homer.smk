from pathlib import Path

DATA_DIR = config['local']['dataDir']

data_dir = Path(DATA_DIR) / 'intermediate' / 'dmc' / 'p80'

class_groups = ['hypo.LUAD', 'hyper.LUAD', 'hypo.LUSC', 'hyper.LUSC', 'hypo.LCC', 'hyper.LCC', 'hypo.SCLC',
                'hyper.SCLC']
# TODO: fix-bug: Only one can be run at a time
class_groups = ['hypo.LUSC']

rule all:
    input:
        expand(data_dir / "homer" / "{sample}",sample=class_groups),


rule homer:
    input:
        data_dir / "one2rest80.{sample}.top500.dmr.bed"
    output:
        directory(data_dir / "homer" / "{sample}")
    shell:
        """
        findMotifsGenome.pl {input} hg38 {output} -mask
        """
