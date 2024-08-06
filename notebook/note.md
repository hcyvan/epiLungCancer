# Generate genomic signal file

```shell
methytools signal -i /mnt/d/data/epiLungCancer/raw/merge.group.3x.all.bed -e tss -o ./tmp/signal.tss.matrix.bed
methytools signal -i /mnt/d/data/epiLungCancer/raw/merge.group.3x.all.bed -e cgi -o ./tmp/signal.cgi.matrix.bed
methytools signal -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e cgi -o ./tmp/signal.cgi.matrix.bed
```

# Extract the mean methylation in a genomic region

```shell
methytools region-methy -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e cgi -o ./tmp/signal.cgi_n50_p50.matrix.bed -u -50 -d 50
methytools region-methy -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e tss -o ./tmp/signal.tss_n50_p50.matrix.bed -u -50 -d 50
methytools region-methy -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e cgi -o ./tmp/signal.cgi_p3900_p4000.matrix.bed -u 3900 -d 4000
methytools region-methy -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e tss -o ./tmp/signal.tss_p3900_p4000.matrix.bed -u 3900 -d 4000
methytools region-methy -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e cgi -o ./tmp/signal.cgi_n3900_n4000.matrix.bed -u -4000 -d -3900
methytools region-methy -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e tss -o ./tmp/signal.tss_n3900_n4000.matrix.bed -u -4000 -d -3900
```
# Extract methylation from the bed format methylation matrix

```shell
methytools extract -i /mnt/d/data/epiLungCancer/intermediate/one2rest.dmc.bed -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -o /mnt/d/data/epiLungCancer/intermediate/one2rest.dmc.beta.bed
```
```shell
methytools extract -i ../CTL.bed -m /mnt/d/data/epiLungCancer/raw/moabs/merge.group.3x.all.bed.gz -o CTL.smv.beta.bed &
methytools extract -i ../LUAD.bed -m /mnt/d/data/epiLungCancer/raw/moabs/merge.group.3x.all.bed.gz -o LUAD.smv.beta.bed &
methytools extract -i ../LUSC.bed -m /mnt/d/data/epiLungCancer/raw/moabs/merge.group.3x.all.bed.gz -o LUSC.smv.beta.bed &
methytools extract -i ../LCC.bed -m /mnt/d/data/epiLungCancer/raw/moabs/merge.group.3x.all.bed.gz -o LCC.smv.beta.bed &
methytools extract -i ../SCLC.bed -m /mnt/d/data/epiLungCancer/raw/moabs/merge.group.3x.all.bed.gz -o SCLC.smv.beta.bed &

```

# Further filtering of mcomp dmc results (mcomppost percentile)
```shell
mcomppost percentile -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.CTL.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.percentile.txt
mcomppost percentile -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUAD.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.LUAD.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUAD.percentile.txt
mcomppost percentile -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUSC.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.LUSC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUSC.percentile.txt
```
# Reverse the direction of comparison (mcomppost reverse)
```shell
mcomppost reverse -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.CTL.vs.LUAD_LUSC_LCC_SCLC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.txt
```
# Merge DMCs into DMRs (mcomppost dmc2dmr)
```shell
mcomppost dmc2dmr -i /mnt/d/data/epiLungCancer/intermediate/dmc/p80/one2rest80.dmc.bed -o /mnt/d/data/epiLungCancer/intermediate/dmc/p80/one2rest80.dmr.py.bed
mcomppost dmc2dmr -i /mnt/d/data/epiLungCancer/intermediate/dmc/p85/one2rest85.dmc.bed -o /mnt/d/data/epiLungCancer/intermediate/dmc/p85/one2rest85.dmr.py.bed
mcomppost dmc2dmr -i /mnt/d/data/epiLungCancer/intermediate/dmc/p90/one2rest90.dmc.bed -o /mnt/d/data/epiLungCancer/intermediate/dmc/p90/one2rest90.dmr.py.bed
```
# Homer
## subtype-specific dmr
```shell
nohup findMotifsGenome.pl ../one2rest80.CTL.hypo.dmr.bed hg38 ./CTL.hypo -mask > nohup.CTL.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.LUAD.hypo.dmr.bed hg38 ./LUAD.hypo -mask > nohup.LUAD.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.LUSC.hypo.dmr.bed hg38 ./LUSC.hypo -mask > nohup.LUSC.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.LCC.hypo.dmr.bed hg38 ./LCC.hypo -mask > nohup.LCC.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.SCLC.hypo.dmr.bed hg38 ./SCLC.hypo -mask > nohup.SCLC.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.CTL.hyper.dmr.bed hg38 ./CTL.hyper -mask > nohup.CTL.hyper.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.LUAD.hyper.dmr.bed hg38 ./LUAD.hyper -mask > nohup.LUAD.hyper.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.LUSC.hyper.dmr.bed hg38 ./LUSC.hyper -mask > nohup.LUSC.hyper.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.LCC.hyper.dmr.bed hg38 ./LCC.hyper -mask > nohup.LCC.hyper.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.SCLC.hyper.dmr.bed hg38 ./SCLC.hyper -mask > nohup.SCLC.hyper.log  2>&1 &
```
```shell
nohup findMotifsGenome.pl ../one2rest85.CTL.hypo.dmr.bed hg38 ./CTL.hypo -mask > nohup.CTL.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest85.LUAD.hypo.dmr.bed hg38 ./LUAD.hypo -mask > nohup.LUAD.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest85.LUSC.hypo.dmr.bed hg38 ./LUSC.hypo -mask > nohup.LUSC.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest85.LCC.hypo.dmr.bed hg38 ./LCC.hypo -mask > nohup.LCC.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest85.SCLC.hypo.dmr.bed hg38 ./SCLC.hypo -mask > nohup.SCLC.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest85.CTL.hyper.dmr.bed hg38 ./CTL.hyper -mask > nohup.CTL.hyper.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest85.LUAD.hyper.dmr.bed hg38 ./LUAD.hyper -mask > nohup.LUAD.hyper.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest85.LUSC.hyper.dmr.bed hg38 ./LUSC.hyper -mask > nohup.LUSC.hyper.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest85.LCC.hyper.dmr.bed hg38 ./LCC.hyper -mask > nohup.LCC.hyper.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest85.SCLC.hyper.dmr.bed hg38 ./SCLC.hyper -mask > nohup.SCLC.hyper.log  2>&1 &
```
```shell
nohup findMotifsGenome.pl ../one2rest90.CTL.hypo.dmr.bed hg38 ./CTL.hypo -mask > nohup.CTL.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest90.LUAD.hypo.dmr.bed hg38 ./LUAD.hypo -mask > nohup.LUAD.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest90.LUSC.hypo.dmr.bed hg38 ./LUSC.hypo -mask > nohup.LUSC.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest90.LCC.hypo.dmr.bed hg38 ./LCC.hypo -mask > nohup.LCC.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest90.SCLC.hypo.dmr.bed hg38 ./SCLC.hypo -mask > nohup.SCLC.hypo.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest90.CTL.hyper.dmr.bed hg38 ./CTL.hyper -mask > nohup.CTL.hyper.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest90.LUAD.hyper.dmr.bed hg38 ./LUAD.hyper -mask > nohup.LUAD.hyper.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest90.LUSC.hyper.dmr.bed hg38 ./LUSC.hyper -mask > nohup.LUSC.hyper.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest90.LCC.hyper.dmr.bed hg38 ./LCC.hyper -mask > nohup.LCC.hyper.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest90.SCLC.hyper.dmr.bed hg38 ./SCLC.hyper -mask > nohup.SCLC.hyper.log  2>&1 &
```
## specific methylation vectors window
```shell
nohup findMotifsGenome.pl ../CTL.bed hg38 ./CTL -mask > nohup.CTL.log  2>&1 &
nohup findMotifsGenome.pl ../LUAD.bed hg38 ./LUAD -mask > nohup.LUAD.log  2>&1 &
nohup findMotifsGenome.pl ../LUSC.bed hg38 ./LUSC -mask > nohup.LUSC.log  2>&1 &
nohup findMotifsGenome.pl ../LCC.bed hg38 ./LCC -mask > nohup.LCC.log  2>&1 &
nohup findMotifsGenome.pl ../SCLC.bed hg38 ./SCLC -mask > nohup.SCLC.log  2>&1 &

```

```shell
nohup findMotifsGenome.pl ../CTL.onco.bed hg38 ./CTL -mask > nohup.CTL.log  2>&1 &
nohup findMotifsGenome.pl ../LUAD.onco.bed hg38 ./LUAD -mask > nohup.LUAD.log  2>&1 &
nohup findMotifsGenome.pl ../LUSC.onco.bed hg38 ./LUSC -mask > nohup.LUSC.log  2>&1 &
nohup findMotifsGenome.pl ../LCC.onco.bed hg38 ./LCC -mask > nohup.LCC.log  2>&1 &
nohup findMotifsGenome.pl ../SCLC.onco.bed hg38 ./SCLC -mask > nohup.SCLC.log  2>&1 &

```
# pattools
```shell
pattools vector-diff -i /mnt/d/data/epiLungCancer/intermediate/vector/20240612/merge.motif.gz -g 4 > SCLC.txt
pattools region-file -t cpg2genome --column col2 --offset-col2-start-and-end 3 --out-format bed -c /mnt/d/project/wgbs_tools/references/hg38/CpG.bed.gz -i SCLC.txt -o SCLC.bed
```