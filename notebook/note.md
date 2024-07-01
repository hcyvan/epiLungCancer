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