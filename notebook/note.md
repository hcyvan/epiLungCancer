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
methytools extract -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.percentile.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz
methytools extract -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.percentile.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.bed
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
