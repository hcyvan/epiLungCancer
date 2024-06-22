# Generate genomic signal file

```shell
methytools signal -i /mnt/d/data/epiLungCancer/raw/merge.group.3x.all.bed -e tss -o ./tmp/signal.tss.matrix.bed
methytools signal -i /mnt/d/data/epiLungCancer/raw/merge.group.3x.all.bed -e cgi -o ./tmp/signal.cgi.matrix.bed
methytools signal -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e cgi -o ./tmp/signal.cgi.matrix.bed
```

# Extract the mean methylation in a region

```shell
methytools region-methy -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e cgi -o ./tmp/signal.cgi_n50_p50.matrix.bed -u -50 -d 50
methytools region-methy -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e tss -o ./tmp/signal.tss_n50_p50.matrix.bed -u -50 -d 50
methytools region-methy -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e cgi -o ./tmp/signal.cgi_p3900_p4000.matrix.bed -u 3900 -d 4000
methytools region-methy -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e tss -o ./tmp/signal.tss_p3900_p4000.matrix.bed -u 3900 -d 4000
methytools region-methy -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e cgi -o ./tmp/signal.cgi_n3900_n4000.matrix.bed -u -4000 -d -3900
methytools region-methy -i /mnt/d/data/epiLungCancer/raw/merge.d3.all.bed.gz -e tss -o ./tmp/signal.tss_n3900_n4000.matrix.bed -u -4000 -d -3900
```

# Further filtering of mcomp dmc results (mcomppost dmc-filter)
## 80 percentile
```shell
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.CTL.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.filtered_80.txt -p 0.8
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUAD.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.LUAD.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUAD.filtered_80.txt -p 0.8
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUSC.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.LUSC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUSC.filtered_80.txt -p 0.8
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LCC.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.LCC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LCC.filtered_80.txt -p 0.8
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.SCLC.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.SCLC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.SCLC.filtered_80.txt -p 0.8
```
## 85 percentile
```shell
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.CTL.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.filtered_85.txt -p 0.85
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUAD.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.LUAD.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUAD.filtered_85.txt -p 0.85
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUSC.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.LUSC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUSC.filtered_85.txt -p 0.85
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LCC.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.LCC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LCC.filtered_85.txt -p 0.85
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.SCLC.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.SCLC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.SCLC.filtered_85.txt -p 0.85
```
## 90 percentile
```shell
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.CTL.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.filtered_90.txt -p 0.90
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUAD.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.LUAD.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUAD.filtered_90.txt -p 0.90
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUSC.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.LUSC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LUSC.filtered_90.txt -p 0.90
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LCC.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.LCC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.LCC.filtered_90.txt -p 0.90
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.SCLC.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.gz -t ./tmp/sample.SCLC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.SCLC.filtered_90.txt -p 0.90
```
# Reverse the direction of comparison (mcomppost reverse)
```shell
mcomppost reverse -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.CTL.vs.LUAD_LUSC_LCC_SCLC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.txt
```