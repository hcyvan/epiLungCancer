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
```shell
mcomppost dmc-filter -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.CTL.vs.LUAD_LUSC_LCC_SCLC.txt -m /mnt/d/data/epiLungCancer/raw/moabs/merge.d3.all.bed.g
z -t CTL_4,CTL_10,CTL_16,CTL_41,CTL_47,CTL_58,CTL_59,CTL_64,CTL_72,CTL_73,CTL_74,CTL_78,CTL_100,CTL_111,CTL_125,CTL_126,CTL_133,CTL_141,CTL_149,CTL_154,CTL_157,CTL_161,CTL_174,CTL_178,CTL_181,CTL_183,CTL_186,CTL_225,CTL_228,CTL_289,CTL_290,CTL_291

```
# Reverse the direction of comparison (mcomppost reverse)
```shell
mcomppost reverse -i /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.CTL.vs.LUAD_LUSC_LCC_SCLC.txt -o /mnt/d/data/epiLungCancer/raw/moabs/dmc/dmc.Rest.vs.CTL.txt
```