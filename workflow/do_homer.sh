cd /mnt/d/data/epiLungCancer/intermediate/vector/w4.2/lungWithGSE186458.leukocytes/homer || exit
nohup findMotifsGenome.pl ../LUAD_0.95_0.2.top500.hypo.smvr.bed hg38 LUAD_0.95_0.2.top500.hypo.smvr -mask > nohup.LUAD_0.95_0.2.top500.hypo.smvr.log  2>&1 &
nohup findMotifsGenome.pl ../LUAD_0.95_0.2.top500.hyper.smvr.bed hg38 LUAD_0.95_0.2.top500.hyper.smvr -mask > nohup.LUAD_0.95_0.2.top500.hyper.smvr.log  2>&1 &
nohup findMotifsGenome.pl ../LUSC_0.95_0.2.top500.hypo.smvr.bed hg38 LUSC_0.95_0.2.top500.hypo.smvr -mask > nohup.LUSC_0.95_0.2.top500.hypo.smvr.log  2>&1 &
nohup findMotifsGenome.pl ../LUSC_0.95_0.2.top500.hyper.smvr.bed hg38 LUSC_0.95_0.2.top500.hyper.smvr -mask > nohup.LUSC_0.95_0.2.top500.hyper.smvr.log  2>&1 &
nohup findMotifsGenome.pl ../LCC_0.95_0.2.top500.hypo.smvr.bed hg38 LCC_0.95_0.2.top500.hypo.smvr -mask > nohup.LCC_0.95_0.2.top500.hypo.smvr.log  2>&1 &
nohup findMotifsGenome.pl ../LCC_0.95_0.2.top500.hyper.smvr.bed hg38 LCC_0.95_0.2.top500.hyper.smvr -mask > nohup.LCC_0.95_0.2.top500.hyper.smvr.log  2>&1 &
nohup findMotifsGenome.pl ../SCLC_0.95_0.2.top500.hypo.smvr.bed hg38 SCLC_0.95_0.2.top500.hypo.smvr -mask > nohup.SCLC_0.95_0.2.top500.hypo.smvr.log  2>&1 &
nohup findMotifsGenome.pl ../SCLC_0.95_0.2.top500.hyper.smvr.bed hg38 SCLC_0.95_0.2.top500.hyper.smvr -mask > nohup.SCLC_0.95_0.2.top500.hyper.smvr.log  2>&1 &

cd /mnt/d/data/epiLungCancer/intermediate/dmc/p80/homer || exit
nohup findMotifsGenome.pl ../one2rest80.LUAD.top500.hypo.dmr.bed hg38 one2rest80.LUAD.top500.hypo.dmr -mask > nohup.one2rest80.LUAD.top500.hypo.dmr.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.LUAD.top500.hyper.dmr.bed hg38 one2rest80.LUAD.top500.hyper.dmr -mask > nohup.one2rest80.LUAD.top500.hyper.dmr.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.LUSC.top500.hypo.dmr.bed hg38 one2rest80.LUSC.top500.hypo.dmr -mask > nohup.one2rest80.LUSC.top500.hypo.dmr.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.LUSC.top500.hyper.dmr.bed hg38 one2rest80.LUSC.top500.hyper.dmr -mask > nohup.one2rest80.LUSC.top500.hyper.dmr.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.LCC.top500.hypo.dmr.bed hg38 one2rest80.LCC.top500.hypo.dmr -mask > nohup.one2rest80.LCC.top500.hypo.dmr.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.LCC.top500.hyper.dmr.bed hg38 one2rest80.LCC.top500.hyper.dmr -mask > nohup.one2rest80.LCC.top500.hyper.dmr.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.SCLC.top500.hypo.dmr.bed hg38 one2rest80.SCLC.top500.hypo.dmr -mask > nohup.one2rest80.SCLC.top500.hypo.dmr.log  2>&1 &
nohup findMotifsGenome.pl ../one2rest80.SCLC.top500.hyper.dmr.bed hg38 one2rest80.SCLC.top500.hyper.dmr -mask > nohup.one2rest80.SCLC.top500.hyper.dmr.log  2>&1 &
