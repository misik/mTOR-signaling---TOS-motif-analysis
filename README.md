## mTOR-signaling---TOS-motif-analysis

This repository contains the code for TOS (mTOR signaling) motif analysis in 6 different proteomes.

**Predicted TOS motif:** 
```
"F[DEAVP][FMILV][DEVRL][MLIEFYAR]" 
```

**Selected proteomes:**
```
Homo_sapiens.GRCh38.pep.all.fa.gz
Mus_musculus.GRCm38.pep.all.fa.gz
Danio_rerio.GRCz11.pep.all.fa.gz
Drosophila_melanogaster.BDGP6.22.pep.all.fa.gz
Caenorhabditis_elegans.WBcel235.pep.all.fa.gz
Xenopus_tropicalis_v9.1.pep.all.fa.gz
```

**Step 1:**

Run find_all_motifs.sh to find all TOS motifs in the selected proteomes

```
bash find_all_motifs.sh 
```

**Step 2:**

Find all biomart orthologs of human genes and merge them with the TOS motif containing proteins. Calculate number of TOS motif containing orthologs and total orthologs. Annotate TOS motif as 5', middle and 3' motif based on its locations within the sequence.

```
Rscript biomart-ortholog-finder.R
```


