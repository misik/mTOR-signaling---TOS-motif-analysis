# Download protein sequence files for all 6 genomes
wget ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-98/fasta/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-98/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.22.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-98/fasta/caenorhabditis_elegans/pep/Caenorhabditis_elegans.WBcel235.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-98/fasta/xenopus_tropicalis/pep/Xenopus_tropicalis.Xenopus_tropicalis_v9.1.pep.all.fa.gz

# Format the protein sequence file to make it searchable by the motif finder script
for file in *all.fa; do zcat  $file | tr '\n' '#' | sed 's/#>/\n>/g' |\
 sed 's/#/@/' | tr -d '#' > $(basename $file)"_edited".fa; done

# Run motif finder to find all proteins with TOS motifs
cat *_edited.fa |python find_all_motifs.py > allspecies_TOS_motifs.txt
cat allspecies_TOS_motifs.txt |tr ',' '\t'|tr -d "[]' " | \
sed 's/#\t/\n/g' | awk NF | awk -F'\t' '{gsub("\\..*$","",$5)}1' OFS='\t'> allspecies_TOS_motifs_final.txt

# Select the human genes with TOS motifs
cat allspecies_TOS_motifs_final.txt|grep GRCh38|awk -F'\t' '{print $5,$9}'|awk '!seen[$0]++' > all_human_genes_with_TOS_motif.txt
