library('biomaRt')
library('dplyr')
library('stringr')
library('splitstackshape')
datasets <- listDatasets(ensembl)
grep('drerio', datasets, value=TRUE)
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl")
mart3 = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
mart4 = useMart("ensembl", dataset="drerio_gene_ensembl")
mart5 = useMart("ensembl", dataset="celegans_gene_ensembl")
mart6 = useMart("ensembl", dataset="xtropicalis_gene_ensembl")

human_ensembl_ids=as.list(getBM(attributes = c( 'ensembl_gene_id'),mart = mart1))

# human / human
human_human=getLDS(attributes=c("ensembl_gene_id"),
                       filters="ensembl_gene_id", mart=mart1,
                       attributesL=c("ensembl_gene_id"), martL=mart1, values = human_ensembl_ids$ensembl_gene_id)
# human / mmusculus
human_mmusculus=getLDS(attributes=c("ensembl_gene_id"),
       filters="ensembl_gene_id", mart=mart1,
       attributesL=c("ensembl_gene_id"), martL=mart2, values = human_ensembl_ids$ensembl_gene_id)
# human / dmelanogaster
human_dmelanogaster=getLDS(attributes=c("ensembl_gene_id"),
                   filters="ensembl_gene_id", mart=mart1,
                   attributesL=c("ensembl_gene_id"), martL=mart3, values = human_ensembl_ids$ensembl_gene_id)
# human / drerio
human_drerio=getLDS(attributes=c("ensembl_gene_id"),
                        filters="ensembl_gene_id", mart=mart1,
                        attributesL=c("ensembl_gene_id"), martL=mart4, values = human_ensembl_ids$ensembl_gene_id)
# human / celegans
human_celegans=getLDS(attributes=c("ensembl_gene_id"),
                        filters="ensembl_gene_id", mart=mart1,
                        attributesL=c("ensembl_gene_id"), martL=mart5, values = human_ensembl_ids$ensembl_gene_id)
# human / xtropicalis
human_xtropicalis=getLDS(attributes=c("ensembl_gene_id"),
                        filters="ensembl_gene_id", mart=mart1,
                        attributesL=c("ensembl_gene_id"), martL=mart6, values = human_ensembl_ids$ensembl_gene_id)

# Combine all human orthologs
all_human_orthologs=rbind(human_human,human_mmusculus,human_dmelanogaster,human_drerio,human_celegans,human_xtropicalis)

# Read human TOS genes and allspecies_TOS_motifs files
human_TOS_genes=read.csv('all_human_genes_with_TOS_motif.txt', sep=' ', header=FALSE, stringsAsFactors = FALSE)
allspecies_TOS_motifs=read.csv('allspecies_TOS_motifs_final.txt', sep='\t',header=FALSE, stringsAsFactors = FALSE)

# Select all human orthologs, which have at least 1 TOS motif in human protein sequence
all_human_orthologs_withhumanTOS=all_human_orthologs %>% filter(str_detect(Gene.stable.ID,paste0(trimws(as.list(human_TOS_genes$V1)),collapse="|")))

# Merge allspecies_TOS_motifs and all_human_orthologs_withhumanTOS
final_TOS_table=merge(allspecies_TOS_motifs, all_human_orthologs_withhumanTOS, by.x = 5, by.y='Gene.stable.ID.1', all = TRUE)
final_TOS_table$V2=NULL
final_TOS_table$V15=NULL
colnames(final_TOS_table)=c('ensembl_gene_id','genome','peptide_id','chromosome','ensembl_transcript_id','gene_biotype','transcript_biotype','gene_symbol','TOS_motif','motif_start','motif_end','peptide_length','motif_location','orthologous_human_ensembl_gene_id')

# Merge human TOS genes with the TOS table that contains all orthologs
final_TOS_table=merge(final_TOS_table,human_TOS_genes,by.x='orthologous_human_ensembl_gene_id',by.y=1,all=TRUE)
colnames(final_TOS_table)[15]='orthologous_human_gene_symbol'
final_TOS_table$genome_=stringr::str_extract(final_TOS_table$ensembl_gene_id, "^.{4}")

# Annotate TOS motif as 5', middle and 3' motif
final_TOS_table$motif_location_category=ifelse(final_TOS_table$motif_location<=0.25,"5' motif",ifelse(final_TOS_table$motif_location>=0.75,"3' motif", "middle"))

# Calculate TOS count based on number of orthologs that have TOS motif
TOS_counts=as.data.frame(final_TOS_table %>% select(orthologous_human_ensembl_gene_id,genome) %>% filter(genome!=NaN) %>% distinct() %>% count(orthologous_human_ensembl_gene_id))

# Calculate ortholog count based on the number of orthologs
ortholog_counts=as.data.frame(final_TOS_table %>% select(orthologous_human_ensembl_gene_id,genome_) %>% distinct() %>% count(orthologous_human_ensembl_gene_id))
colnames(TOS_counts)[2]="TOS_count"
colnames(ortholog_counts)[2]="ortholog_count"
final_TOS_table_=merge(final_TOS_table, TOS_counts, by='orthologous_human_ensembl_gene_id', all = TRUE)
final_TOS_table__=merge(final_TOS_table_, ortholog_counts, by='orthologous_human_ensembl_gene_id', all = TRUE)
final_TOS_table__=cSplit(indt = final_TOS_table__,splitCols ='peptide_id',sep=".",direction = 'wide')
final_TOS_table__$gene_id=paste0(final_TOS_table__$peptide_id_1,".",final_TOS_table__$peptide_id_2)
final_TOS_table__$peptide_id_1=NULL
final_TOS_table__$peptide_id_2=NULL
final_TOS_table__$peptide_id_3=NULL

# Write the results
write.table(final_TOS_table__,'final_TOS_conservation_table.txt', sep='\t', quote=F,row.names = F)
