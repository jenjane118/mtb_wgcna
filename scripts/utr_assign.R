## function for assigning categories to utrs
## Written by: Jennifer J. Stiens, 1 May, 2021

## update v2: 15 Dec, 2021


utr_assign <- function(gene_list, refseq_name, annot_file, tss_file) {
  # returns a dataframe with each predicted utr given an assigment
  # 5', 3' and related gene or 'NA'
  # gene_list: list of genes including predicted utrs from baerhunter
  # refseq_name: name of reference sequence that will be on granges, ex: "AL123456.3"
  # ref_seq: reference annotation file for making granges object
  
  library(GenomicRanges)
  library(rtracklayer)
  library(stringr)
  
  #gene_list = geneInfo_df$gene_ID
  #refseq_name = "AL123456.3"
  #annot_file= here("ref_seqs/Mtb_h37rv.ASM19595v2_AL123456.3.gff3")
  #tss_file= here("/WGCNA_12_2021/Data/shell_cortes_srna_tss.txt")
  
  # read in gff for ref genome (use same one as for BH prediction) as genomic ranges object
  ref_granges <- import(annot_file)
  # filter for type='gene'
  ref_granges <- ref_granges[mcols(ref_granges)$type=="gene"]
  ref_genes   <- as.data.frame(ref_granges)
  # col'ID' has 'gene:Rv0001' and 'gene_id' has just ORF name 'Rv0001'
  # parse list for UTRs predicted from baerhunter
  
  UTRs <- gene_list[grepl("UTR", gene_list)]
  if (length(UTRs)==0){
    return(NULL)
  }else{
    utrs <- str_split_fixed(UTRs, "UTR:", 2)[,2]
    UTR_df           <- data.frame(matrix(0, nrow=length(utrs), ncol=4),stringsAsFactors = F)
    colnames(UTR_df) <- c("pred_utr", "start", "stop", "strand")
    UTR_df$pred_utr  <- UTRs
  }
  # parse coordinates and strand from utr name
  for (i in 1:length(utrs)){
    strand <- substr(utrs[i], 1, 1)
    start  <- str_split_fixed(utrs[i], "_", 2)[,1]
    start  <- sub(".", "", start)
    stop   <- str_split_fixed(utrs[i], "_", 2)[,2]
    UTR_df$strand[i] <- strand
    UTR_df$start[i]  <- as.numeric(start)
    UTR_df$stop[i]   <- as.numeric(stop)
  }
  UTR_df$strand <- sub('p', '+', UTR_df$strand)
  UTR_df$strand <- sub('m', '-', UTR_df$strand)
  
  # create grange object of utr coordinates
  utr_Granges<-GRanges(seqnames = refseq_name,
                         ranges   = IRanges(UTR_df$start, 
                                            end = UTR_df$stop),
                         strand   = Rle(strand(UTR_df$strand))
                         )

  # searching whether utr co-locates with gene?
  # precede is gene that "is preceded" by the UTR (comes after UTR==downstream)
  pg<-GenomicRanges::precede(utr_Granges, ref_granges,ignore.strand=F)
  utr_precede_genes<-ref_genes$gene_id[pg]
  # precede doesn't pay attention to circular genome, use first gene in gff
  utr_precede_genes[is.na(utr_precede_genes)] <- ref_genes$gene_id[1]
  # follow is gene that 'is followed' by UTR (comes before UTR==upstream)
  ug<-GenomicRanges::follow(utr_Granges, ref_granges, ignore.strand=F)
  utr_follow_genes<-ref_genes$gene_id[ug]
  utr_follow_genes[is.na(utr_follow_genes)] <- ref_genes$gene_id[nrow(ref_genes)]
  ## nearest(utrs, refseq, select="all", ignore.strand = F) will give me nearest gene feature.
  # select=all means returns both nearest features in case of a tie, arbitrary chooses one.
  # see how many ties
  a<-nearest(utr_Granges, ref_granges, select="all", ignore.strand=F)
  b<-as.matrix(a)
  #show frequency of hits for each utr (1 hit or 2 hits, 2=tie)
  c<-as.data.frame(table(b[,1]))
  utr_tied<-c$Freq
  # see nearest gene (arbitrarily choose in ties)
  ng<-nearest(utr_Granges, ref_granges, select="arbitrary", ignore.strand=F)
  utr_nearest_genes<-ref_genes$gene_id[ng]

  #include preceding (downstream gene), following (upstream gene), and nearest genes in df 
  # (and indicate if there is a tie)
  UTR_df$downstream <-utr_precede_genes 
  UTR_df$upstream   <-utr_follow_genes  
  UTR_df$nearest    <-utr_nearest_genes
  UTR_df$tied       <-utr_tied
  
  # if nearest == precede/downstream (gene preceded by/downstream of UTR), 5' UTR; 
  # if nearest == follow/upstream (gene followed by UTR / UTR upstream of gene), 3'UTR
  
  # if tied==look for tss
  # find TSS in utrs
  # subset utr_granges by strand 
  utr_fwd_Granges<-utr_Granges[strand(utr_Granges) == "+"]
  fwd_df<-as.data.frame(utr_fwd_Granges)
  utr_rev_Granges<-utr_Granges[strand(utr_Granges) =="-"]
  rev_df<-as.data.frame(utr_rev_Granges)
  # read in tss and make into granges obj
  tss<-read.delim(tss_file, sep=" ")
  tss_gr<-GRanges(seqnames = refseq_name,
                  ranges   = IRanges(tss$Genome.position, 
                                     end = tss$Genome.position),
                  strand   = Rle(strand(tss$Strand))
  )
  
  ## find overlap of each to tss within 10 bp of start  
  # (this needs to be "end" for negative strand)
  # forward strand tss
  fwd_tss_overlaps <-
    findOverlaps(
      utr_fwd_Granges,
      tss_gr,
      type = "start",
      maxgap = 10,
      select = "first",
      ignore.strand = FALSE
    )
  qh<-as.matrix(fwd_tss_overlaps)
  fwd_df$tss<-ifelse(is.na(qh), FALSE, tss[qh])
 
  # reverse strand tss
  rev_tss_overlaps<-
    findOverlaps(
      utr_rev_Granges,
      tss_gr,
      type = "end",
      maxgap = 10,
      select = "first",
      ignore.strand = FALSE
    )
  qh<-as.matrix(rev_tss_overlaps)
  rev_df$tss<-ifelse(is.na(qh), FALSE, tss[qh])
  #put back together
  df_tss<-rbind(fwd_df, rev_df)
  # add TSS into utr df
  UTR_df <- left_join(UTR_df, df_tss, 
                           by=c("start"="start", "stop"="end", "strand"="strand"), 
                           keep=FALSE) %>%
                select(-c("width", "seqnames"))
  
  # limit identified 5' UTRs to those with TSS (5' UTRs only)
  # 3' UTRs as those that are upstream of a gene (and finish 20-40 nt before downstream gene)
  
  # only need start of downstream gene to see if 3' UTR (end of UTR > 20nts from downstream start)
  for (i in 1:nrow(UTR_df)){
    if (UTR_df$strand[i] == "+"){
      UTR_df$downstream_start[i]  <-
        ref_genes$start[match(UTR_df$downstream[i], ref_genes$gene_id)]
      # consider edge case of final utr on strand
      UTR_df$dist_to_start[i] <- ifelse(UTR_df$downstream_start[i] - UTR_df$stop[i] > 0,
        UTR_df$downstream_start[i] - UTR_df$stop[i], NA)
    }else{
      # for minus strand, need end coordinate of downstream gene
      UTR_df$downstream_start[i]  <-
        ref_genes$end[match(UTR_df$downstream[i], ref_genes$gene_id)]
      UTR_df$dist_to_start[i] <- ifelse(UTR_df$start[i] - UTR_df$downstream_start[i] > 0,
        UTR_df$start[i] - UTR_df$downstream_start[i], NA) }
  }
  # assign 5', 3' or btwn to each utr 
  
  for (i in 1:nrow(UTR_df)){
    # if tied (same dist between two genes, located between two genes)
    if (UTR_df$tied[i] == 2){
      UTR_df$utr[i] <- paste("UTR_BTWN", UTR_df$downstream[i], sep="_")
    #not tied (closer to either upstream or downstream gene)
    }else{
      if (UTR_df$nearest[i] == UTR_df$upstream[i]) {
        # check distance from end of UTR to start of downstream gene
        if (is.na(UTR_df$dist_to_start[i]) | UTR_df$dist_to_start[i] > 40){
          UTR_df$utr[i] <- paste("3UTR", UTR_df$upstream[i], sep="_")
        }else{
          UTR_df$utr[i] <- paste("UTR_BTWN", UTR_df$downstream[i], sep="_") }
      }else{
          # if nearest==downstream, but no TSS: UTR precedes gene, likely 5' UTR
          UTR_df$utr[i] <- paste("5UTR", UTR_df$downstream[i], sep="_") }
    }
  }  
  return(UTR_df)
}


