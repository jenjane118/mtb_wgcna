## function for verifying srnas with presence of tss
# Jennifer J. Stiens
# j.j.stiens@gmail.com
# v2 16/12/2021

tss_srnas <- function(gene_list, refseq_name, tss_file) {
  # search for tss in sequence of srna and make df
  
  library(GenomicRanges)
  library(stringr)
  
  # gene_list <- all_genes
  # refseq_name = "AL123456.3"
  # annot_file  = here("ref_seqs/Mtb_h37rv.ASM19595v2_AL123456.3.gff3")
  # tss_file    = here("Data/shell_cortes_srna_tss.txt")
  
  # create grange object of srna coordinates
  # parse list for sRNAs predicted from baerhunter
  og_SRNAs <- gene_list[grepl("sRNA", gene_list)]
  if (length(og_SRNAs)==0){
    return(NULL)
    }else{
    srnas <- str_split_fixed(og_SRNAs, "sRNA:", 2)[,2]
    SRNA_df <- data.frame(matrix(0, nrow=length(srnas), ncol=4),stringsAsFactors = F)
    colnames(SRNA_df) <- c("pred_srna", "start", "stop", "strand")
    SRNA_df$pred_srna <- og_SRNAs
  }
  # parse coordinates and strand from snra name
  for (i in 1:length(srnas)){
    strand <- substr(srnas[i], 1, 1)
    start  <- str_split_fixed(srnas[i], "_", 2)[,1]
    start  <- sub(".", "", start)
    stop   <- str_split_fixed(srnas[i], "_", 2)[,2]
    SRNA_df$strand[i] <- strand
    SRNA_df$start[i]  <- as.numeric(start)
    SRNA_df$stop[i]   <- as.numeric(stop)
  }
  SRNA_df$strand <- sub('p', '+', SRNA_df$strand)
  SRNA_df$strand <- sub('m', '-', SRNA_df$strand)

  # create grange object of srna coordinates for each strand
  # subset srna_df by strand 
  fwd_SRNA_df <- SRNA_df[SRNA_df$strand=="+",]
  rev_SRNA_df <- SRNA_df[SRNA_df$strand=="-",]
  
  srna_fwd_Granges<-GRanges(seqnames = refseq_name,
                            ranges   = IRanges(fwd_SRNA_df$start, 
                                               end = fwd_SRNA_df$stop),
                            strand   = Rle(strand(fwd_SRNA_df$strand)))
  fwd_df<-as.data.frame(srna_fwd_Granges)
  srna_rev_Granges<-GRanges(seqnames = refseq_name,
                            ranges   = IRanges(rev_SRNA_df$start, 
                                               end = rev_SRNA_df$stop),
                            strand   = Rle(strand(rev_SRNA_df$strand)))
  rev_df<-as.data.frame(srna_rev_Granges)
  # read in tss and make into granges obj
  tss<-read.delim(tss_file, sep=" ")
  tss_gr<-GRanges(seqnames = refseq_name,
                  ranges   = IRanges(start = tss$Genome.position, 
                                     end = tss$Genome.position),
                  strand   = Rle(strand(tss$Strand))
  )
  # find TSS in srnas
  ## find overlap of each to tss within 10 bp of start 
  # forward strand tss
  fwd_tss_overlaps <- findOverlaps(
    srna_fwd_Granges,
    tss_gr,
    type = "start",
    maxgap = 20,
    select = "first",
    ignore.strand = FALSE
  )
  qh<-as.matrix(fwd_tss_overlaps)
  fwd_df$tss<-ifelse(is.na(qh), FALSE, tss[qh])
  
  # reverse strand tss
  rev_tss_overlaps<- findOverlaps(
    srna_rev_Granges,
    tss_gr,
    type = "end",
    maxgap = 20,
    select = "first",
    ignore.strand = FALSE
  )
  qh<-as.matrix(rev_tss_overlaps)
  rev_df$tss<-ifelse(is.na(qh), FALSE, tss[qh])
  #put back together
  df_tss<-rbind(fwd_df, rev_df)
  # add TSS into utr df
  SRNA_df <- left_join(SRNA_df, df_tss, 
                      by=c("start"="start", "stop"="end", "strand"="strand"), 
                      keep=FALSE) %>%
             select(-c("width", "seqnames"))
  return(SRNA_df)
}

name_srnas <- function(srna_df, refseq_name, annot_file){
  ## function to assign name to each putative sRNA based on Lamichhane, 2012
  ## if not overlapping any ORF on either strand: Rv1XXXX(c)
  ## XXXX = locus of preceding ORF (ORF coming before sRNA)
  ## if overlapping ORF on either strand: RvXXXX(c)
  ## XXXX = locus of overlapping ORF on either strand
  ## 'c' = sRNA on minus strand
  
  # use separate gRanges for + and - strand sRNAs 
  # use 'overlap' to find gene overlapping sRNA, ignore.strand=T
  # can use 'follow' to find gene preceding sRNA, ignore.strand=T
  # assign name and put 'c' on - strand named sRNAs
  
  library(GenomicRanges)
  library(rtracklayer)
  library(stringr)
  
  #test data
  #ref_seq_name = "AL123456.3"
  #annot_file=here("ref_seqs/Mtb_h37rv.ASM19595v2_AL123456.3.gff3")
  #srna_df <- SRNA_df
  
  # read in gff for ref genome (use same one as for BH prediction) as genomic ranges object
  ref_granges <- import(annot_file)
  # filter for type='gene'
  ref_granges <- ref_granges[mcols(ref_granges)$type=="gene"]
  ref_genes   <- as.data.frame(ref_granges)
  
  # parse list from module genes, or module/trait df
  # uses dataframe generated from tss_rna function
  
  # create grange object of srna coordinates for each strand
  ## we do need to know start and stop because if not overlapping, 
  ## we need to find gene that precedes start of srna 
  ## for negative strand, use 'precede' rather than follow
  pos_df<-srna_df[srna_df$strand=="+",]
  
  if (nrow(pos_df) != 0){
    pos_Granges <- GRanges(seqnames = refseq_name,
                           ranges   = IRanges(pos_df$start, 
                                              end = pos_df$stop),
                           strand   = Rle(strand(pos_df$strand))
    )
    ## does srna overlap any gene? 
    #ignore strand=T, because want overlap on either strand
    # if srna is on positive strand:
    for (i in 1:nrow(pos_df)){
      fwd_overlaps <- findOverlaps(
        pos_Granges[i],
        ref_granges,
        type = "any",
        minoverlap = 1,
        maxgap = -1,
        select = "last",  #want most downstream gene overlap for naming
        ignore.strand = TRUE
      )
      overlap_gene<-ref_genes$gene_id[fwd_overlaps]
      pos_df$ov_orf[i]<-overlap_gene
      if (is.na(overlap_gene)==FALSE){
        overlap_name <- substr(overlap_gene, 3, 6)
        pos_df$srna_name[i] <- paste("ncRv", overlap_name, sep="")
      }else{
        ## use follow for the genes not already named
        # follow is locus that 'is followed' by sRNA (comes before sRNA==upstream)
        ug<-GenomicRanges::follow(pos_Granges[i], ref_granges, ignore.strand=T)
        follow_locus <- ref_genes$gene_id[ug]
        follow_name <- substr(follow_locus, 3, 6)
        pos_df$srna_name[i] <- paste("ncRv1", follow_name, sep="")
      }
    }
  }
  
  # if srna is on negative strand
  rev_df <- srna_df[srna_df$strand == "-",]
  if (nrow(rev_df) != 0){
    rev_Granges <- GRanges(seqnames = refseq_name,
                           ranges = IRanges(rev_df$start,
                                            end = rev_df$stop),
                           strand = Rle(strand(rev_df$strand))
    )
    
    for (i in 1:nrow(rev_df)){
      rev_overlaps <- findOverlaps(
        rev_Granges[i],
        ref_granges,
        type = "any",
        minoverlap = 1,
        maxgap = -1,
        select = "first",
        ignore.strand = TRUE
      )
      overlap_gene<-ref_genes$gene_id[rev_overlaps]
      rev_df$ov_orf[i]<-overlap_gene
      if (is.na(overlap_gene)==FALSE){
        overlap_name <- substr(overlap_gene, 3, 6)
        rev_df$srna_name[i] <- paste("ncRv", overlap_name, "c", sep="")
      }else{
        ## use follow for the genes not already named
        # follow is gene that 'is followed' by sRNA (comes before sRNA==upstream)
        # this is reversed for negative strand: use precede
        ug<-GenomicRanges::precede(rev_Granges[i], ref_granges, ignore.strand=T)
        follow_locus <- ref_genes$gene_id[ug]
        follow_name <- substr(follow_locus, 3, 6)
        rev_df$srna_name[i] <- paste("ncRv1", follow_name, "c", sep="")
      }
    }
  }
  # join back together
  new_srna_df <- rbind(pos_df, rev_df)
  new_srna_df <- new_srna_df[order(new_srna_df$start),]
  return(new_srna_df)
}
