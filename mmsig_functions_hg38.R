#' Print message
#'
#' @param dbg boolean whether to print or not
#' @param mess message to print
#' @param ... additional parameters
#'
#' @return message

spit = function (dbg, mess, ...) {
  if (dbg) { print( sprintf(mess,...) ) }
}

#' Clean input VCF-like file
#'
#' @param snvs input mutation catalog in 5-columns format
#'
#' @return clean mutation catalog in 5-columns format
#'
vcfCleaner <- function(snvs){
  
  snvs <- snvs[snvs$ref %in% c("A","C","G","T") & snvs$alt %in% c("A","C","G","T"),] # remove indels
  snvs$chr <- as.character(snvs$chr)
  snvs <- snvs[snvs$chr %in% paste0("chr", c(1:22, "X")),]
  snvs$pos <- as.numeric(snvs$pos)
  snvs <- snvs[stats::complete.cases(snvs),]
  return(snvs)
}

#' Statistical tests for transcriptional strand bias
#'
#' @param data_5cols input mutational catalog in 5-column format
#' @return list of results from strand bias in each trinucleotide context and pooled analysis for SBS-MM1 and SBS35
#' @importFrom dplyr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom tidyr separate
#' @importFrom dplyr filter
#' @importFrom MutationalPatterns strand_bias_test
#' @importFrom reshape2 dcast
#' @importFrom dplyr rowwise
#' @importFrom stats poisson.test
#' @importFrom TxDb.Hsapiens.UCSC.hg38.knownGene TxDb.Hsapiens.UCSC.hg38.knownGene
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @import BSgenome.Hsapiens.UCSC.hg38
#'

getStrandBias <- function(data_5cols){
  
  # Generate GRange list from data in vcf-like format
  sample_list <- list()
  samples <- unique(data_5cols$sample)
  
  for(i in (1:length(samples))) {
    temp <- data_5cols[data_5cols$sample== samples[i],]
    temp_gr <- with(temp,
                    GRanges(chr,
                            IRanges(start=pos, end=pos),
                            REF= ref,
                            ALT=alt))
    genome(temp_gr) <- "hg38"
    
    sample_list[[i]]<- (temp_gr)
  }
  
  names(sample_list) <- samples
  
  # transcriptional strand annotation
  genes_hg19 <- suppressMessages(GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
  
  mut_mat_stranded <- MutationalPatterns::mut_matrix_stranded(vcf_list = sample_list,
                                                              ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
                                                              ranges = genes_hg19,
                                                              mode = "transcription")
  
  mut_df_stranded <- as.data.frame(mut_mat_stranded) %>%
    rownames_to_column(var = "type") %>%
    melt(id.var = "type", variable.name = 'group', value.name = 'no_mutations') %>%
    separate(col = 'type', into = c('type', 'strand'), sep = "-")
  
  # Test for all transcriptional strand bias by 3-nt and sample
  mut_strandtest <- suppressMessages(strand_bias_test(mut_df_stranded))
  
  # Test for mm1
  mm1 <- c("C[C>T]A","G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T")
  
  mut_mm1 <- mut_df_stranded %>%
    filter(type %in% mm1) %>%
    dplyr::group_by(group, strand) %>%
    dplyr::summarise(count = sum(no_mutations)) %>%
    as.data.frame() %>%
    dcast(group ~ strand, value.var = "count") %>%
    rowwise() %>%
    dplyr::mutate(ratio = transcribed/untranscribed) %>%
    as.data.frame()
  
  # Applying the same poisson test as the strand_bias_test function
  
  p_poisson <- list()
  for(i in 1:nrow(mut_mm1)){
    transcribed <- mut_mm1[i,"transcribed"]
    untranscribed <- mut_mm1[i,"untranscribed"]
    pval <- poisson.test(c(transcribed, untranscribed), r = 1)$p.value
    
    p_poisson[[i]] <- ifelse(is.numeric(pval), pval, NA)
  }
  
  mut_mm1$p_poisson <- unlist(p_poisson)
  
  mut_mm1 <- dplyr::mutate(mut_mm1, MM1_flag = ifelse(p_poisson < 0.05 & ratio > 1, "*", ""))
  
  # Test for platinum / SBS35
  
  SBS35 <- c("C[C>A]C", "C[C>T]C")
  
  mut_SBS35 <- mut_df_stranded %>%
    filter(type %in% SBS35) %>%
    dplyr::group_by(group, strand) %>%
    dplyr::summarise(count = sum(no_mutations)) %>%
    as.data.frame() %>%
    dcast(group ~ strand, value.var = "count") %>%
    rowwise() %>%
    dplyr::mutate(ratio = transcribed/untranscribed) %>%
    as.data.frame()
  
  # Applying the same poisson test as the strand_bias_test function
  
  p_poisson <- list()
  for(i in 1:nrow(mut_SBS35)){
    transcribed <- mut_SBS35[i,"transcribed"]
    untranscribed <- mut_SBS35[i,"untranscribed"]
    pval <- poisson.test(c(transcribed, untranscribed), r = 1)$p.value
    
    p_poisson[[i]] <- ifelse(is.numeric(pval), pval, NA)
  }
  
  mut_SBS35$p_poisson <- unlist(p_poisson)
  
  mut_SBS35 <- dplyr::mutate(mut_SBS35, SBS35_flag = ifelse(p_poisson < 0.05 & ratio > 1, "*", ""))
  
  # output
  
  output <- list(all_3nt = mut_strandtest,
                 mm1 = mut_mm1,
                 SBS35 = mut_SBS35)
  
  return(output)
}


#' Print message and variable
#'
#' @param dbg boolean whether to print or not
#' @param name message to print
#' @param var variable to print
#'
#' @return name and variable
#'
spat = function (dbg, name, var) {
  if (dbg) {print(name); print(var)}
}


#' Mutational signature reference
#'
#' Mutational signature reference containing the signatures previously identified in multiple myeloma
#' and other hematological cancers, including relevant chemotherapeutic agents
#'
#' @format A data frame with 96 rows and 13 variables:
#' \describe{
#'   \item{class}{nucleotide change and trinucleotide context}
#'   \item{sub}{nucleotide change}
#'   \item{tri}{trinucleotide context}
#'   \item{SBS}{proportion of each mutational class in each of the SBS references included}
#'   ...
#' }
#' @source \url{https://cancer.sanger.ac.uk/cosmic/signatures/SBS/index.tt}
"signature_ref"



#' Check if the input mutation catalog is aligned to the correct reference genome
#'
#' @param snvs input mutation catalog in 5-columns format
#' @param genome BSgenome object with the full reference genome used for analysis
#'
#' @return reference genome match: TRUE or FALSE
#'
refCheck <- function(snvs, genome){
  # number of mutations to check
  nCheck <- ifelse(nrow(snvs) > 100, 100, nrow(snvs))
  # subset mutations to check
  sub <- snvs[sample(nrow(snvs), nCheck),]
  ## Get the reference genome positions
  newRef <- BSgenome::getSeq(genome, as.character(sub$chr), sub$pos, sub$pos)
  # Check reference genome positions
  out <- all(newRef == sub$ref)
  return(out)
}

#' Example data in 96 classes format
#'
#' Complete mutational catalogs from whole genome sequencing of 82 multiple myeloma samples
#'
#' @format Mutational count matrix with 96 rows and 82 variables:
#' \describe{
#'   \item{row names}{mutational class}
#'   \item{column names}{sample ID}
#'   ...
#' }
#' @source \url{https://www.nature.com/articles/s41467-020-15740-9}
"mm_96_classes"


#' Example data in 5-columns format
#'
#' Complete mutational catalogs from whole genome sequencing of 82 multiple myeloma samples
#'
#' @format A data frame with 445813 rows and 5 variables:
#' \describe{
#'   \item{sample}{sample ID}
#'   \item{chr}{chromosome}
#'   \item{pos}{position on the hg19 reference genome}
#'   \item{ref}{nucleotide in the hg19 reference genome}
#'   \item{alt}{alternative nucleotide in the sample}
#'   ...
#' }
#' @source \url{https://www.nature.com/articles/s41467-020-15740-9}
"mm_5_col"



#' Mutational signature fitting with mmsig
#'
#' @param muts.input mutation input data frame as specified under input.format
#' @param sig.input mutational signature reference with mutational classes as rows and signature exposures as columns: Substitution.Type, Trinucleotide, signature.1, ... signature.n
#' @param input.format vcf: five column vcf-like data frame with the following columns: sample, chr (e.g. chr1), pos, ref, alt. classes: samples as columns and the 96 mutational classes as rows
#' @param sample.sigt.profs NULL = use all signatures provided in the reference. Optionally provide list with signatures to consider for each sample
#' @param bootstrap TRUE/FALSE for whether bootstrapping is to be performed
#' @param iterations number of bootstrapping iterations to perform (only if bootstrap == TRUE)
#' @param strandbias TRUE/FALSE for whether transcriptional strand bias should be tested for (only for vcf-like input format)
#' @param refcheck check that input mutational catalog (if vcf-format) is aligned to hg19
#' @param cos_sim_threshold cosine similarity threshold below which signatures are removed from the final profile
#' @param force_include vector with the names of signatures to always keep in the final profile of every sample
#' @param dbg FALSE = silent; TRUE = verbose
#' @importFrom dplyr left_join
#' @importFrom deconstructSigs mut.to.sigs.input
#' @import BSgenome.Hsapiens.UCSC.hg38
#'
#' @return mutational signature fitting results for all samples
#' @export
#'
mm_fit_signatures = function(muts.input,
                             sig.input,
                             input.format = "vcf",
                             sample.sigt.profs=NULL,
                             bootstrap=FALSE,
                             iterations=1000,
                             strandbias=FALSE,
                             refcheck=TRUE,
                             cos_sim_threshold=0.01,
                             force_include=c("SBS1", "SBS5"),
                             dbg=FALSE) {
  "

    "
  
  options(scipen = 999)
  
  #####################################################
  ########        Read/assign input data       ########
  #####################################################
  
  # Mutational signature reference
  consigts.defn <- sig.input
  
  # Input mutation data
  if(input.format == "vcf") {
    
    if(sum(c("sample", "chr", "pos", "ref", "alt") %in% names(muts.input)) != 5){
      stop("ERROR: Please provide the following data columns: sample, chr, pos, ref, alt")
    }
    
    muts.input <- vcfCleaner(muts.input)
    
    if(nrow(muts.input) == 0){
      stop("ERROR: Inappropriate input data format")
    }
    
    if(refcheck){
      if(!refCheck(muts.input, BSgenome.Hsapiens.UCSC.hg38)){
        stop("ERROR: Wrong reference genome, please provide mutational data aligned to hg19/GRCh37")
      }
    }
    
    # Input sample data
    samples.muts <- mut.to.sigs.input(mut.ref = muts.input,
                                      sample.id = "sample",
                                      chr = "chr",
                                      pos = "pos",
                                      ref = "ref",
                                      alt = "alt",
                                      bsg = BSgenome.Hsapiens.UCSC.hg38)
    
    samples.muts <- as.data.frame(t(samples.muts))
    
  } else if(input.format == "classes"){
    samples.muts <- muts.input
    samples.muts <- samples.muts[names(samples.muts) != "Total"]                       # remove totals column
    samples.muts <- samples.muts[1:96,,drop=FALSE]                                     # remove additional rows (e.g. "Total)"
  } else{
    stop("ERROR: Invalid input format, please provide a vcf-like format ('vcf') or 96 mutational classes ('classes')")
  }
  
  # Process signature reference
  tcons.defn <- t(consigts.defn[,3:ncol(consigts.defn)])
  mutlist = paste(consigts.defn[,2],paste(substr(consigts.defn[,2],1,1),substr(consigts.defn[,1],3,3),substr(consigts.defn[,2],3,3),sep=""),sep=">")
  consigts.defn <- sapply(consigts.defn[,3:ncol(consigts.defn)], as.numeric)  #n be cautious in the original the first two columns are included
  rownames(consigts.defn)<-mutlist
  ref_signatures <- colnames(consigts.defn) # names of signatures in reference
  
  # Process sample mutational profiles
  samples <- colnames(samples.muts)
  
  # Assign which signatures to fit to each sample
  if(is.list(sample.sigt.profs)){
    spit(dbg, "using mm signature profiles from input argument")
    sigt.profs <- sample.sigt.profs
  } else {
    spit(dbg, "defaulting to use all signatures in the provided reference")
    mm.sigts <- ref_signatures
    
    sigt.profs <- list()
    for (i in 1:length(samples)) {
      sigt.profs[[samples[i]]] <- mm.sigts
    }
  }
  
  # Define output variable
  output <-list()
  
  #####################################################
  ########       Fit mutational signatures     ########
  #####################################################
  
  sigfit <- fit_signatures(samples.muts=samples.muts,
                           consigts.defn=consigts.defn,
                           sigt.profs=sigt.profs,
                           cos_sim_threshold=cos_sim_threshold,
                           force_include=force_include,
                           dbg=dbg)
  
  output$estimate <- sigfit
  
  #####################################################
  ########             Bootstrapping           ########
  #####################################################
  
  if(bootstrap){
    
    # Point estimates
    sig_est <- sigfit
    sig_est$sample <- row.names(sig_est)
    sig_est <- melt(sig_est[names(sig_est) != "mutations"],
                    id.vars = "sample",
                    variable.name = 'signature',
                    value.name = 'estimate')
    
    sig_est$signature <- as.character(sig_est$signature)
    
    sigboot <- bootstrap_fit_signatures(samples.muts = samples.muts,
                                        consigts.defn = consigts.defn,
                                        sigt.profs = sigt.profs,
                                        iterations = iterations,
                                        cos_sim_threshold = cos_sim_threshold,
                                        force_include = force_include)
    
    sigboot <- left_join(sigboot, sig_est, by = c("sample", "signature"))
    
    output$bootstrap <- sigboot
  }
  
  #####################################################
  ########              Strand Bias            ########
  #####################################################
  
  if(strandbias & input.format == "vcf"){
    
    strand_bias_out <- getStrandBias(muts.input)
    
    output$strand_bias_all_3nt <- strand_bias_out$all_3nt
    output$strand_bias_mm1 <- strand_bias_out$mm1
    output$strand_bias_SBS35 <- strand_bias_out$SBS35
    
  } else if(strandbias & input.format != "vcf") {
    warning("Transcriptional strand bias cannot be estimated from 96 classes input. \n Please provide vcf-like input format.")
  }
  
  output$mutmatrix <- samples.muts
  
  return(output)
}


#' Signature fitting function for mmsig
#'
#' @param samples.muts 96 classes mutational profile for all samples
#' @param consigts.defn 96 classes mutational profile for mutational signature reference
#' @param sigt.profs which signatures to fit for each sample
#' @param cos_sim_threshold cosine similarity threshold to remove signatures
#' @param force_include vector with the names of signatures to always keep in the final profile of every sample
#' @param dbg boolean whether to print or not
#'
#' @return estimated mutational signature profile of all samples
#'
fit_signatures = function(samples.muts,
                          consigts.defn,
                          sigt.profs,
                          cos_sim_threshold,
                          force_include,
                          dbg=dbg) {
  max.em.iter=2000
  consigt.names <- colnames(consigts.defn)
  samples <- colnames(samples.muts)
  
  # Mutational signature fitting procedure for each individual sample
  sigt.fraction = array(NA,dim=c(length(consigt.names), length(samples)))
  rownames(sigt.fraction) <- consigt.names
  colnames(sigt.fraction) <- samples
  
  for (j in 1:length(samples)) {
    sample.mut.freqs = as.numeric(samples.muts[,j])
    sample.mut.freqs[is.na(sample.mut.freqs)] = 0
    sample.sigts <- unique(sigt.profs[[ samples[j] ]])
    sample.sigts <- sample.sigts[match(consigt.names[consigt.names %in% sample.sigts], sample.sigts)]
    sample.consigts.defn  <- consigts.defn[, colnames(consigts.defn) %in% sample.sigts]
    spat(dbg, "colnames sample.consigts.defn (before em)", colnames(sample.consigts.defn))
    alpha <- em_signatures(sigts.defn=sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
    spat(dbg, "alpha", alpha)
    sample.consigts.defn <- sample.consigts.defn[, colnames(sample.consigts.defn) %in% names(alpha)]   # output sample.... should be identical to input sample....
    sampleAlpha <- alpha[match(colnames(sample.consigts.defn), names(alpha))]
    
    if (!all(alpha==sampleAlpha)) {stop("non-identical alphas")}
    
    spat(dbg, "colnames: sample.consigts.defn (after em (and reduction))", colnames(sample.consigts.defn))
    reconstructed <- sample.consigts.defn %*% alpha * sum(sample.mut.freqs)
    sample.cos.sim.meas <- MutationalPatterns::cos_sim_matrix(reconstructed, matrix(sample.mut.freqs, ncol=1))
    spat(dbg, "sample.cos.sim.meas", sample.cos.sim.meas)
    
    rem.alpha <- sampleAlpha                     # holds the final result
    rem.sample.consigts.defn <- sample.consigts.defn
    spit(dbg, "length of rem.alpha: %d", length(rem.alpha))
    reducing = TRUE
    
    # Signature profile shrinkage by cosine similarity (removing signatures that are not necessary to explain profile)
    while (reducing) {
      spat(dbg, "in the while, rem.alpha: ", rem.alpha)
      cosReduction <- NULL
      rem.names <- setdiff(names(rem.alpha), force_include)
      if(length(rem.names) == 0){ ## Avoiding script crash when only the forced signatures are present.
        spit(dbg, "removed all signatures except forced inclusion: exiting while...")
        break
      }
      
      for(c in rem.names){
        spit(dbg, "doing c: %s", c)
        red.sample.consigts.defn <- rem.sample.consigts.defn[,colnames(rem.sample.consigts.defn)!=c]
        red.alpha <- em_signatures(sigts.defn=red.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
        red.reconstructed <- red.sample.consigts.defn %*% red.alpha * sum(sample.mut.freqs)
        red.cos.sim.meas <- MutationalPatterns::cos_sim_matrix(red.reconstructed, matrix(sample.mut.freqs, ncol=1))
        cosReduction <- c(cosReduction, sample.cos.sim.meas-red.cos.sim.meas)
      }
      names(cosReduction) <- rem.names
      if (min(cosReduction) < cos_sim_threshold) {
        spit(dbg, "removing: %s", names(cosReduction)[which.min(cosReduction)])
        rem.sample.consigts.defn <- rem.sample.consigts.defn[,- which(colnames(rem.sample.consigts.defn)==names(which.min(cosReduction)))]
        rem.alpha <-  em_signatures(sigts.defn=rem.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
        reducing = TRUE
      }
      else {
        spit(dbg, "exiting while...")
        reducing = FALSE
      }
    }
    
    spit(dbg,"... while exited")
    rem.alpha.names <- names(rem.alpha)
    
    for (n in 1:length(consigt.names)) {
      if (consigt.names[n] %in% rem.alpha.names) {
        sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
      }
      else {
        sigt.fraction[n,j] <- 0
      }
    }
  }
  spat(dbg, "sigt.fraction", sigt.fraction)
  tdf.sigt.fraction <- as.data.frame(t(sigt.fraction))
  colsums.samples.muts <- colSums(samples.muts)
  sig <- cbind(tdf.sigt.fraction, "mutations"=colsums.samples.muts)
  colnames(sig)[colnames(sig)=="SBS.MM1"]<-"SBS-MM1"
  
  return(sig)
}



#' Expectation maximization algorithm to estimate the signature contribution
#'
#' @param sigts.defn mutational signature definitions
#' @param mut.freqs mutational profile of sample
#' @param max.iter maximum number of iterations
#' @param dbg boolean whether to print or not
#'
#' @return estimated mutational signature contributions
#'
em_signatures = function(sigts.defn, mut.freqs, max.iter, dbg) {
  nof.sigts.defn <- ncol(sigts.defn)
  alpha = stats::runif(nof.sigts.defn); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (i in 1:max.iter) {
    contr = t(array(alpha, dim=c(nof.sigts.defn,96))) * sigts.defn
    probs = contr/array(rowSums(contr), dim=dim(contr))
    probs[is.na(probs)] = 0
    probs = probs * mut.freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  spit(dbg, "em: exit iteration: %d", i)
  
  return( alpha )
}



#' Bootstrapping function for mmsig
#'
#' @param samples.muts 96 classes mutational profile for all samples
#' @param consigts.defn 96 classes mutational profile for mutational signature reference
#' @param sigt.profs which signatures to fit for each sample
#' @param cos_sim_threshold cosine similarity threshold to remove signatures
#' @param force_include vector with the names of signatures to always keep in the final profile of every sample
#' @param iterations number of mutational profiles to draw from the multinomial distribution of each sample
#' @importFrom plyr create_progress_bar
#' @importFrom plyr progress_text
#' @importFrom dplyr bind_rows
#'
#' @return For each mutational signature, returns its relative contribution as point estimate and bootstrapping mean with 95 % CI
#'
bootstrap_fit_signatures <- function(samples.muts,
                                     consigts.defn,
                                     sigt.profs,
                                     cos_sim_threshold,
                                     force_include,
                                     iterations = 1000){
  
  # Setup
  samples <- names(samples.muts)
  classes <- row.names(samples.muts)
  
  # List to populate with signatures
  mutSigs <- list()
  
  # Setup progress bar
  pbar <- create_progress_bar('text')
  pbar$init(length(samples))
  
  for(i in 1:length(samples)){
    # Loop through samples, generating a data frame of signature contributions for each
    sub <- as.integer(samples.muts[classes,i])
    total <- sum(sub)
    # sample new 96-classes profiles from the multinomial distribution
    bootMat <- data.frame(rmultinom(n = iterations, size = total, prob = sub/total))
    row.names(bootMat) <- classes
    
    # prepare the signatures to fit for each sample
    sig.prof <- list()
    
    for(s in 1:ncol(bootMat)){
      sig.prof[[names(bootMat)[s]]] <- sigt.profs[[samples[i]]]
    }
    
    ### Run mmSig
    sig_out <- fit_signatures(samples.muts=bootMat,
                              consigts.defn=consigts.defn,
                              sigt.profs=sig.prof,
                              cos_sim_threshold=cos_sim_threshold,
                              force_include=force_include,
                              dbg=FALSE)
    
    mutSigs[[i]] <- sig_out
    pbar$step()
  }
  
  names(mutSigs) <- samples
  
  # Generate final summary data frame
  
  ## Summary statistics
  my_summary <- function(x){
    c(mean(x), quantile(x, probs = 0.025), quantile(x, probs = 0.975))
  }
  
  mutSigsSummary <- list()
  for(i in 1:length(mutSigs)){
    s <- names(mutSigs)[i]
    temp <- mutSigs[[i]]
    out <- data.frame(t(sapply(temp[names(temp) != "mutations"], my_summary)))
    names(out) <- c('mean', 'CI025', 'CI975')
    out$signature <- row.names(out)
    out$sample <- s
    out <- out[c('sample', 'signature', 'mean', 'CI025', 'CI975')]
    mutSigsSummary[[i]] <- out
  }
  
  mutSigsSummary <- bind_rows(mutSigsSummary)
  
  return(mutSigsSummary)
}







