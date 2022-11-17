#' Build the SCE object from the counts matrix
#'
#' @param counts Counts matrix with rownames corresponding to gene names and colnames corresponding to cell barcodes.
#'   Can be a counts matrix with intron counts specified by -I or an alevin-fry "USA" matrix,
#'   with intron counts marked by "-U" and ambiguous counts "-A".
#' @param intron_mode Logical indicating if the files included alignment to intronic regions.
#'   If TRUE, the main "counts" assay will contain counts aligned to introns and an additional "spliced"
#'   assay will contain counts aligned to spliced cDNA only.
#'   Default is FALSE.
#' @param round_counts Logical indicating in the count matrix should be rounded to integers on import.
#'   Default is TRUE.
#'
#' @return SingleCellExperiment object containing either just a counts assay with spliced cDNA only if
#'   `intron_mode` is FALSE. If `intron_mode` is TRUE, the counts assay will contain both spliced and unspliced
#'    counts and the spliced assay will contain the counts for just the spliced cDNA.
#'
#' @examples
#' \dontrun{
#'
#' # build a SCE object with only spliced cDNA as the main counts assay
#' build_sce(counts)
#'
#' # build a SCE object with unspliced cDNA as the main counts assay and spliced
#' # counts as a second assay
#' build_sce(counts,
#'           intron_mode = TRUE)
#'
#'}
build_sce <- function(counts,
                      intron_mode = FALSE,
                      round_counts= TRUE){

  if(!is.logical(intron_mode)){
    stop("intron_mode must be set as TRUE or FALSE")
  }

  if(!is.logical(round_counts)){
    stop("round_counts must be set as TRUE or FALSE")
  }

  # if using intron mode and want to include the unspliced data get counts for both unspliced and spliced
  if(intron_mode) {

    unspliced <- collapse_intron_counts(counts, which_counts = c("unspliced"))
    spliced <- collapse_intron_counts(counts, which_counts = c("spliced"))

    # before creating the SCE object we need to check that the dimensions match
    # cells will be the same
    num_cells <- ncol(unspliced)
    # genes may not match so grab the genes that are found as unspliced only or spliced only
    unspliced_only_genes <- rownames(unspliced)[!rownames(unspliced) %in% rownames(spliced)]
    spliced_only_genes <- rownames(spliced)[!rownames(spliced) %in% rownames(unspliced)]

    if(length(unspliced_only_genes) > 0){
      # if any unspliced only genes, fill up the spliced only matrix with empty rows with those genes
      unspliced_only_mat <- as(matrix(0, nrow = length(unspliced_only_genes), ncol = num_cells), "sparseMatrix")
      rownames(unspliced_only_mat) <- unspliced_only_genes
      spliced <- rbind(spliced, unspliced_only_mat)
    }

    if(length(spliced_only_genes) > 0){
      # if any spliced only genes, fill up the unspliced only matrix with empty rows with those genes
      spliced_only_mat <- as(matrix(0, nrow = length(spliced_only_genes), ncol = num_cells), "sparseMatrix")
      rownames(spliced_only_mat) <- spliced_only_genes
      unspliced <- rbind(unspliced, spliced_only_mat)
    }

    # create list of assays to use as input to create SCE object
    assay_list <- list(counts = unspliced,
                       spliced = spliced)

  } else {

    # no intron mode then just collapse
    # spliced <- collapse_intron_counts(counts, which_counts = c("spliced"))
    assay_list <- list(counts = counts)

  }

  if (round_counts){

    assay_list <- lapply(assay_list, round)

  }

  sce <- SingleCellExperiment(assays = assay_list)

  return(sce)

}
