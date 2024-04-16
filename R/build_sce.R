#' Build the SCE object from the counts matrix
#'
#' @param counts Counts matrix with rownames corresponding to gene names and colnames corresponding to cell barcodes.
#'   Can be a counts matrix with intron counts specified by -I or an alevin-fry "USA" matrix,
#'   with intron counts marked by "-U" and ambiguous counts "-A".
#' @param include_unspliced Whether or not to include the unspliced reads in the counts matrix.
#'   If TRUE, the main "counts" assay will contain unspliced reads and spliced reads and an additional "spliced"
#'   assay will contain spliced reads only. If TRUE, requires that data has been aligned to a reference containing
#'   spliced and unspliced reads.
#'   Default is TRUE.
#' @param round_counts Logical indicating in the count matrix should be rounded to integers on import.
#'   Default is TRUE.
#'
#' @return SingleCellExperiment object. If `include_unspliced` is TRUE (default), the counts assay will contain
#'   both spliced and unspliced counts and the spliced assay will contain the counts for just the spliced cDNA.
#'   If `include_unspliced` is FALSE, the counts assay will contain spliced cDNA counts only.
#'
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
#'   include_unspliced = TRUE
#' )
#' }
build_sce <- function(counts,
                      include_unspliced = TRUE,
                      round_counts = TRUE) {
  if (!is.logical(include_unspliced)) {
    stop("include_unspliced must be set as TRUE or FALSE")
  }

  if (!is.logical(round_counts)) {
    stop("round_counts must be set as TRUE or FALSE")
  }

  # define if counts matrix has any unspliced reads
  has_unspliced <- any(grep("-[IU]$", rownames(counts)))

  if (include_unspliced && !has_unspliced) {
    stop("No counts corresponding to intronic reads detected.
          If `include_unspliced` is TRUE a reference with spliced and unspliced reads must be used.")
  }

  # if has unspliced data get counts for both unspliced and spliced
  if (include_unspliced && has_unspliced) {
    total <- collapse_intron_counts(counts, which_counts = c("total"))
    spliced <- collapse_intron_counts(counts, which_counts = c("spliced"))

    # before creating the SCE object we need to make sure that the dimensions of
    # spliced and unspliced counts matrix match up
    # first get spliced only genes, unspliced only genes and common genes
    spliced_genes <- rownames(spliced)
    total_genes <- rownames(total)
    common_genes <- intersect(spliced_genes, total_genes)
    spliced_only_genes <- setdiff(spliced_genes, total_genes)
    total_only_genes <- setdiff(total_genes, spliced_genes)

    # build similar matrices that will all have common genes, spliced only genes, unspliced only genes
    spliced <- rbind(
      spliced[common_genes, ],
      spliced[spliced_only_genes, ],
      Matrix::Matrix(0,
        nrow = length(total_only_genes),
        ncol = ncol(spliced),
        dimnames = list(total_only_genes, colnames(spliced)),
        sparse = TRUE, doDiag = FALSE
      )
    )

    total <- rbind(
      total[common_genes, ],
      Matrix::Matrix(0,
        nrow = length(spliced_only_genes),
        ncol = ncol(total),
        dimnames = list(spliced_only_genes, colnames(total)),
        sparse = TRUE, doDiag = FALSE
      ),
      total[total_only_genes, ]
    )

    # create list of assays to use as input to create SCE object
    assay_list <- list(
      counts = total,
      spliced = spliced
    )
  } else if (!include_unspliced && has_unspliced) {
    # still aligned to introns, but want to collapse and just return spliced
    spliced <- collapse_intron_counts(counts, which_counts = c("spliced"))
    assay_list <- list(counts = spliced)
  } else {
    # aligned to spliced only so just return the counts, no introns to collapse
    assay_list <- list(counts = counts)
  }

  if (round_counts) {
    assay_list <- lapply(assay_list, round)
  }

  sce <- SingleCellExperiment(assays = assay_list)

  return(sce)
}
