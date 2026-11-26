taxa_rank_counts <- function(ps, rank) {
  sample_totals <- phyloseq::sample_sums(ps)
  counts <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) counts <- t(counts)
  labels <- rownames(counts)
  if (rank != "ASV") {
    taxonomy <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
    shiny::validate(shiny::need(rank %in% colnames(taxonomy), "The selected taxonomic rank is unavailable."))
    taxonomy <- as(taxonomy, "matrix")[labels, , drop = FALSE]
    rank_labels <- trimws(taxonomy[, rank])
    known <- !is.na(rank_labels) & nzchar(rank_labels)
    # Aggregate by lineage and retain each unclassified feature separately.
    lineage <- as.data.frame(taxonomy[, seq_len(match(rank, colnames(taxonomy))), drop = FALSE])
    groups <- vctrs::vec_group_id(lineage)
    groups[!known] <- nrow(counts) + seq_len(sum(!known))
    representative <- !duplicated(groups)
    counts <- rowsum(counts, groups, reorder = FALSE)
    rownames(counts) <- labels[representative]
    labels <- ifelse(known, rank_labels, labels)[representative]
    duplicated_labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
    labels[duplicated_labels] <- paste0(labels[duplicated_labels], " [", rownames(counts)[duplicated_labels], "]")
    labels <- make.unique(labels)
  }
  list(counts = counts, labels = labels, sample_totals = sample_totals)
}
