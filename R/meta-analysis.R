library(metap)

# Documentation for P‑Value Combination Methods using metap
#
# These helper functions combine a vector of individual p‑values into a single
# summary p‑value using different statistical methods implemented in the metap
# package.
#
# General Notes:
# – All functions return a single numeric p‑value.
# – Assumes independence of input tests; violations can inflate type I error.
# – The metap package must be installed (`install.packages("metap")`) and loaded
# – Choice of method depends on expected signal pattern:
#     • Tippett: one very strong effect,
#     • Fisher: many moderately significant effects,
#     • Stouffer: balanced moderate effects,
#     • Wilkinson: flexibly targets a specified order statistic.


#------------------------------------------------------------------------------
# Wilkinson Method (Order‑Statistic Combination)
#
#    Function: wilkinson_combined_pvalue(pvalues)
#    • Input:
#        – pvalues: numeric vector of independent p‑values (0 < p ≤ 1).
#    • Statistical Principle:
#        – Based on the k-th smallest p‑value in the vector.
#        – For k = 1 (the smallest), equivalent to Tippett’s method; for k = n 
#          (the largest), equivalent to the maximum p.
#        – More generally, you can test whether at least k p‑values are below 
#          a threshold by comparing the k-th order statistic to its beta 
#          distribution under the null hypothesis of uniformity.
#    • Operation in metap:
#        1. metap::wilkinsonp() defaults to combining all p‑values via their
#           full order distribution (equally weighting every order-statistic).
#        2. Internally computes the probability of observing the ordered set
#           under the null.
#    • Output:
#        – Combined p‑value reflecting the extremity of the ensemble of p‑values.
#
#------------------------------------------------------------------------------

wilkinson_combined_pvalue <- function(pvalues) {
  wilkinson_result <- metap::wilkinsonp(pvalues)
  combined_pvalue <- wilkinson_result$p
  return(combined_pvalue)}

#------------------------------------------------------------------------------
# Tippett Method (Minimum P‑Value Combination)
#
#    Function: tippett_combined_pvalue(pvalues)
#    • Input:
#        – pvalues: numeric vector of independent p‑values.
#    • Statistical Principle:
#        – Uses only the minimum p‑value (p_min = min(pvalues)).
#        – Under the null, p_min follows Beta(1, n) distribution (for n tests);
#          small p_min indicates at least one strong signal.
#        – Most sensitive to a single strong effect but ignores the rest.
#    • Operation in metap:
#        1. metap::minimump() computes P(min(P_i) ≤ p_min | null) 
#           = 1 − (1 − p_min)^n.
#        2. Returns this adjusted probability as the combined p‑value.
#    • Output:
#        – Combined p‑value that flags if any single p is unexpectedly small.
#
#------------------------------------------------------------------------------

tippett_combined_pvalue <- function(pvalues) {
  tippett_result <- metap::minimump(pvalues)
  combined_pvalue <- tippett_result$p
  return(combined_pvalue)}


#------------------------------------------------------------------------------
# Fisher Method (Sum of Log P‑Values)
#
#    Function: fisher_combined_pvalue(pvalues)
#    • Input:
#        – pvalues: numeric vector of independent p‑values.
#    • Statistical Principle:
#        – Transforms each p_i to −2·ln(p_i), then sums these quantities:
#            X² = −2 ∑ ln(p_i).
#        – Under the null, X² ∼ χ² with 2n degrees of freedom.
#        – Sensitive to moderate signals across many tests.
#    • Operation in metap:
#        1. metap::sumlog() calculates the X² statistic.
#        2. Computes the tail probability P(χ²_2n ≥ X²).
#    • Output:
#        – Combined p‑value reflecting cumulative evidence across all p-values.
#
#------------------------------------------------------------------------------

fisher_combined_pvalue <- function(pvalues) {
  fisher_result <- metap::sumlog(pvalues)
  combined_pvalue <- fisher_result$p
  return(combined_pvalue)}

#------------------------------------------------------------------------------
# Stouffer Method (Sum of Z‑Scores)
#
#    Function: stouffer_combined_pvalue(pvalues)
#    • Input:
#        – pvalues: numeric vector of independent p‑values.
#    • Statistical Principle:
#        – Converts each p_i to a Z‑score via Z_i = Φ⁻¹(1 − p_i) (standard normal
#          quantile).
#        – Computes the weighted sum (here equal weights):
#            Z_combined = (1/√n) ∑ Z_i.
#        – Under the null, Z_combined ∼ N(0,1).
#        – Balances contributions: strong but isolated signals are 
#          down-weighted, while consistent moderate signals accumulate.
#    • Operation in metap:
#        1. metap::sumz() transforms p-values to Z’s and sums them.
#        2. Calculates the survival function P(N(0,1) ≥ Z_combined).
#    • Output:
#        – Combined p‑value reflecting average signal strength across p-values.
#
#------------------------------------------------------------------------------

stouffer_combined_pvalue <- function(pvalues) {
  stouffer_result <- metap::sumz(pvalues)
  combined_pvalue <- stouffer_result$p
  return(combined_pvalue)}
