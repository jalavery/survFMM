# Clean subgroup labels for survFMM objects A common problem in finite mixture modeling is that subgroup labels are superficial, e.g., the group labeled 'subgroup 1' does not always correspond to the subgroup with the smallest treatment effect. The purpose of this function is to update subgroup labels in order of increasing magnitude of the treatment term such that subgroup 1 always has the smallest treatment term and subgroup 2 always has the largest treatment estimate.

Clean subgroup labels for survFMM objects A common problem in finite
mixture modeling is that subgroup labels are superficial, e.g., the
group labeled 'subgroup 1' does not always correspond to the subgroup
with the smallest treatment effect. The purpose of this function is to
update subgroup labels in order of increasing magnitude of the treatment
term such that subgroup 1 always has the smallest treatment term and
subgroup 2 always has the largest treatment estimate.

## Usage

``` r
clean_subgroup_labels(survFMM_object, tx_term)
```

## Arguments

- survFMM_object:

  Object returned from the survFMM function (may be either AFT-FMM or
  IPCW-FMM model)

- tx_term:

  The name of the variable corresponding to the treatment

## Value

aft_fmm_clean A survFMM object with the same list elements as the input
object, but with subgroup labels ordered by increasing magnitude of the
treatment term
