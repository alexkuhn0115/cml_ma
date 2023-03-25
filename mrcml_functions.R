
# Requires library(TwoSampleMR)

# Add extra arguments for interval functions? Or inheritance? What is 
# the syntax?

# This is based on the data processing states described in Xue et al. 2021
exposure_outcome_processing <- function(exp_id, out_id, final_data = FALSE, ...) {
  
  # Extract IVs for exposure
  exposure_dat <- extract_instruments(exp_id)
  
  # LD-clumping of extracted IVs
  exposure_dat <- clump_data(exposure_dat)
  
  # Extract summary statistics of outcome, do not use proxy
  outcome_dat <- extract_outcome_data(exposure_dat$SNP, out_id, proxies = 0)
  
  # Harmonize data: correct strand for non-palindromic SNPs, drop all palindromic SNPs
  dat <- harmonise_data(exposure_dat, outcome_dat, action = 3)
  final_dat <- dat[dat$mr_keep,]
  
  # Get final estimates and standard errors for exposure and outcome
  b_exp <- final_dat$beta.exposure
  b_out <- final_dat$beta.outcome
  se_exp <- final_dat$se.exposure
  se_out <- final_dat$se.outcome
  
  if (final_data) {
    return(list(b_exp = b_exp, b_out = b_out, se_exp = se_exp, se_out = se_out, 
                final_dat = final_dat))
  }
  
  return(list(b_exp = b_exp, b_out = b_out, se_exp = se_exp, se_out = se_out))
}
