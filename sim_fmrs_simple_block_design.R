# Code to simulate an example fMRS dataset for a block design fMRS experiment
# and export a BIDS structure

library(spant)

seq_tr      <- 2    # set the sequence TR
N_scans     <- 448  # just under 15 mins with TR = 2s, but still divisible by 64
bz_inhom_lb <- 4    # Gaussian line-broadening to simulate B0 inhomogeneity
bold_lb_hz  <- 0.0  # linewidth differences in Hz from BOLD T2* effect
ss_spec_snr <- 10   # single shot spectral SNR
subjects    <- 2    # number of subject scans to generate in BIDS format
runs        <- 2    # number of runs to generate in BIDS format
set.seed(1)         # random number generator seed used for noise samples

# Make a data frame containing a single row of basis signal amplitudes.
# Metabolite values are for visual cortex listed in Bednarik et al 2015 Table 1.
# Note Alanine and Glycine are not listed in the table and therefore set to 0.
basis_amps <- data.frame("ala"    = 0.00, "asc"    = 0.96, "asp"   = 3.58,
                         "cr"     = 4.22, "gaba"   = 1.03, "glc"   = 0.62,
                         "gln"    = 2.79, "gly"    = 0.00, "gsh"   = 1.09,
                         "glu"    = 8.59, "gpc"    = 0.54, "ins"   = 6.08,
                         "lac"    = 1.01, "naa"    = 11.9, "naag"  = 1.32,
                         "pch"    = 0.40, "pcr"    = 3.34, "peth"  = 0.93,
                         "sins"   = 0.27, "tau"    = 1.27, "lip09" = 0.00,
                         "lip13a" = 0.00, "lip13b" = 0.00, "lip20" = 0.00,
                         "mm09"   = 4.00, "mm12"   = 4.00, "mm14"  = 4.00,
                         "mm17"   = 4.00, "mm20"   = 4.00)

# Duplicate the row N_scans times to make a table of values
basis_amps <- basis_amps[rep(1, N_scans),]

# simulate two 120 second blocks of stimulation starting at 100 and 500 seconds
onsets    <- c(100, 500)
durations <- rep(120, 2)

# generate a dummy mrs_data object for generating the regressors
mrs_data_dummy <- sim_zero(dyns = N_scans) |> set_tr(seq_tr) |>
                  set_Ntrans(N_scans)

# generate metabolite response functions assuming simple trapezoidal shapes
glu_rf <- gen_trap_reg(onsets, durations, mrs_data = mrs_data_dummy,
                       rise_t = 120, fall_t = 150)

lac_rf <- gen_trap_reg(onsets, durations, mrs_data = mrs_data_dummy,
                       rise_t = 120, fall_t = 150)

asp_rf <- gen_trap_reg(onsets, durations, mrs_data = mrs_data_dummy,
                       rise_t = 120, fall_t = 150)

glc_rf <- gen_trap_reg(onsets, durations, mrs_data = mrs_data_dummy,
                       rise_t = 120, fall_t = 150)

bold_rf <- gen_bold_reg(onsets, durations, mrs_data = mrs_data_dummy)

# update the amplitudes according to predicted changes from Bednarik et at 2015
# Table 1
lac_perc_change <-  29.6
glu_perc_change <-  3.3
glc_perc_change <- -16.0
asp_perc_change <- -5.4

# update metabolite data frame to have dynamic metabolite values
basis_amps$glu <- basis_amps$glu * (glu_rf$stim * glu_perc_change / 100 + 1)
basis_amps$lac <- basis_amps$lac * (lac_rf$stim * lac_perc_change / 100 + 1)
basis_amps$asp <- basis_amps$asp * (asp_rf$stim * asp_perc_change / 100 + 1)
basis_amps$glc <- basis_amps$glc * (glc_rf$stim * glc_perc_change / 100 + 1)

# simulate a typical basis for TE=28ms semi-LASER acquisition at 3T
acq_paras <- def_acq_paras(ft = 127.8e6)
basis     <- sim_basis(names(basis_amps), pul_seq = seq_slaser_ideal,
                       TE1 = 0.008, TE2 = 0.011, TE3 = 0.009)

# apply basis amplitudes to the basis set to generate a simulated fMRS dataset
mrs_dyn_orig <- basis2dyn_mrs_data(basis, basis_amps, seq_tr)

# broaden basis to simulate B0 inhomogeneity, apply any addition BOLD related 
# narrowing
bold_lb_dyn <- (1 - bold_rf$stim_bold) * bold_lb_hz
mrs_dyn     <- mrs_dyn_orig |> lb(bz_inhom_lb) |> lb(bold_lb_dyn, 0) 

# duplicate the data to generate multiple subjects and runs with different noise
# samples
mrs_dyn_list <- rep(list(mrs_dyn), subjects * runs)

# add noise
mrs_dyn_list <- mrs_dyn_list |> add_noise_spec_snr(ss_spec_snr)

# export to BIDS structure
mrs_data_list2bids(mrs_dyn_list, "~/fmrs_block_bids", runs = runs)