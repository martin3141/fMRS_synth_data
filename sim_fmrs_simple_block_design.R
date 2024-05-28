# Code to simulate an example fMRS dataset for a block design fMRS
# experiment.

library(spant)

seq_tr      <- 2    # set the sequence TR
N_scans     <- 448  # just under 15 mins with TR = 2s, but still divisible by 64
bz_inhom_lb <- 4    # Gaussian line-broadening to simulate B0 inhomogeneity
bold_lb_hz  <- 0.0  # linewidth differences in Hz from BOLD T2* effect
noise_level <- 10   # frequency domain noise standard deviation added to taste
set.seed(1)         # random number generator seed

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

# Duplicate the row N_scans times
basis_amps <- basis_amps[rep(1, N_scans),]

# Update the amplitudes according to predicted dynamics

# From Bednarik et at 2015 Table 1.
lac_perc_change <-  29.6
glu_perc_change <-  3.3
glc_perc_change <- -16.0
asp_perc_change <- -5.4

# simulate two 120 second blocks of stimulation starting at 100 and 500 seconds
onsets    <- c(100, 500)
durations <- rep(120, 2)

# generate a dummy mrs_data object for generating the regressors
mrs_data_dummy <- sim_resonances() |> set_tr(seq_tr) |> set_Ntrans(N_scans) |>
                  rep_dyn(N_scans)

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

# update metabolite data frame to have dynamic metabolite values
basis_amps$glu <- basis_amps$glu * (glu_rf$stim * glu_perc_change / 100 + 1)
basis_amps$lac <- basis_amps$lac * (lac_rf$stim * lac_perc_change / 100 + 1)
basis_amps$asp <- basis_amps$asp * (asp_rf$stim * asp_perc_change / 100 + 1)
basis_amps$glc <- basis_amps$glc * (glc_rf$stim * glc_perc_change / 100 + 1)

# aimulate a typical basis for TE=28ms semi-LASER acquisition at 3T
acq_paras <- def_acq_paras(ft = 127.8e6)
basis     <- sim_basis(names(basis_amps), pul_seq = seq_slaser_ideal,
                       TE1 = 0.008, TE2 = 0.011, TE3 = 0.009)

# apply basis amplitudes to the basis set to generate a simulated fMRS dataset
mrs_dyn_orig <- basis2dyn_mrs_data(basis, basis_amps, seq_tr)

# broaden basis to simulate B0 inhomogeneity, apply any addition BOLD related 
# narrowing and add noise
bold_lb_dyn <- (1 - bold_rf$stim_bold) * bold_lb_hz
mrs_dyn <- mrs_dyn_orig |> lb(bz_inhom_lb) |>  lb(bold_lb_dyn, 0) |>
  add_noise(noise_level)

# plots
mrs_dyn |> lb(4) |> sub_mean_dyns() |> image(xlim = c(4, 0.5))

mrs_dyn |> lb(4) |> mean_dyn_blocks(32) |> sub_mean_dyns() |>
                    stackplot(xlim = c(4, 0.5), y_offset = 110)

# extract the task and rest spectra and plot the mean difference
# generate a Boxcar function to describe the task
# this isn't optimal, but maybe useful in the case where only minimal
# assumptions can be made about the metabolite response function
boxcar_rf <- gen_trap_reg(onsets, durations, mrs_data = mrs_data_dummy,
                          rise_t = 0, fall_t = 0)

plot(boxcar_rf$time, boxcar_rf$stim, type = "l")
lines(boxcar_rf$time, glu_rf$stim, col = "red")
lines(boxcar_rf$time, bold_rf$stim_bold, col = "blue")

task_bool <- boxcar_rf$stim > .5
task_inds <- which(task_bool)
rest_inds <- which(!task_bool)

mean_rest <- mrs_dyn |> get_dyns(rest_inds) |> mean_dyns()
mean_task <- mrs_dyn |> get_dyns(task_inds) |> mean_dyns()

# plot the mean task and rest spectra and subtract
subtracted <- (mean_task - mean_rest)  |> lb(4) |> scale_mrs_amp(50)
list(subtracted, mean_rest, mean_task) |> lb(4) |>
  stackplot(xlim = c(4, 0.5), y_offset = 20, labels = 
                       c("(task-rest) x 50", "rest", "task"),
            mar = c(3.5, 1, 1, 6.5))

# export as nifti MRS
write_mrs(mrs_dyn, "fmrs_block.nii.gz", force = TRUE)