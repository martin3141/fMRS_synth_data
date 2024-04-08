# Code to simulate an example fMRS dataset for a block design fMRS
# experiment. Only simple line broadening and additive noise are considered.

library(spant)

seq_tr      <- 2   # set the sequence TR
N_scans     <- 448 # just under 15 mins with TR = 2s, but still divisible by 64
basis_lb    <- 4   # Gaussian line-broadening to simulate typical shimming
noise_level <- 10  # frequency domain noise standard deviation added to taste
set.seed(1)        # random number generator seed

# Simulate a typical basis for TE=28ms semi-LASER acquisition at 3T
acq_paras <- def_acq_paras() # field strength could be adjusted here
brain_sim <- sim_brain_1h(acq_paras, full_output = TRUE, TE1 = 0.008,
                          TE2 = 0.011, TE3 = 0.009, pul_seq = seq_slaser_ideal)

# Use more accurate baseline concentration estimates for visual cortex from 
# Bednarik et al 2015 Table 1. Note Ascorbate is not simulated here, but likely
# small enough to ignore at 3T.
amps  <- brain_sim$amps
amps["Ala"]  <- 0.00 # not listed in Bednarik et al, so set to zero
amps["Asp"]  <- 3.58
amps["Cr"]   <- 4.22
amps["GABA"] <- 1.03
amps["Glc"]  <- 0.62
amps["Gln"]  <- 2.79
amps["Gly"]  <- 0.00 # not listed in Bednarik et al, so set to zero
amps["GSH"]  <- 1.09
amps["Glu"]  <- 8.59
amps["GPC"]  <- 0.54
amps["Ins"]  <- 6.08
amps["Lac"]  <- 1.01
amps["NAA"]  <- 6.08
amps["NAAG"] <- 1.32
amps["PCh"]  <- 0.40
amps["PCr"]  <- 3.34
amps["PEth"] <- 0.93
amps["sIns"] <- 0.27
amps["Tau"]  <- 1.27

# From Bednarik et at 2015 Table 1.
lac_perc_change <-  29.6
glu_perc_change <-  3.3
glc_perc_change <- -16.0
asp_perc_change <- -5.4

# generate a dataframe of baseline metabolite level for each dynamic scan
amps_df <- data.frame(t(amps))
amps_df <- amps_df[rep(1, N_scans),]

# generate a dummy mrs_data object to aid metabolite response fn generation
mrs_data_dummy <- brain_sim$mrs_data |> set_tr(seq_tr) |> set_Ntrans(N_scans) |>
                  rep_dyn(N_scans)

# generate metabolite response functions assuming simple trapezoidal shapes with
# 30 second rise time and 120 second fall time
# simulate two 120 second blocks of stimulation starting at 100 and 500 seconds
onsets    <- c(100, 500)
durations <- rep(120, 2)
labels    <- rep("x", 2)

glu_rf <- gen_trap_rf(onsets, durations, labels, mrs_data_dummy,
                      rise_t = 120, fall_t = 150)

lac_rf <- gen_trap_rf(onsets, durations, labels, mrs_data_dummy,
                      rise_t = 120, fall_t = 150)

asp_rf <- gen_trap_rf(onsets, durations, labels, mrs_data_dummy,
                      rise_t = 120, fall_t = 150)

glc_rf <- gen_trap_rf(onsets, durations, labels, mrs_data_dummy,
                      rise_t = 120, fall_t = 150)

# update metabolite dataframe to have dynamic metabolite values
amps_df$Glu <- amps_df$Glu * (glu_rf$x * glu_perc_change / 100 + 1)
amps_df$Lac <- amps_df$Lac * (lac_rf$x * lac_perc_change / 100 + 1)
amps_df$Asp <- amps_df$Asp * (asp_rf$x * asp_perc_change / 100 + 1)
amps_df$Glc <- amps_df$Glc * (glc_rf$x * glc_perc_change / 100 + 1)

# generate a list of spectra based on dynamic metabolite values
mrs_list <- vector(mode = "list", length = N_scans)
for (n in 1:N_scans) {
  amps_n <- amps_df[n,] |> as.numeric()
  mrs_list[[n]] <- basis2mrs_data(brain_sim$basis, sum_elements = TRUE,
                                  amps = amps_n)
}

# convert list of spectra to a single dynamic scan and set timing parameters
mrs_dyn_orig <- append_dyns(mrs_list) |> set_tr(seq_tr) |> set_Ntrans(N_scans) 

# broaden basis add noise
mrs_dyn <- mrs_dyn_orig |> lb(basis_lb) |> add_noise(noise_level)

# plots
mrs_dyn |> lb(4) |> sub_mean_dyns() |> image(xlim = c(4, 0.5))

mrs_dyn |> lb(4) |> mean_dyn_blocks(32) |> sub_mean_dyns() |>
                    stackplot(xlim = c(4, 0.5), y_offset = 110)

# extract the task and rest spectra and plot the mean difference
# generate a Boxcar function to describe the task
# this isn't optimal, but maybe useful in the case where only minimal
# assumptions can be made about the metabolite response function
boxcar_rf <- gen_trap_rf(onsets, durations, labels, mrs_data_dummy, rise_t = 0,
                         fall_t = 0)

plot(boxcar_rf$time, boxcar_rf$x, type = "l")
lines(boxcar_rf$time, glu_rf$x, col = "red")
task_bool <- boxcar_rf$x > .5
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
