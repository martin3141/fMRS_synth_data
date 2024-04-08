# Code to simulate an example dataset for an event-related design fMRS
# experiment. One second duration stimuli are modeled at 20 second
# intervals. The glutamate response function (GRF) is modeled as a smoothed
# trapezoid lasting approximately 2.5 seconds (Mullins 2018). Only simple line
# broadening and additive noise are considered.

library(spant)

seq_tr      <- 2   # set the sequence TR
N_scans     <- 448 # just under 15 mins with TR = 2s, but still divisible by 64
basis_lb    <- 4   # Gaussian line-broadening to simulate typical shimming
noise_level <- 100 # frequency domain noise standard deviation added to taste
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

# From Yakovlev et al 2022
glu_perc_change <- 14.0

# generate a dataframe of baseline metabolite level for each dynamic scan
amps_df <- data.frame(t(amps))
amps_df <- amps_df[rep(1, N_scans),]

# generate a dummy mrs_data object to aid metabolite response fn generation
mrs_data_dummy <- brain_sim$mrs_data |> set_tr(seq_tr) |> set_Ntrans(N_scans) |>
                  rep_dyn(N_scans)

onsets    <- seq(from = 19.3, to = 860, by = 20)
durations <- rep(1.0, length(onsets))
labels    <- rep("x", length(onsets))

# high temporal resolution GRF plot
gen_trap_rf(onsets, durations, labels, mrs_data_dummy,
          rise_t = 0.25, fall_t = 1.0, match_tr = FALSE, smo_sigma = 0.2) |> 
  plot(type = "l", xlim = c(0, 90), ylab = "GRF")

glu_rf <- gen_trap_rf(onsets, durations, labels, mrs_data_dummy, rise_t = 0.25,
                      fall_t = 1.0, smo_sigma = 0.2)$x

# update metabolite dataframe to have dynamic metabolite values
amps_df$Glu <- amps_df$Glu * (glu_rf * glu_perc_change / 100 + 1)

# generate a list of spectra based on dynamic metabolite values
mrs_list <- vector(mode = "list", length = N_scans)
for (n in 1:N_scans) {
  amps_n <- amps_df[n,] |> as.numeric()
  mrs_list[[n]] <- basis2mrs_data(brain_sim$basis, sum_elements = TRUE,
                                  amps = amps_n)
}

# convert list of spectra to a single dynamic scan and set timing parameters
mrs_dyn_orig <- append_dyns(mrs_list) |> set_tr(seq_tr) |> set_Ntrans(N_scans) 

# broaden basis and add noise
mrs_dyn <- mrs_dyn_orig |> lb(basis_lb) |> add_noise(noise_level)

# plots
mrs_dyn |> lb(4) |> sub_mean_dyns() |> image(xlim = c(4, 0.5))

# get a vector of dynamic indices of "active" spectra
task_bool <- glu_rf > .5

task_inds <- which(task_bool)
rest_inds <- which(!task_bool)

mean_task <- mrs_dyn |> mean_dyns(task_inds)
mean_rest <- mrs_dyn |> mean_dyns(rest_inds)

# plot the mean task and rest spectra and subtract
subtracted <- (mean_task - mean_rest)  |> lb(4) |> scale_mrs_amp(5)
list(subtracted, mean_rest, mean_task) |> lb(4) |> 
  stackplot(xlim = c(4, 0.5), y_offset = 20, labels = 
              c("(task-rest) x 5", "rest", "task"), mar = c(3.5, 1, 1, 6.5))

# export as nifti MRS
write_mrs(mrs_dyn, "fmrs_event_related.nii.gz", force = TRUE)