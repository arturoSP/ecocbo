#' Calculate Simulated Effect Sizes
#' @keywords internal
#' @noRd

sim_es_nested <- function(
  data,
  steps,
  type,
  Sest.method,
  cases,
  N,
  M,
  n,
  m,
  k,
  transformation,
  method,
  dummy,
  useParallel,
  model,
  jitter.base
) {
  # get values for size limits
  data[, 1] <- as.factor(data[, 1])
  data[, 2] <- as.factor(data[, 2])
  factSect <- data[, 1]
  nSect <- nlevels(factSect) # number of treatments
  nivel <- levels(factSect)

  ## Simulated data for H0
  # Create a list to store the results
  ListSim0 <- vector(mode = "list", length = nSect)
  names(ListSim0) <- nivel

  ListParam0 <- vector(mode = "list", length = nSect)
  names(ListParam0) <- nivel
  ListParamA <- vector(mode = "list", length = nSect)
  names(ListParamA) <- nivel

  # helper: format one simulated case list

  format_sim_nested <- function(dataSim, sector_label, dummy = FALSE) {
    out <- lapply(dataSim, function(df) {
      df <- dplyr::mutate(df, sector = as.factor(sector_label))
      new_names <- gsub(
        "unseen\\.species\\s*(\\d+)",
        paste0("unseen.species.", sector_label, ".\\1"),
        names(df)
      )
      names(df) <- new_names
      df
    })

    if (dummy == TRUE) {
      out <- lapply(out, dplyr::mutate, dummy = 1)
    }

    out
  }

  # helper: bind all sectors for each case

  assemble_cases_nested <- function(
    sim_list_by_sector,
    cases,
    nivel,
    dummy = FALSE
  ) {
    out <- vector(mode = "list", length = cases)
    names(out) <- seq_len(cases)

    for (j in seq_len(cases)) {
      tmp_j <- lapply(nivel, function(i) {
        as.data.frame(sim_list_by_sector[[i]][[j]])
      })

      out[[j]] <- dplyr::bind_rows(tmp_j)

      if (dummy == TRUE) {
        out[[j]] <- dplyr::mutate(
          out[[j]],
          dplyr::across(dplyr::everything(), ~ tidyr::replace_na(.x, 0))
        ) |>
          dplyr::relocate(
            c(dplyr::starts_with("unseen"), dummy, sector, sites, N),
            .after = dplyr::last_col()
          )
      } else {
        out[[j]] <- dplyr::mutate(
          out[[j]],
          dplyr::across(dplyr::everything(), ~ tidyr::replace_na(.x, 0))
        ) |>
          dplyr::relocate(
            c(dplyr::starts_with("unseen"), sector, sites, N),
            .after = dplyr::last_col()
          )
      }
    }

    out
  }

  # run the simulation for the different sectors we already have.
  for (i in nivel) {
    # Prepare data by setting replicates to just one value
    dataTrimmedA <- data[data[, 1] == i, -1]
    dataTrimmed0 <- dataTrimmedA
    dataTrimmed0[, 1] <- "zero"

    # Calculate simulation parameters
    dataParameter0 <- assempar(
      dataTrimmed0,
      type = type,
      Sest.method = Sest.method
    )
    dataParameter0$par <- dataParameter0$par |>
      dplyr::mutate(
        Species = sub(
          "unseen\\.species\\s*(\\d+)",
          "unseen.species.\\1",
          Species
        )
      )

    ListParam0[[i]] <- dataParameter0$par[, c(1, 3)]

    dataParameterA <- assempar(
      dataTrimmedA,
      type = type,
      Sest.method = Sest.method
    )
    dataParameterA$par <- dataParameterA$par |>
      dplyr::mutate(
        Species = sub(
          "unseen\\.species\\s*(\\d+)",
          "unseen.species.\\1",
          Species
        )
      )

    ListParamA[[i]] <- dataParameterA

    # Calculate simulated communities
    dataSim0 <- simdata(
      dataParameter0,
      cases = cases,
      N = N,
      sites = M,
      jitter.base = jitter.base
    )

    ListSim0[[i]] <- format_sim_nested(
      dataSim = dataSim0,
      sector_label = i,
      dummy = dummy
    )
  }

  ListParam0 <- pool_fw_across_sectors(ListParam0)

  # Organize the data from different simulations in cases
  # Create a list to store results
  simH0 <- assemble_cases_nested(
    sim_list_by_sector = ListSim0,
    cases = cases,
    nivel = nivel,
    dummy = dummy
  )

  # Calculate species similarity percentage from H0
  ListContribution <- use_simper_nested(simH0)

  ## design and fill the results matrix ----
  NN <- cases * k * (m - 1) * (n - 1)
  resultsHa <- matrix(nrow = NN, ncol = 6)
  colnames(resultsHa) <- c(
    "dat.sim",
    "k",
    "m",
    "n",
    "pseudoFH0",
    "pseudoFHa"
  )

  resultsHa[, 1] <- rep(seq(cases), times = 1, each = (k * (m - 1) * (n - 1)))
  resultsHa[, 2] <- rep(1:k, times = (n - 1) * (m - 1) * cases)
  resultsHa[, 3] <- rep(seq(2, m), times = (n - 1) * cases, each = k)
  resultsHa[, 4] <- rep(seq(2, n), times = cases, each = k * (m - 1))

  ## design the arrays to store the lists ----
  # Factors matrix and data matrix for Ha
  factEnv <- as.data.frame(simH0[[1]][,
    (dim(simH0[[1]])[2] - 2):dim(simH0[[1]])[2]
  ])
  names(factEnv) <- c("sector", "site", "N")

  rm(simH0, ListSim0)
  gc()

  ## Set of parameters for using balancedtwostage ----
  Y <- cbind(1:(N * M))
  YPU <- as.numeric(gl(M, N))
  mm <- resultsHa[, 3]
  nn <- resultsHa[, 3] * resultsHa[, 4]

  Y1 <- cbind(Y, YPU)
  mn <- cbind(mm, nn)

  ## Output container ----
  resultOut <- vector("list", length = NN * (steps + 1))
  pcoaOut <- vector("list", length = NN * (steps + 1))

  cl <- NULL

  if (useParallel) {
    parabar::configure_bar(type = "basic", style = 3)

    cl <- parabar::start_backend(
      cores = parallelly::availableCores() / 2,
      cluster_type = "psock",
      backend_type = "async"
    )

    on.exit(
      {
        if (!is.null(cl)) {
          try(parabar::stop_backend(cl), silent = TRUE)
        }
      },
      add = TRUE
    )

    parabar::export(
      cl,
      variables = c(
        "balanced_sampling_es_nested",
        "safe_balanced_sampling_es_nested",
        "dbmanova_nested",
        "calc_dist_nested"
      ),
      environment = asNamespace("ecocbo")
    )
  }

  for (st in 1:(steps + 1)) {
    paso = st - 1

    stepSim <- vector(mode = "list", length = nSect)
    names(stepSim) <- nivel

    if (st == 1) {
      for (i in nivel) {
        # Calculate simulated communities using Ha param as basis
        dataSimA <- simdata(
          ListParamA[[i]],
          cases = cases,
          N = N,
          sites = M,
          jitter.base = jitter.base
        )

        stepSim[[i]] <- format_sim_nested(
          dataSim = dataSimA,
          sector_label = i,
          dummy = dummy
        )
      }
    } else {
      # Calculate proportion of species to alter in each sector
      alpha = paso / steps

      for (i in nivel) {
        # Hace selección progresiva por ListContribution
        contrib_i <- ListContribution[
          ListContribution$species %in% ListParamA[[i]]$par$Species,
        ]
        adjustN <- min(nrow(contrib_i), ceiling(ListParamA[[i]]$Sest * alpha))
        adjusting <- dplyr::pull(contrib_i[seq_len(adjustN), "species"])

        parHaTemp <- ListParamA[[i]]

        idx_adj <- which(parHaTemp$par$Species %in% adjusting)
        fw_pool_adj <- ListParam0$fw_pool[match(
          parHaTemp$par$Species[idx_adj],
          ListParam0$Species
        )]

        parHaTemp$par[idx_adj, "fw"] <- (1 - alpha) *
          parHaTemp$par[idx_adj, "fw"] +
          alpha * fw_pool_adj

        # Calculate simulated communities using Ha param as basis
        dataSimA <- simdata(
          parHaTemp,
          cases = cases,
          N = N,
          sites = M,
          jitter.base = jitter.base
        )

        stepSim[[i]] <- format_sim_nested(
          dataSim = dataSimA,
          sector_label = i,
          dummy = dummy
        )
      }
    }

    simHaTemp <- assemble_cases_nested(
      sim_list_by_sector = stepSim,
      cases = cases,
      nivel = nivel,
      dummy = dummy
    )

    # Simulation arguments ---
    xH0 <- dim(simHaTemp[[1]])[1]
    yH0 <- dim(simHaTemp[[1]])[2]
    casesHa <- length(simHaTemp)

    HaSim <- array(unlist(simHaTemp), dim = c(xH0, yH0, casesHa))
    HaSim <- HaSim[, 1:(dim(HaSim)[2] - 3), ]

    rm(stepSim, simHaTemp)
    gc()

    if (useParallel) {
      result1 <- parabar::par_lapply(
        cl,
        x = 1:NN,
        fun = safe_balanced_sampling_es_nested,
        Y1,
        mn,
        nSect,
        M,
        N,
        HaSim,
        resultsHa,
        factEnv,
        transformation,
        method,
        model
      )
    } else {
      pb <- txtProgressBar(max = NN, style = 3)
      result1 <- vector("list", length = NN)

      for (i in seq_len(NN)) {
        result1[[i]] <- safe_balanced_sampling_es_nested(
          i = i,
          Y1 = Y1,
          mn = mn,
          nSect = nSect,
          M = M,
          N = N,
          HaSim = HaSim,
          resultsHa = resultsHa,
          factEnv = factEnv,
          transformation = transformation,
          method = method,
          model = model
        )

        setTxtProgressBar(pb, i)
      }

      close(pb)
    }

    idx <- (paso * NN + 1):((paso + 1) * NN)

    resultOut[idx] <- lapply(seq_len(NN), function(i) {
      cur <- result1[[i]]
      data.frame(
        reduction_level = paso / steps,
        step = paso,
        dat_sim = resultsHa[i, "dat.sim"],
        k = resultsHa[i, "k"],
        m = resultsHa[i, "m"],
        n = resultsHa[i, "n"],
        ecological_effect_A = cur$ecological_effect_A,
        ecological_effect_BA = cur$ecological_effect_BA,
        omega2_A = cur$omega2_A,
        omega2_BA = cur$omega2_BA,
        R2_A = cur$R2_A,
        R2_BA = cur$R2_BA,
        pseudoF_A = cur$pseudoF_A,
        pseudoF_BA = cur$pseudoF_BA,
        SS_A = cur$SS_A,
        SS_BA = cur$SS_BA,
        SS_total = cur$SS_total,
        df_A = cur$df_A,
        df_BA = cur$df_BA,
        MS_residual = cur$MS_residual,
        n_groups_A = cur$n_groups_A,
        n_groups_BA = cur$n_groups_BA,
        ok = cur$ok,
        error_message = cur$error_message,
        stringsAsFactors = FALSE
      )
    })

    pcoaOut[idx] <- lapply(seq_len(NN), function(i) {
      cur <- result1[[i]]
      pts <- cur$pcoa_points
      pts$reduction_level <- paso / steps
      pts$step <- paso
      pts
    })

    rm(HaSim, result1)
    gc()

    cat("Step ", paso, ": Done! - ", steps - paso, " remaining")
  }
  if (useParallel && !is.null(cl)) {
    parabar::stop_backend(cl)
    cl <- NULL
  }

  resultOut <- dplyr::bind_rows(resultOut)
  pcoaOut <- dplyr::bind_rows(pcoaOut)

  es_obj <- new_effect_size_data(
    data = resultOut,
    pcoa = pcoaOut,
    call = match.call(),
    model = model
  )

  return(es_obj)
}


###########
###########
#' Calculate Simulated Effect Sizes
#' @keywords internal
#' @noRd

sim_es_single <- function(
  data,
  steps = 10,
  type = "counts",
  Sest.method = "average",
  cases = 5,
  N = 100,
  M = NULL,
  n,
  m = NULL,
  k = 50,
  transformation = "none",
  method = "bray",
  dummy = FALSE,
  useParallel = TRUE,
  model = "single.factor",
  jitter.base = 0.5
) {
  # read data and store it in two objects, one for H0 and one for Ha ----
  datH0 <- data
  datH0[, 1] <- as.factor("zero")
  datHa <- data
  datHa[, 1] <- as.factor(data[, 1])
  a <- nlevels(datHa[, 1])
  M <- a

  ## Helper matrix to store labels ----
  NN <- cases * k * (n - 1)
  resultsHa <- matrix(NA_real_, nrow = NN, ncol = 6)
  colnames(resultsHa) <- c(
    "dat.sim",
    "k",
    "m",
    "n",
    "FobsHa",
    "FobsH0"
  )

  resultsHa[, 1] <- rep(seq(cases), times = 1, each = (k * (n - 1)))
  resultsHa[, 2] <- rep(1:k, times = (n - 1) * cases)
  resultsHa[, 3] <- M
  resultsHa[, 4] <- rep(seq(2, n), times = 1, each = k)

  Y <- cbind(1:(N * M))
  YPU <- as.numeric(gl(M, N))
  mm <- resultsHa[, 3]
  nn <- resultsHa[, 3] * resultsHa[, 4]

  # Calculate H0 simulation parameters
  parH0 <- SSP::assempar(data = datH0, type = type, Sest.method = Sest.method)
  parHa <- SSP::assempar(data = datHa, type = type, Sest.method = Sest.method)

  ## Calculate species similarity percentage from Ha ----
  sppContribution <- use_simper(datHa)

  ## Output container ----
  resultOut <- vector("list", length = NN * (steps + 1))
  pcoaOut <- vector("list", length = NN * (steps + 1))

  simH0 <- SSP::simdata(
    parH0,
    cases = cases,
    N = N,
    sites = M,
    jitter.base = jitter.base
  )

  xH0 <- dim(simH0[[1]])[1]
  yH0 <- dim(simH0[[1]])[2]
  casesHa <- length(simH0)

  H0Sim <- simH0

  if (dummy == TRUE) {
    yH0 <- yH0 + 1
    for (i in seq_len(casesHa)) {
      H0Sim[[i]] <- cbind(simH0[[i]], dummy = 1)
      H0Sim[[i]] <- H0Sim[[i]][, c(1:(yH0 - 3), (yH0), (yH0 - 2):(yH0 - 1))]
    }
  }

  rm(simH0)
  H0Sim <- array(unlist(H0Sim), dim = c(xH0, yH0, casesHa))

  # Loop for species contribution attenuation ----
  n_spp <- dim(sppContribution)[1]
  propSpp <- 1 / steps

  for (st in 1:(steps + 1)) {
    paso = st - 1
    # Calculate number of species to alter
    removing <- paso * propSpp
    adjustN <- ceiling(nrow(sppContribution) * removing)

    if (adjustN != 0) {
      adjusting <- sppContribution[1:adjustN, 1]

      # Obtain H0 parameters from null hypothesis
      zeroParameter <- parH0$par[
        parH0$par$Species %in% adjusting,
        c(1, 3)
      ]

      # Alter the Ha parameters with the corresponding null hypothesis version
      parHaTemp <- parHa
      parHaTemp$par[
        parHaTemp$par$Species %in% adjusting,
        c(3)
      ] <- zeroParameter[, 2]
    } else {
      parHaTemp <- parHa
    }

    simHa <- SSP::simdata(
      parHaTemp,
      cases = cases,
      N = N,
      sites = M,
      jitter.base = jitter.base
    )

    # Simulation arguments ---
    xH0 <- dim(simHa[[1]])[1]
    yH0 <- dim(simHa[[1]])[2]
    casesHa <- length(simHa)

    HaSim <- simHa

    if (dummy == TRUE) {
      yH0 <- yH0 + 1
      for (i in seq_len(casesHa)) {
        HaSim[[i]] <- cbind(simHa[[i]], dummy = 1)
        HaSim[[i]] <- HaSim[[i]][, c(1:(yH0 - 3), (yH0), (yH0 - 2):(yH0 - 1))]
      }
    }

    rm(simHa)

    HaSim <- array(unlist(HaSim), dim = c(xH0, yH0, casesHa))

    # Loop to calculate ecological and inferential effect sizes ----

    if (useParallel) {
      # Registering the cluster of workers with parabar
      parabar::configure_bar(type = "basic", style = 3)
      cl <- parabar::start_backend(
        cores = parallelly::availableCores() / 2,
        cluster_type = "psock",
        backend_type = "async"
      )

      # Exporing functions needed for the parallel iterations
      parabar::export(
        cl,
        variables = c(
          "balanced_sampling_es",
          "dbmanova_oneway",
          "calc_dist"
        ),
        environment = asNamespace("ecocbo")
      )

      # Executing the loop in parallel
      result1 <- parabar::par_lapply(
        cl,
        x = 1:NN,
        fun = balanced_sampling_es,
        Y,
        mm,
        nn,
        YPU,
        HaSim,
        resultsHa,
        transformation,
        method
      )

      parabar::stop_backend(cl)
    } else {
      pb <- txtProgressBar(max = NN, style = 3)
      result1 <- vector("list", length = NN)

      for (i in seq_len(NN)) {
        # Performs the operation iteratively in a for loop
        result1[[i]] <- balanced_sampling_es(
          i,
          Y,
          mm,
          nn,
          YPU,
          HaSim,
          resultsHa,
          transformation,
          method
        )

        # Updating the progress bar
        setTxtProgressBar(pb, i)
      }
      close(pb)
    }

    idx <- (paso * NN + 1):((paso + 1) * NN)

    resultOut[idx] <- lapply(seq_len(NN), function(i) {
      cur <- result1[[i]]
      data.frame(
        reduction_level = paso / steps,
        step = paso,
        dat_sim = resultsHa[i, "dat.sim"],
        k = resultsHa[i, "k"],
        m = resultsHa[i, "m"],
        n = resultsHa[i, "n"],
        ecological_effect = cur$ecological_effect,
        omega2 = cur$omega2,
        R2 = cur$R2,
        pseudoF = cur$pseudoF,
        SS_between = cur$SS_between,
        SS_total = cur$SS_total,
        df_between = cur$df_between,
        MS_residual = cur$MS_residual,
        n_groups = cur$n_groups,
        stringsAsFactors = FALSE
      )
    })

    pcoaOut[idx] <- lapply(seq_len(NN), function(i) {
      cur <- result1[[i]]
      pts <- cur$pcoa_points
      pts$reduction_level <- paso / steps
      pts$step <- paso
      pts
    })

    cat("Step ", paso, ": Done! - ", steps - paso, " remaining")
  }

  resultOut <- dplyr::bind_rows(resultOut)
  pcoaOut <- dplyr::bind_rows(pcoaOut)

  es_obj <- new_effect_size_data(
    data = resultOut,
    pcoa = pcoaOut,
    call = match.call(),
    model = model
  )

  return(es_obj)
}
