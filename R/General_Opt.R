#' @name General_Opt
#' @title General Optimizer Routine (nesting many). See example notation by running data("spec")
#' @description Computes a chosen portfolio optimizer with full investment and weight constraints.
#' @details The argument \code{sigma} is a covariance matrix.
#' The MV solution is calculated using \code{quadprog}, others use NLOPTR
#' @param Optim Optimization techniques used: posibilities: c("All", "MV", "MaxSharpe", "MinVol", "MaxDiv", "ERC"). All binds rows if feasible.
#' @param sigma a \eqn{(N \times N)}{(N x N)} covariance matrix.
#' @param mu a \eqn{(N \times 1)}{(N x 1)} vector of estimated returns.
#' @param UB scalar or \eqn{(N\times 1)}{(N x 1)} vector of upper bound weight constraint.
#' @param LB scalar or \eqn{(N\times 1)}{(N x 1)} vector of lower bound weight constraint.
#' @param groups vector of group IDs. The names of the vector must be identical to the asset names.
#' @param group.UB scalar or \eqn{(N_groups\times 1)}{(N_groups x 1)} vector of upper bound group constraints.
#' @param group.LB scalar or \eqn{(N_groups\times 1)}{(N_groups x 1)} vector of lower bound group constraints.
#' @param gamma risk aversion parameter. Higher means more risk averse - this is multiplied by cov. Default: \code{gamma = 0}.
#' @param rf optional risk free rate (annualised).
#' @param maxeval used in nloptr.
#' @param xtol_rel used in nloptr.
#' @param ftol_abs used in nloptr.
#' @param tol used in nloptr.
#' @return A \eqn{(N \times 1)}{(N x 1)} vector of optimal portfolio weights.
#' @author Nico Katzke
#' @examples
#' data("spec")
#' spec$n -> n; spec$gamma -> gamma; spec$sigma -> sigma; spec$mu -> mu; spec$group.LB->group.LB; spec$group.UB -> group.UB;
#' spec$groups_mat -> groups_mat; spec$groups -> groups; spec$LB -> LB; spec$UB -> UB; spec$assets -> assets; rf = 0;
#' maxeval = 5000; xtol_rel = 1e-10; ftol_abs = 1e-10; tol = 1e-10
#' all_opts <- General_Opt( Optim = "All", sigma, mu = mu, UB = UB, LB = LB, groups = groups, group.UB = group.UB, group.LB = group.LB, groups_mat = groups_mat, gamma = gamma)
#' MV_Roll_gamma <- seq(0.1, 150, length.out = 100) %>% as.list() %>% map_df(~General_Opt( Optim = "MV", sigma, mu = mu, UB = UB, LB = LB, groups = groups, group.UB = group.UB, group.LB = group.LB, groups_mat = groups_mat, gamma = .))
#' @export

General_Opt <- function(
    Optim = c("All", "MV", "MaxSharpe", "MinVol", "MaxDiv", "ERC"),
    sigma,
    mu = NULL,
    UB = NULL,
    LB = NULL,
    groups = NULL,
    group.UB = NULL,
    group.LB = NULL,
    groups_mat = NULL,
    gamma = 0, DivRatio = DR,
    rf = 0,
    maxeval = 5000, xtol_rel = 1e-10, ftol_abs = 1e-10, tol = 1e-10
) {

  library(quadprog)
  library(nloptr)


  if(length(Optim) > 1 | !Optim %in% c("All", "MV", "MaxSharpe", "MinVol", "MaxDiv", "ERC")) stop("Provide valid Optim of either: All, MV, MaxSharpe, MinVol, MaxDiv or ERC)")

  n <- dim(sigma)[1]
  asset_names <- colnames(sigma)

  if (!is.null(mu)) {
    if (length(mu)!=n) {
      stop("Different dimensions implied by 'sigma' and 'mu'")
    }
  }

  # Fetch constraints
  if (is.null(UB)) {
    UB <- rep(1, n)
  } else if (length(UB) == 1) {
    # Check constraint
    if (UB * n < 1) stop("Inconsistent constraint (increase UB)")
    UB <- rep(UB, n)
  } else {
    # Check constraint
    if (length(UB) != n) stop("Inconsistent contraint (incorrect elements in UB)")
    UB <- UB
  }
  if (is.null(LB)) {
    LB <- rep(0, n)
  } else if (length(LB) == 1) {
    # Check constraint
    if (LB * n > 1) stop("Inconsistent constraint (decrease LB)")
    LB <- rep(LB, n)
  } else {
    # Check constraint
    if (length(LB) != n) stop("Inconsistent contraint (incorrect elements in LB)")
    LB <- LB
  }
  # Check constraint
  if (!all(pmax(UB, LB) == UB) || !all(pmin(UB, LB) == LB))
    stop("Inconsistent constraint (UB smaller than LB)")


  # GROUPS


  if (!is.null(groups)) {

    if(!is.null(groups_mat) & !is.null(colnames(groups_mat))) groups <- colnames(groups_mat)
    n_groups <- length(unique(groups))

    # Fetch constraints
    if (is.null(group.UB)) {
      group.UB <- rep(1, n_groups)
      names(group.UB) <- unique(groups)
    } else if (length(group.UB) == 1) {
      # Check constraint
      if (group.UB * n_groups < 1) stop("Inconsistent constraint (increase group.UB)")
      group.UB <- rep(group.UB, n_groups)
      names(group.UB) <- unique(groups)
    } else {
      # Check constraint
      if (length(group.UB) != n_groups) stop("Inconsistent contraint (incorrect elements in group.UB)")
      group.UB <- group.UB
    }
    if (is.null(group.LB)) {
      group.LB <- rep(0, n_groups)
      names(group.LB) <- unique(groups)
    } else if (length(group.LB) == 1) {
      # Check constraint
      if (group.LB * n_groups > 1) stop("Inconsistent constraint (decrease group.LB)")
      group.LB <- rep(group.LB, n_groups)
      names(group.LB) <- unique(groups)
    } else {
      # Check constraint
      if (length(group.LB) != n_groups) stop("Inconsistent contraint (incorrect elements in group.LB)")
      group.LB <- group.LB
    }

    if (!all(groups %in% names(group.UB)) | !all(groups %in% names(group.UB))) stop("Inconsistent constraint (missing group names in 'group.UB' or 'group.LB')")

    # Reordering messes with Amat later. remove:
    # group.UB <- group.UB[unique(groups)]
    # group.LB <- group.LB[paste0("LB_", unique(groups))]

    if (!all(pmax(group.UB, group.LB) == group.UB) || !all(pmin(group.UB, group.LB) == group.LB))
      stop("Inconsistent constraint (group.UB smaller than group.LB)")

    #Previous - deprecated, as this does not allow securities to be in different groups at the same time.
    # E.g. a Equity and Global label to an aset class
    #         groups_mat <- sapply(unique(groups), function(x) x==groups)
    #         groups_mat <- cbind(-groups_mat, groups_mat)

    if(is.null(groups_mat)) {
      groups_mat <- sapply(unique(groups), function(x) x==groups)
      groups_mat <- cbind(groups_mat, -groups_mat)[, 1:ncol(groups_mat)]
    }

    # groups_mat <- cbind(-groups_mat, groups_mat)

  } else {

    groups_mat <- NULL
  }


  if (all(dim(sigma) == 1)) { return(1) }

  if (!is.null(mu)) {
    dvec <- mu
  } else {
    dvec <- rep(0, n)
  }


  # Auxiliary helper functions:
  .build_constraints <- function(n, dvec, UB, LB, groups_mat, group.UB, group.LB) {

    # Amat columns: budget | return fence | -I (UB) | +I (LB) | -groups (UB) | +groups (LB)
    Amat <- cbind(
      1,                     # budget equality
      -diag(n),              # upper bounds  (Amat'w >= -UB  =>  w <= UB)
      diag(n),               # lower bounds  (Amat'w >= LB)
      -groups_mat,           # group UB      (Amat'w >= -group.UB)
      groups_mat             # group LB      (Amat'w >= group.LB)
    )

    bvec <- c(
      1,                     # sum(w) = 1  (equality, meq = 1)
      -UB,
      LB,
      -group.UB,
      group.LB
    )

    list(Amat = Amat, bvec = bvec)
  }
  # -----------------------------------------------------------------------------
  # Internal helper: build Amat / bvec WITHOUT the return fence column
  # -----------------------------------------------------------------------------
  .build_constraints_no_ret <- function(n, UB, LB, groups_mat, group.UB, group.LB) {

    Amat <- cbind( 1, -diag(n), diag(n), -groups_mat, groups_mat)
    bvec <- c( 1,-UB,LB,-group.UB,group.LB)

    list(Amat = Amat, bvec = bvec)
  }

  # Aux weight constraint meet checker:
  optim_w_checkr <- function(opt_weights, groups_mat, asset_names, group.UB, group.LB){
    if(!is.null(opt_weights)) {names(opt_weights) <- asset_names}
    if(!is.null(opt_weights) & !is.null(groups_mat)) {
      if(!all(t(groups_mat) %*% opt_weights <= group.UB+1e-8)) stop("Violation of UB")
      if(!all(t(groups_mat) %*% opt_weights >= group.LB-1e-8)) stop("Violation of LB")
    }
  }

  risk_contribs <- function(w) {
   ( w * as.numeric(sigma %*% w) )/ (as.numeric(t(w) %*% sigma %*% w) + 1e-10)
  }


  # ------------------------------------------------------------
  # Warm starts for initial w0
  # ------------------------------------------------------------
  if(!is.null(groups_mat)){
    Amat_warm <- cbind(rep(1, n), -diag(n), diag(n), -groups_mat,groups_mat)
    bvec_warm <- c(1,-UB,LB,-group.UB,group.LB)
  } else {
    Amat_warm <- cbind(rep(1, n), -diag(n), diag(n))
    bvec_warm <- c(1,-UB,LB)
  }

  mv_warm <- tryCatch(
    quadprog::solve.QP(
      Dmat = 2 * sigma,
      dvec = rep(0, n),
      Amat = Amat_warm,
      bvec = bvec_warm,
      meq  = 1
    )$solution,
    error = function(e) rep(1/n, n)
  )

  Res <- list()





  if(Optim == "All" | Optim == "MV"){

    safeOpt <- purrr::safely(quadprog::solve.QP)

    # max μ′w − gamma x w′Σ w
    # gamma applied to covariance - so sensitivity is to cov, not return capping.
    # Gamma now a proper Risk aversion parameter, as opposed to a coinstraint on returns.
    # Small gamma → aggressive / return-seeking, big gamma: conservative
    # Thus: U(w) = μ′w − γ x σ x 2(w)

    con  <- .build_constraints(n, dvec, UB, LB, groups_mat, group.UB, group.LB)

  if (gamma <= 0) stop("gamma must be > 0")
      # Dmat <- 2 * gamma * sigma
      # This normalized version is numerically more stable and equivalent:
    Dmat <- 2 * sigma
    dvec_use <- dvec / gamma

    # This is gives us: max (mu′w − gamma w ′ Σ w)

    # QP minimises 0.5 x'Dmat x - dvec'x:

    opt <- safeOpt(Dmat      = Dmat,
                   dvec      = dvec_use,
                   Amat      = con$Amat,
                   bvec      = con$bvec,
                   meq       = 1)

    if(!is.null(opt$error)) {

      warning("\n\nNo convergence for MV Portfolio...\n\n{opt$error}\n\n")
      if(Optim == "All") Res$MV <- NULL
      result <- NULL

    } else {

    opt_weights  <- opt$result$solution

    optim_w_checkr(opt_weights,  groups_mat, asset_names, group.UB, group.LB)

    w <- opt_weights
    ret <- if (!is.null(dvec)) as.numeric(dvec %*% w) else NA
    vol <- sqrt(as.numeric(t(w) %*% sigma %*% w))
    rc  <- risk_contribs(w)
    asset_vol <- sqrt(diag(sigma))
    DR  <- sum(w * asset_vol) / (vol + 1e-10)

    result  <-  tibble::tibble(Assets = asset_names, Holdings = opt_weights,
                               Exp_return   = ret,
                               Exp_volatility = vol,
                               optimizer = "MV",
                               sharpe = (ret - rf) / (vol + 1e-10),
                               rf = rf,
                               risk_contribs    = setNames(rc, if (!is.null(dvec)) names(dvec) else paste0("A", seq_len(n))),
                               rc_dispersion    = sd(rc),
                               gamma = gamma, DivRatio = DR
    )



    if(Optim == "All" & !is.null(opt$error)) { Res$MV <- result }

  }

  }


  if(Optim == "All" | Optim == "MaxSharpe"){

    # =============================================================================

    excess <- dvec - rf

    # fallback if QP warm start fails
    if (is.null(mv_warm)) {
      # start from LB and distribute remaining weight across available headroom
      w0 <- LB
      w0 <- mv_warm / sum(mv_warm)
      residual <- 1 - sum(w0)

      if (residual < -1e-12) {
        stop("Infeasible: sum(LB) > 1, so no valid starting point exists.")
      }

      room <- UB - LB
      if (sum(room) < residual - 1e-12) {
        stop("Infeasible: sum(UB) < 1, so no valid starting point exists.")
      }

      if (residual > 0) {
        if (sum(room) > 0) {
          w0 <- w0 + residual * room / sum(room)
        }
      }
    } else {
      w0 <- mv_warm
    }

    # ------------------------------------------------------------
    # 2) x0 fix: ensure strictly within lb/ub numerically
    # ------------------------------------------------------------
    #
    # This is the fix for:
    # Error in is.nloptr(ret) : at least one element in x0 < lb
    #
    LB_adj <- LB + tol
    UB_adj <- UB - tol

    if (any(LB_adj > UB_adj)) {
      # If bounds are too tight for tol, fall back to original bounds
      LB_adj <- LB
      UB_adj <- UB
    }

    w0 <- pmax(LB_adj, pmin(UB_adj, w0))

    # Renormalize carefully if needed
    if (abs(sum(w0) - 1) > 1e-10) {
      # push weights back to sum to 1 without violating bounds too much
      for (iter in 1:200) {
        s <- sum(w0)
        if (abs(s - 1) < 1e-12) break

        if (s < 1) {
          room_up <- UB_adj - w0
          if (sum(room_up) <= 1e-14) break
          w0 <- w0 + (1 - s) * room_up / sum(room_up)
        } else {
          room_dn <- w0 - LB_adj
          if (sum(room_dn) <= 1e-14) break
          w0 <- w0 - (s - 1) * room_dn / sum(room_dn)
        }
        w0 <- pmax(LB_adj, pmin(UB_adj, w0))
      }
    }

    # final defensive check
    if (any(w0 < LB_adj - 1e-12)) {
      stop("Initial point x0 still violates lower bounds after adjustment.")
    }
    if (any(w0 > UB_adj + 1e-12)) {
      stop("Initial point x0 still violates upper bounds after adjustment.")
    }

    # ------------------------------------------------------------
    # objective: maximize Sharpe => minimize negative Sharpe
    # ------------------------------------------------------------
    obj_fn <- function(w) {
      er <- sum(excess * w)
      pv <- as.numeric(t(w) %*% sigma %*% w)

      if (!is.finite(pv) || pv <= 0) return(1e12)
      if (!is.finite(er)) return(1e12)
      -er / sqrt(pv)
    }

    grad_fn <- function(w) {
      Sw <- as.vector(sigma %*% w)
      er <- sum(excess * w)
      pv <- as.numeric(t(w) %*% sigma %*% w)

      if (!is.finite(pv) || pv <= 0) {
        return(rep(0, n))
      }

      psd <- sqrt(pv)

      # grad of -(er / sqrt(pv))
      # = -( excess / sqrt(pv) - er * (Sigma w) / (pv^(3/2)) )
      -(excess / psd - er * Sw / (pv^(3/2)))
    }

    # ------------------------------------------------------------
    # equality: sum(w) = 1
    # ------------------------------------------------------------
    heq_fn <- function(w) {
      sum(w) - 1
    }

    heq_jac <- function(w) {
      rep(1, n)
    }

    # ------------------------------------------------------------
    # inequalities for raw nloptr MUST be <= 0
    # ------------------------------------------------------------
    #
    # group lower: group.LB <= G'w  ->  group.LB - G'w <= 0
    # group upper: G'w <= group.UB  ->  G'w - group.UB <= 0
    #
    g_ineq_fn <- function(w) {

      if (!is.null(groups_mat)) {
        gw <- as.numeric(t(groups_mat) %*% w)
      } else {
        gw <- numeric(0)
      }


      if (!is.null(groups_mat)) {
        c(group.LB - gw, gw - group.UB)
      } else {
        numeric(0)
      }

    }


    g_ineq_jac <- function(w) {
      if (!is.null(groups_mat)) {
        rbind(-t(groups_mat), t(groups_mat))
      } else {
        matrix(0, 0, length(w))
      }
    }


    # ------------------------------------------------------------
    # solve with SLSQP
    # ------------------------------------------------------------

    safenloptr <- purrr::safely(nloptr::nloptr)
    sol <- safenloptr(
      x0              = w0,
      eval_f          = obj_fn,
      eval_grad_f     = grad_fn,
      lb              = LB_adj,
      ub              = UB_adj,
      eval_g_eq       = heq_fn,
      eval_jac_g_eq   = heq_jac,
      eval_g_ineq     = g_ineq_fn,
      eval_jac_g_ineq = g_ineq_jac,
      opts = list(
        algorithm   = "NLOPT_LD_SLSQP",
        maxeval     = maxeval,
        xtol_rel    = xtol_rel,
        ftol_abs    = ftol_abs,
        print_level = 0
      )
    )

    if(!is.null(sol$error) ){
      warning("\n\nNo convergence for Max Sharpe Portfolio...\n\n")
      if(Optim == "All") Res$Shrp <- NULL
      result <- NULL
    } else {

    w   <- sol$result$solution
    optim_w_checkr(opt_weights = w,  groups_mat, asset_names, group.UB, group.LB)

    ret <- if (!is.null(dvec)) as.numeric(dvec %*% w) else NA
    vol <- sqrt(as.numeric(t(w) %*% sigma %*% w))
    rc  <- risk_contribs(w)
    asset_vol <- sqrt(diag(sigma))
    DR  <- sum(w * asset_vol) / (vol + 1e-10)

    result  <-  tibble::tibble(Assets = asset_names, Holdings = w,
                               Exp_return   = ret,
                               Exp_volatility = vol,
                               sharpe = (ret - rf) / (vol + 1e-10),
                               rf = rf,
                               rc_dispersion    = sd(rc),
                               risk_contribs    = setNames(rc, if (!is.null(dvec)) names(dvec) else paste0("A", seq_len(n))),
                               optimizer = "Max_Sharpe",
                               gamma = 0, DivRatio = DR
    )

    if(Optim == "All") Res$Shrp <- result

    }

  }


  if(Optim == "All" | Optim == "MinVol"){

    # =============================================================================
    # Solves:  min  w' Sigma w
    #   s.t.   sum(w) = 1,  LB <= w <= UB,  group bounds
    # =============================================================================

      Dmat <- 2 * sigma
      dqp  <- rep(0, n)

      con  <- .build_constraints_no_ret(n, UB, LB, groups_mat, group.UB, group.LB)

      safeOpt <- purrr::safely(quadprog::solve.QP)

      sol  <- safeOpt(Dmat = Dmat,
                       dvec = dqp,
                       Amat = con$Amat,
                       bvec = con$bvec,
                       meq  = 1)

      if(!is.null(sol$error) ){
        warning("\n\nNo convergence for Minvol Portfolio...\n\n")
        if(Optim == "All") Res$Minvol <- NULL
        result <- NULL

      } else {

    w   <- sol$result$solution
    optim_w_checkr(opt_weights = w,  groups_mat, asset_names, group.UB, group.LB)

    ret <- if (!is.null(dvec)) as.numeric(dvec %*% w) else NA
    vol <- sqrt(as.numeric(t(w) %*% sigma %*% w))
    rc  <- risk_contribs(w)
    asset_vol <- sqrt(diag(sigma))
    DR  <- sum(w * asset_vol) / (vol + 1e-10)

    result  <-  tibble::tibble(Assets = asset_names, Holdings = w,
                               Exp_return   = ret,
                               Exp_volatility = vol,
                               sharpe = (ret - rf) / (vol + 1e-10),
                               rf = rf,
                               rc_dispersion    = sd(rc),
                               risk_contribs    = setNames(rc, if (!is.null(dvec)) names(dvec) else paste0("A", seq_len(n))),
                               optimizer = "MinVol",
                               gamma = 0, DivRatio = DR
    )

    if(Optim == "All") Res$Minvol <- result

  }

  }


  if(Optim == "All" | Optim == "MaxDiv"){

    # =============================================================================
    # Maximises the diversification ratio:
    #   DR(w) = (w' sigma_diag) / sqrt(w' Sigma w)
    # where sigma_diag = sqrt(diag(Sigma))  (individual asset vols)
    #
    # Equivalent to minimising -DR(w).  Here we are using nloptr SLSQP.
    # =============================================================================

      asset_vol <- sqrt(diag(sigma))

      sigma <- (sigma + t(sigma)) / 2
      diag(sigma) <- diag(sigma) + 1e-10

      groups_mat <- as.matrix(groups_mat)

      # FEASIBLE WARM START (min variance QP)  as with Ma x Sharpe

      LB_adj <- LB + tol
      UB_adj <- UB - tol

      w0 <- pmax(LB_adj, pmin(UB_adj, mv_warm))


      obj_fn <- function(w) {
        port_vol <- sqrt(as.numeric(t(w) %*% sigma %*% w))
        if (port_vol <= 0) return(1e12)
        -sum(w * asset_vol) / port_vol
      }

      grad_fn <- function(w) {
        Sw       <- as.numeric(sigma %*% w)
        port_var <- as.numeric(t(w) %*% sigma %*% w)
        port_vol <- sqrt(port_var)
        wv       <- sum(w * asset_vol)

        if (port_var <= 0) return(rep(0, n))

        -(asset_vol * port_vol - wv * Sw / port_vol) / port_var
      }


      heq_fn  <- function(w) sum(w) - 1
      heq_jac <- function(w) rep(1, n)


      g_ineq_fn <- function(w) {
        if (!is.null(groups_mat)) {
          gw <- as.numeric(t(groups_mat) %*% w)
        } else {
          gw <- numeric(0)
        }


        if (!is.null(groups_mat)) {
          c(group.LB - gw, gw - group.UB)
        } else {
          numeric(0)
        }
      }


      g_ineq_jac <- function(w) {
        if (!is.null(groups_mat)) {
          rbind(-t(groups_mat), t(groups_mat))
        } else {
          matrix(0, 0, length(w))
        }
      }

      sol <- safenloptr(
        x0              = w0,
        eval_f          = obj_fn,
        eval_grad_f     = grad_fn,
        lb              = LB_adj,
        ub              = UB_adj,
        eval_g_eq       = heq_fn,
        eval_jac_g_eq   = heq_jac,
        eval_g_ineq     = g_ineq_fn,
        eval_jac_g_ineq = g_ineq_jac,
        opts = list(
          algorithm = "NLOPT_LD_SLSQP",
          maxeval   = maxeval,
          xtol_rel  = xtol_rel,
          ftol_abs  = ftol_abs,
          print_level = 0
        )
      )

      if(!is.null(sol$error)) {

        warning("\n\nNo convergence for MV Portfolio...\n\n{sol$error}\n\n")
        if(Optim == "All") Res$Maxdiv <- NULL
        result <- NULL

      } else {

        w <- sol$result$solution
        optim_w_checkr(opt_weights = w,  groups_mat, asset_names, group.UB, group.LB)

      vol <- sqrt(as.numeric(t(w) %*% sigma %*% w))
      rc  <- risk_contribs(w)
      ret <- if (!is.null(dvec)) as.numeric(dvec %*% w) else NA
      asset_vol <- sqrt(diag(sigma))
      DR  <- sum(w * asset_vol) / (vol + 1e-10)

      result  <-  tibble::tibble(Assets = asset_names, Holdings = w,
                                 Exp_return   = ret,
                                 Exp_volatility = vol,
                                 sharpe = (ret - rf) / (vol + 1e-10),
                                 rf = rf,
                                 rc_dispersion    = sd(rc),
                                 risk_contribs    = setNames(rc, if (!is.null(dvec)) names(dvec) else paste0("A", seq_len(n))),
                                 optimizer = "MaxDiv",
                                 gamma = 0,
                                 DivRatio = DR
      )

      if(Optim == "All") Res$Maxdiv <- result

  }

  }

    if(Optim == "All" | Optim == "ERC"){

      # =============================================================================
      # Finds weights such that each asset contributes equally to portfolio variance:
      #   RC_i = w_i * (Sigma w)_i / (w' Sigma w)  is equal for all i
      #
      # Objective (Maillard, Roncalli & Teiletche 2010):
      #   min  sum_i sum_j ( RC_i - RC_j )^2
      # Gives us the best feasible ERC approximation
      # =============================================================================


        LB_adj <- LB + tol
        UB_adj <- UB - tol

        w0 <- pmax(LB_adj, pmin(UB_adj, mv_warm))

        obj_fn <- function(w) {
          Sw <- as.numeric(sigma %*% w)
          port_var <- as.numeric(t(w) %*% sigma %*% w)

          if (port_var <= 0) return(1e12)

          rc <- w * Sw / port_var
          mean_rc <- mean(rc)
          sum((rc - mean_rc)^2)
        }


        grad_fn <- function(w) {
          Sw <- as.numeric(sigma %*% w)
          port_var <- as.numeric(t(w) %*% sigma %*% w)

          if (port_var <= 0) return(rep(0, n))

          rc <- w * Sw / port_var
          mean_rc <- mean(rc)
          diff_rc <- rc - mean_rc

          drc_dwi <- (Sw + w * diag(sigma)) / port_var -
            2 * w * Sw^2 / port_var^2

          2 * diff_rc * drc_dwi
        }


        heq_fn  <- function(w) sum(w) - 1
        heq_jac <- function(w) rep(1, n)


        g_ineq_fn <- function(w) {

          if (!is.null(groups_mat)) {
            gw <- as.numeric(t(groups_mat) %*% w)
          } else {
            gw <- numeric(0)
          }


          if (!is.null(groups_mat)) {
            c(group.LB - gw, gw - group.UB)
          } else {
            numeric(0)
          }
        }


        g_ineq_jac <- function(w) {
          if (!is.null(groups_mat)) {
            rbind(-t(groups_mat), t(groups_mat))
          } else {
            matrix(0, 0, length(w))
          }
        }

        opt <- safenloptr(
          x0              = w0,
          eval_f          = obj_fn,
          eval_grad_f     = grad_fn,
          lb              = LB_adj,
          ub              = UB_adj,
          eval_g_eq       = heq_fn,
          eval_jac_g_eq   = heq_jac,
          eval_g_ineq     = g_ineq_fn,
          eval_jac_g_ineq = g_ineq_jac,
          opts = list(
            algorithm = "NLOPT_LD_SLSQP",
            maxeval   = maxeval,
            xtol_rel  = xtol_rel,
            ftol_abs  = ftol_abs,
            print_level = 0
          )
        )


        if(!is.null(opt$error)) {

          warning("\n\nNo convergence for ERC Portfolio...\n\n{opt$error}\n\n")
          if(Optim == "All") Res$ERC <- NULL
          result <- NULL

        } else {

      w   <- opt$result$solution
      optim_w_checkr(opt_weights = w,  groups_mat, asset_names, group.UB, group.LB)

      rc  <- risk_contribs(w)
      vol <- sqrt(as.numeric(t(w) %*% sigma %*% w))
      ret <- if (!is.null(dvec)) as.numeric(dvec %*% w) else NA
      asset_vol <- sqrt(diag(sigma))
      DR  <- sum(w * asset_vol) / (vol + 1e-10)


      result  <-  tibble::tibble(Assets = asset_names, Holdings = w,
                                 Exp_return   = ret,
                                 Exp_volatility = vol,
                                 sharpe = (ret - rf) / (vol + 1e-10),
                                 rf = rf,
                                 optimizer = "ERC",
                                 risk_contribs    = setNames(rc, if (!is.null(dvec)) names(dvec) else paste0("A", seq_len(n))),
                                 rc_dispersion    = sd(rc),      # should be ~0 at optimum
                                 gamma = 0, DivRatio = DR
      )




      if(Optim == "All") Res$ERC <- result

    }

    }


    if(Optim == "All") {
      result <- Res %>% bind_rows()
    }

    result

  }
