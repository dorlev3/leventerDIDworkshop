#' Simulate a mean shift DID data-set
#'
#' @param nobs Number of observations
#' @param ngroups Number of groups
#' @param group_treat_order Order of group treatment. If groups have same treatment timing, then omit from argument. For example, if three groups (ngroups = 3), and groups 1 and 2 are both treated at same time, input is: 'group_treat_timing = c(1,3)'
#' @param group_treat_time Timing each group gets treated. From above example, if inputted 'group_treat_timing = c(1,3)' -- i.e. groups 1:2 are treated first and groups 3: are treated after, the need to input two time periods, one for first treatment and one for later treatment. If later are never-treated, input a time period after maximum time period in the data, for example 'group_treat_time=c(1995, 9999)'
#' @param years Minimum and maximum time periods
#' @param tau_vec A vector of treatment effects. For hetergenous data-sets, simply input different numbers, e.g. 'tau_vec = c(1,.5)'
#' @param seed A random seed
#'
#' @returns A simulated data frame for DID analysis.
#' @import data.table
#' @export
#'
#' @examples
#' df <- sim_data_mean_shift()
sim_data_mean_shift <- function(nobs = 300,
                                ngroups = 3,
                                group_treat_order = c(1, 2, 3),
                                years = c(1980, 2010),
                                group_treat_time = c(1990, 2000, 9999),
                                tau_vec = c(1, 1, 1),
                                seed = 1) {

  set.seed(seed)
  dat = CJ(id = 1:nobs, year = years[1]:years[2])[
    , state := sample(1:ngroup, 1), by = "id"][
      , year_fe := rnorm(1, mean = 0, sd = 1), by = "year"][
        , unit_fe := rnorm(1, mean = state, sd = 1), by = "id"]

  setkey(dat, state, id, year)

  treatment_groups = data.table(
    state = group_treat_order,
    cohort = group_treat_time,
    tau = tau_vec)
  dat = treatment_groups[dat, roll = TRUE, on = "state"]

  dat                                                [
    , treat  := as.numeric(year >= cohort)          ][
      , tau_single := fifelse(treat == 1, tau, 0)       ][
        , error  := rnorm(.N, 0, 1)                    ][
          , Y0 := unit_fe + year_fe + error       ][
            , Y1 := Y0 + tau_single ][
              , Y := treat*Y1 + (1-treat)*Y0][
                , time_to_treat := year - cohort                ]

  return(dat)
}


#' Simulate a trend shift DID data-set
#'
#' @param nobs Number of observations
#' @param ngroups Number of groups
#' @param group_treat_order Order of group treatment. If groups have same treatment timing, then omit from argument. For example, if three groups (ngroups = 3), and groups 1 and 2 are both treated at same time, input is: 'group_treat_timing = c(1,3)'
#' @param group_treat_time Timing each group gets treated. From above example, if inputted 'group_treat_timing = c(1,3)' -- i.e. groups 1:2 are treated first and groups 3: are treated after, the need to input two time periods, one for first treatment and one for later treatment. If later are never-treated, input a time period after maximum time period in the data, for example 'group_treat_time=c(1995, 9999)'
#' @param years Minimum and maximum time periods
#' @param tau_vec A vector of treatment effects. For hetergenous data-sets, simply input different numbers, e.g. 'tau_vec = c(1,.5)'
#' @param seed A random seed
#'
#' @returns A simulated data frame for DID analysis.
#' @import data.table
#' @export
#'
#' @examples
#' df <- sim_data_trend_shift()
sim_data_trend_shift <- function(nobs = 300,
                                 ngroups = 3,
                                 group_treat_order = c(1, 2, 3),
                                 years = c(1980, 2010),
                                 group_treat_time = c(1990, 2000, 9999),
                                 tau_vec = c(1, 1, 1),
                                 seed = 1) {

  set.seed(seed)
  dat = CJ(id = 1:nobs, year = years[1]:years[2])[
    , state := sample(1:ngroup, 1), by = "id"][
      , year_fe := rnorm(1, mean = 0, sd = 1), by = "year"][
        , unit_fe := rnorm(1, mean = state, sd = 1), by = "id"]

  setkey(dat, state, id, year)

  treatment_groups = data.table(
    state = group_treat_order,
    cohort = group_treat_time,
    tau = tau_vec)
  dat = treatment_groups[dat, roll = TRUE, on = "state"]

  dat                                                [
    , treat  := as.numeric(year >= cohort)          ][
      , tau_single := fifelse(treat == 1, tau, 0)       ][
        , cumtau := cumsum(tau_single), by = "id"            ][
          , error  := rnorm(.N, 0, 1)                    ][
            , Y0 := unit_fe + year_fe + error       ][
              , Y1 := Y0 + cumtau ][
                , Y := treat*Y1 + (1-treat)*Y0][
                  , time_to_treat := year - cohort                ]

  return(dat)
}
