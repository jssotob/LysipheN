#' Perform LysipheN data analysis and vis
#'
#' @param data_path String with the path to locate the Lysiphen file(s).
#' @param file_ext csv by default. Also xls or xlsx available.
#' @param threshold Integer of length 1 with the maximum weight change allowed (in grams). Any weight lower or higher than this offset will be detected as outlier and imputed by the immediate previous one.
#'
#' @return A list
#' @export
#'
#' @examples
#' lysiphen(data_path = getwd(), file_ext = "csv", threshold = 100)
#'
lysiphen <- function(data_path, file_ext = "csv", threshold = 100) {
  # Loading data ------------------------------------------------------------

  file_ext <- tolower(file_ext)

  if (file_ext == "csv") {
    data <- lapply(
      list.files(data_path,
        pattern = paste0("\\.", file_ext, "$"),
        full.names = T
      ),
      function(x) read.csv(x, stringsAsFactors = F)
    ) %>%
      do.call(rbind, .)
    cat(paste0(list.files(data_path,
      pattern = paste0("\\.", file_ext, "$")
    ) %>%
      length(), " files were found and loaded", "\n"))
  } else if (file_ext %in% c("xls", "xlsx")) {
    data <- lapply(
      list.files(data_path,
        pattern = paste0("\\.", file_ext, "$"),
        full.names = T
      ),
      function(x) readxl::read_excel(x, sheet = 1)
    ) %>%
      do.call(rbind, .)
    cat(paste0(list.files(data_path,
      pattern = paste0("\\.", file_ext, "$")
    ) %>%
      length(), " files were found and loaded", "\n"))
  } else {
    stop(paste0(
      "The file extension must be either", " '.csv', ",
      " '.xls', ", "or", " '.xlsx'", "\n"
    ))
  }


  # Naming and VPD ----------------------------------------------------------

  data %<>%
    dplyr::mutate(
      VPD_kPa = (610.78 * 2.71828^(ENV_TEMP / (ENV_TEMP + 238.3) * 17.2694) * (1 - ENV_HUM / 100)) / 100
    ) %>%
    dplyr::rename(
      weight_g = PIN1,
      mainboard_T = PIN2,
      TC1_T = PIN3,
      TC2_T = PIN4,
      irrig_wat_ml = PIN5,
      load_cell_raw_value = PIN6,
      irrig_ID = NFILL,
      cumul_irrig_water = VFILL
    )


  # Detecting and imputing outliers -----------------------------------------


  unit_date <- data %>%
    dplyr::group_split(CROP, LOCATION, DATE) %>%
    as.list()

  outliers <- list()

  for (i in 1:length(unit_date)) {
    detected <- list()
    for (j in 2:nrow(unit_date[[i]])) {
      if (unit_date[[i]]$weight_g[j] >= (unit_date[[i]]$weight_g[j - 1] + threshold) ||
        unit_date[[i]]$weight_g[j] <= (unit_date[[i]]$weight_g[j - 1] - threshold)) {
        detected[[i]] <- unit_date[[i]][c(j - 1, j), ] # just to show, but keep j only
        unit_date[[i]]$weight_g[j] <- unit_date[[i]]$weight_g[j - 1]
      }
    }
    outliers[[i]] <- do.call(rbind, detected)
    rm(detected)
  }

  outliers <- DT::datatable(do.call(rbind, outliers), rownames = F) # Output
  rm(i, j)

  # Irrigation --------------------------------------------------------------

  irrigation <- list()

  for (i in 1:length(unit_date)) {
    if (any(unit_date[[i]]$IRRIGATION_STATUS)) {
      rp <- which(unit_date[[i]]$IRRIGATION_STATUS) %>% max()
      unit_date[[i]]$IRRIGATION_STATUS[rp + 1] <- TRUE
      ini <- which(unit_date[[i]]$IRRIGATION_STATUS) %>% min()
      ini <- unit_date[[i]]$weight_g[ini]

      fin <- which(unit_date[[i]]$IRRIGATION_STATUS) %>% max()
      fin <- unit_date[[i]]$weight_g[fin]

      irrigation[[i]] <- unit_date[[i]] %>%
        dplyr::select(DATE, CROP, LOCATION) %>%
        .[1, ] %>%
        dplyr::mutate(
          ini_irrigation = ini,
          fin_irrigation = fin
        )

      rm(rp, ini, fin)
    } else {
      cat(paste0(i, " not any T \n"))
    }
  }

  irrigation <- do.call(rbind, irrigation)
  rm(i)


  # Transpiration -----------------------------------------------------------

  unit_date <- do.call(rbind, unit_date) %>%
    dplyr::group_split(CROP, LOCATION) %>%
    as.list() %>%
    lapply(., function(x) {
      transpiration <- 0
      for (i in 2:nrow(x)) {
        transpiration[i] <- x$weight_g[i - 1] - x$weight_g[i]
      }
      x$Transpiration <- transpiration
      x %<>%
        dplyr::mutate(Transpiration = dplyr::case_when(
          IRRIGATION_STATUS ~ 0,
          Transpiration <= 0 ~ 0,
          TRUE ~ Transpiration
        ))
    })

  unit_date <- lapply(unit_date, function(x) {
    cumu_tr <- 0
    for (s in 2:nrow(x)) {
      cumu_tr[s] <- x$Transpiration[s] + cumu_tr[s - 1]
    }
    x$cumu_tr <- cumu_tr
    return(x)
  })


  # summary table -----------------------------------------------------------


  summary_table <- lapply(unit_date, function(x) {
    s <- x %>%
      dplyr::mutate(DATE = lubridate::ymd(DATE)) %>%
      dplyr::select(DATE)
    l <- s$DATE[nrow(s)] - s$DATE[1]

    a <- data.frame(
      CROP = x$CROP[1],
      GU = x$LOCATION[1],
      cycle_duration = paste0(as.numeric(l), " days"),
      total_transp_ml = round(x$cumu_tr[nrow(x)], 2)
    )
    return(a)
  }) %>%
    do.call(rbind, .)

  summary_table <- DT::datatable(summary_table, rownames = F) # Output

  # Transpiration vs VPD ----------------------------------------------------

  p_vpd_trans <- do.call(rbind, unit_date) %>%
    dplyr::group_split(DATE, CROP, LOCATION) %>%
    as.list() %>%
    lapply(., function(x) {
      m <- x %>%
        dplyr::filter(!IRRIGATION_STATUS)
      a <- m %>%
        dplyr::arrange(desc(VPD_kPa)) %>%
        .[1, ]
      return(a)
    }) %>%
    do.call(rbind, .) %>%
    tidyr::unite(., "date_time", DATE, TIME) %>%
    dplyr::mutate(date_time = lubridate::ymd_hms(date_time)) %>%
    ggplot2::ggplot(., ggplot2::aes(x = VPD_kPa, y = Transpiration, color = LOCATION)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~CROP) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "VPD (kPa)", y = "Transpiration (ml)")


  ### IT LOOKS LIKE IN SOME DAYS THE IRRIGATION OCCURS AT THE MAX VPD. THEREFORE, TRANSPIRATION IS 0
  ### LOW VPD CORRELATED WITH HIGH HUMIDITY


  # plots -------------------------------------------------------------------

  unit_date <- do.call(rbind, unit_date) # Output

  p_env <- unit_date %>%
    dplyr::select(DATE, TIME, LOCATION, ENV_TEMP, ENV_HUM, VPD_kPa) %>%
    tidyr::unite(., "date_time", DATE, TIME) %>%
    dplyr::mutate(date_time = lubridate::ymd_hms(date_time)) %>%
    tidyr::gather(., key = "Variable", value = "value", -date_time, -LOCATION) %>%
    ggplot2::ggplot(ggplot2::aes(x = date_time)) +
    ggplot2::geom_line(mapping = ggplot2::aes(y = value, color = Variable)) +
    ggplot2::scale_colour_manual(
      values = c("blue", "red", "darkgreen"),
      labels = c(
        "Relative Humidity (%)", "Air Temperature (C)",
        "VPD (kPa)"
      )
    ) +
    ggplot2::facet_wrap(~LOCATION, ncol = 5, scales = "free_x") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12)
    ) +
    ggplot2::labs(x = "Date", y = "", title = "Environmental conditions")


  p_gu <- unit_date %>%
    dplyr::filter(!IRRIGATION_STATUS) %>%
    tidyr::unite(., "date_time", DATE, TIME) %>%
    dplyr::mutate(date_time = lubridate::ymd_hms(date_time)) %>%
    ggplot2::ggplot(ggplot2::aes(x = date_time, y = cumu_tr)) +
    ggplot2::geom_line(color = "#669966", linewidth = 1) +
    ggplot2::facet_wrap(~LOCATION, ncol = 5) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Date", y = "Cumulative Transpiration (ml)",
      title = "Cumulative Transpiration (ml) per LysipheN device"
    )


  p_crop <- unit_date %>%
    dplyr::filter(!IRRIGATION_STATUS) %>%
    tidyr::unite(., "date_time", DATE, TIME) %>%
    dplyr::mutate(date_time = lubridate::ymd_hms(date_time)) %>%
    ggplot2::ggplot(ggplot2::aes(x = date_time, y = cumu_tr, color = LOCATION)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_wrap(~CROP) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(
      x = "Date", y = "Cumulative Transpiration (ml)",
      title = "Cumulative Transpiration (ml) per crop"
    )

  # output ------------------------------------------------------------------

  output <- list(
    data = unit_date,
    outliers = outliers,
    summary_table = summary_table,
    env_plot = p_env,
    cum_transp_plots = list(
      device = p_gu,
      crop = p_crop,
      VPD = p_vpd_trans
    )
  )



  return(output)
}
