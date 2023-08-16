library(pacman)
p_load(
  dplyr, magrittr, lubridate, DT, plotly, RColorBrewer, ggpubr, dygraphs,
  xts, htmltools, broom, readxl
)

data_path <- "G:/My Drive/CIAT_2/Duvan/Exp2/GU_data/new/"
file_ext <- "Csv"



lysiphen <- function(data_path = getwd(), file_ext = "csv") {
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
  }

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


  # Calculating Irrigation --------------------------------------------------

  unit_date <- data %>%
    dplyr::group_split(CROP, LOCATION, DATE) %>%
    as.list()




  return()
}
