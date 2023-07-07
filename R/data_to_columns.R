
# ------------------------------------------------------------------------------------------------ #

### Title: Data formating ####
### Author: Patrón-Rivero, C. #
### Project: "Evaluating Porthidium’s ecological niches and phylogenetic relationships: #
### an approach to understanding intra and interspecific variation" ###

# ------------------------------------------------------------------------------------------------ #

# Convert into different formats #

# ------------------------------------------------------------------------------------------------ #

data_to_columns <- function(dir_csv) {
  result <- list()
  data <- read.csv(dir_csv)
  result$data <- separate(data, Test, into = c("Spp1", "Spp2"), sep = " vs ")

  new_data <- data.frame(Spp1 = character(), Spp2 = character(), D = numeric(), p.val_D = numeric(), I = numeric(), p.val_I = numeric(), stringsAsFactors = FALSE)

  for (species in unique(data$Spp1)) {
    new_data <- rbind(new_data, data.frame(Spp1 = species, Spp2 = species, D = 1, p.val_D = 0, I = 1, p.val_I = 0))
  }

  rownames(new_data) <- NULL
  result$new_data <- rbind(data, new_data)

  return(result)
}

# ------------------------------------------------------------------------------------------------ #

# Arguments #

# ------------------------------------------------------------------------------------------------ #

#' @param dir_csv Directory when the csv resulting from sim_spp_test function

#' @return A list with two files: 1) data = Same dataframe but with two columns from the two species tested, and 2) a new_data = new dataframe as data but including the species comparisons with itselfs for matrix visualization or phhylogenetic analysis.

#' @examples
#'
#' dir_csv <- "D/example/result"
#'
#' data_to_columns <- function(dir_csv)
#'
# ------------------------------------------------------------------------------------------------ #

### EndNotRun
