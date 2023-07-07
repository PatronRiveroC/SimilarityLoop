
# ------------------------------------------------------------------------------------------------ #

### Title: Niche similarity test for multiple species ####
### Author: Patrón-Rivero, C. #
### Terra adaptation: JD Sánchez-Rodriguez #
### Project: "Evaluating Porthidium’s ecological niches and phylogenetic relationships: #
### an approach to understanding intra and interspecific variation" ###

# ------------------------------------------------------------------------------------------------ #

# Similaryty for multiple Species Test #

# ------------------------------------------------------------------------------------------------ #

sim_spp_test <- function(var_dir, spp_dir, dir_save, M_dir, n_back, spp){

  val_1 <- list()
  val_2 <- list()
  val_3 <- list()
  val_4 <- list()
  val_5 <- list()

  histo_D <- function(dat, variable, line, title) {
    ggplot(dat, aes(x = variable, fill = variable)) +
      geom_histogram(fill = "red", alpha = 0.2, bins = 30) +
      geom_vline(xintercept = line, color = "black", linetype = "dashed") +
      ggtitle(title) +
      labs(x = NULL, y = NULL) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0, margin = margin(b = 10)))
  }

  histo_I <- function(dat, variable, line, title) {
    ggplot(dat, aes(x = variable, fill = variable)) +
      geom_histogram(fill = "blue", alpha = 0.2, bins = 30) +
      geom_vline(xintercept = line, color = "black", linetype = "dashed") +
      ggtitle(title) +
      labs(x = NULL, y = NULL) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0, margin = margin(b = 10)))
  }

  var <- list.files(var_dir, pattern = ".tif", full.names = TRUE)
  env <- terra::rast(c(var))
  crs(env) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

  for(i in 1:length(sp)){
    sp_comp <- s_f[-i]
    Ms <- list.files(M_dir, pattern = ".shp")
    Ms_1 <- terra::vect(paste0(M_dir, "/", Ms[[i]]))
    Ms_2f <- terra::crop(env, Ms_1)
    Ms_2f <- terra::mask(Ms_2f, Ms_1)
    Ms_3 <- terra::as.points(Ms_2f)
    if(nrow(Ms_3) < n_back){
      Ms_f <- Ms_3
    } else {
      Ms_f <- Ms_3[sample(nrow(Ms_3), n_back), ]
    }
    Ms_c <- Ms[-i]
    val_1[i] <- NA
    val_2[i] <- NA
    val_3[i] <- NA
    val_4[i] <- NA
    val_5[i] <- NA

    for (j in 1:length(sp_comp)){
      Ms_1 <- terra::vect(paste0(M_dir, "/", Ms_c[[j]]))
      Ms_2c <- terra::crop(env, Ms_1)
      Ms_2c <- terra::mask(Ms_2c, Ms_1)
      Ms_3 <- terra::as.points(Ms_2c)
      if(nrow(Ms_3) < n_back) {
        Ms_comp <- Ms_3
      } else {
        Ms_comp <- Ms_3[sample(nrow(Ms_3), n_back), ]
      }
      setwd(spp_dir)
      sp_f <- read.csv(s_f[[i]])
      sp_f <- terra::vect(sp_f, geom=c("Long", "Lat"), crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", keepgeom=FALSE)
      sp_c <- read.csv(sp_comp[j])
      sp_c <- terra::vect(sp_c, geom=c("Long", "Lat"), crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", keepgeom=FALSE)
      sp_fr <- terra::extract(Ms_2f, sp_f)
      sp_fr <- cbind(sp_f, sp_fr)
      sp_fr <- terra::na.omit(sp_fr)
      sp_cr <- terra::extract(Ms_2c, sp_c)
      sp_cr <- cbind(sp_c, sp_cr)
      sp_cr <- terra::na.omit(sp_cr)
      focal_sp <- enmtools.species()
      comp_sp <- enmtools.species()
      name_f <- substr(s_f[[i]], 1, nchar(s_f[[i]])-4)
      name_c <- substr(sp_comp[j], 1, nchar(sp_comp[j])-4)
      focal_sp$species.name <- name_f
      comp_sp$species.name <- name_c
      focal_sp$presence.points <- sp_fr[, 3:4]
      comp_sp$presence.points <- sp_cr[, 3:4]
      focal_sp$background.points <- Ms_f[, 1:2]
      comp_sp$background.points <- Ms_comp[, 1:2]
      focal_m <- enmtools.ecospat.bg(focal_sp, comp_sp, env[[1:2]], nreps = 1000, bg.source = "points")
      val_1[[i]][j] <- focal_m[[15]]$obs$D
      val_2[[i]][j] <- focal_m[[15]]$p.D
      val_3[[i]][j] <- focal_m[[15]]$obs$I
      val_4[[i]][j] <- focal_m[[15]]$p.I
      q <- as.data.frame(focal_m[[15]]$sim$D)
      colnames(q)[1] <- "simulate"
      p_d <- as.numeric(focal_m[[15]]$obs$D)
      pD <- round(focal_m[[15]]$p.D, 4)
      t <- paste0(name_f, " vs ", name_c, " p = ", pD)
      a <- histo_D(dat = q, variable = q$simulate, line = p_d, title = t)
      w <- as.data.frame(focal_m[[15]]$sim$I)
      colnames(w)[1] <- "simulate"
      p_i <- as.numeric(focal_m[[15]]$obs$I)
      pI <- round(focal_m[[15]]$p.I, 4)
      t <- paste0(name_f, " vs ", name_c, " p = ", pI)
      b <- histo_I(dat = w, variable = w$simulate, line = p_i, title = t)
      qw <- a / b
      p1 <- qw + plot_annotation(theme = theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold")),
                                 tag_levels = 'A', tag_suffix = ')')
      setwd(dir_save)
      ggsave(file = paste0(name_f, "_", name_c, "_600_dpi.jpeg"), plot = p1, width = 20,
             height = 20, dpi = 600, units = "cm", device = "jpeg")
      val_5[[i]][j] <- paste0(name_f, " vs ", name_c)
    }
  }
  D <- as.data.frame(Reduce(function(x, y) c(x, unlist(y)), val_1))
  colnames(D)[1] <- "D"
  p_val <- as.data.frame(Reduce(function(x, y) c(x, unlist(y)), val_2))
  colnames(p_val)[1] <- "p-val_D"
  I <- as.data.frame(Reduce(function(x, y) c(x, unlist(y)), val_3))
  colnames(I)[1] <- "I"
  p_val_I <- as.data.frame(Reduce(function(x, y) c(x, unlist(y)), val_4))
  colnames(p_val_I)[1] <- "p-val_I"
  compa <- as.data.frame(Reduce(function(x, y) c(x, unlist(y)), val_5))
  colnames(compa)[1] <- "Test"
  final <- cbind(compa, D, p_val, I, p_val_I)
  write.csv(final, paste0(dir_save, "/DI_val_", "sim", ".csv"), row.names = FALSE)
}

# ------------------------------------------------------------------------------------------------ #

# Arguments #

# ------------------------------------------------------------------------------------------------ #

#' @param var_dir Directory when the variables in tiff format are contained. A PCA before the analysis is required due to the loop takes the first two variables in the directory

#' @param spp_dir Directory when the csv occurrences are contained

#' @param dir_save Directory when the results are going to be deposited

#' @param M_dir Directory when the shp with the M hypothesis are contained (see Soberon & Peterson, 2005)

#' @param n_back Number of background points. Defauld = 10000. If the analysis area has less than 10000 pixels to sampling, the maximum number of pixels will be the n_back.

#' @param spp A vector containing all the species names as you wanted to be printed in the histogram plots and in the final csv's.

#' @return A csv file containing 5 columns: 1) Test: the species that have tested, 2) D: D values for the overlap, 3) p-val_D: p value of the test, 4) I = I values for the overlap, 3) p-val_I: p value of the test.

#' @examples
#'
#' var_dir <- "D/example/PCA"
#' spp_dir <- "D/example/csv"
#' dir_save <- "D/example/results"
#' M_dir <- "D/example/Ms"
#' n_back <- 10000
#' spp <- c("spp1", "spp2", "spp3", "sppn")
#'
#' sim_spp_test(var_dir = var_dir, spp_dir = spp_dir, dir_save = dir_save, M_dir = M_dir, n_back = n_back, spp = spp)

# ------------------------------------------------------------------------------------------------ #

### EndNotRun
