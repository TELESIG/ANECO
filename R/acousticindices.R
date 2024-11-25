#' Set of acoustic indices
#'
#'The function creates a loop to obtain a series of acoustic indices from a set of audio files. The function calibrates amplitude measurements to obtain indices calculated with absolute values of sound pressure level. Additionally, it has the option to calculate a set of commonly implemented indices using the seewave (Sueur et al., 2016) and soundecology (Villanueva-Rivera & Pijanowski, 2016) packages with default values (without the calibration process), as well as indices useful for detecting audio segments with the presence of noises such as heavy rain and cicadas.
#' @param dir Character. Path of the folder where the audios files are located
#' @param calibparam Calibration parameters. Vector with 3 numeric elements indicating the calibration parameters: micsens, gain, vADC, in that order. Where \strong{micsens:} The sensitivity of the acoustic transducer; \strong{gain:} The user-defined gain settings for the deplyment in decibels; \strong{vADC:} Zero-to-peak voltage of the analogue-to-digital converter. \emph{Note: vADC = V rms * sqrt(2)}. For more information about these parameters see Merchant et al 2018.  Additionally, Numeric values can be replaced by the word \code{'SM4'} to set default values for wildlife acoustics SM4 recorders = c(-35, 48, 1), as well as by the word \code{'SMmini'} and  \code{'SM2+'} to use the default values of SMmini recorders \code{c(-11, 18, 1.5)} and SM2+ \code{c(-36, 48, 1.41)} of the same brand.
#' @param sel.ind Vector containing the index classes to be calculated. Default is \emph{'all'}. See section: \strong{Details}
#' @param exclude Vector containing the index classes that you want to exclude from the calculation. Default is \emph{'none'}. See section: \strong{Details}
#' @param noise.ind Logical. If \code{noise.ind = TRUE} the necessary indices for the prediction of rain and cicadas will be calculated, see \emph{"noise"} in the table of the section \strong{Details}. Default is \code{noise.ind = TRUE}
#' @param bioband Vector with two values, the first value corresponds to the lower limit of the biological frequency band, the second value corresponds to the upper limit, in Hz. Default is \code{bioband = c(2000, 11000)}
#' @param noiseband Vector with two values, the first value corresponds to the lower limit of the anthropic frequency range, the second value corresponds to the upper limit, in Hz. Default is \code{noiseband = c(0, 2000)}
#' @param prefix.format Character. A code indicating the format in which the prefix is set where S represents site letters and G represents letters for the recorder. For example, for files named BAJA08_20210523_050000.wav The prefix is composed of site (BAJA) and recorder (08). Therefore, the prefix format is: \code{"SSSSGG"}. That is, four letters for the site ("S") and two letters for the recorder "G". The default value is \code{"SSGG"}.
#' @param ncores Numeric. Number of cores used in parallel analysis. By default it uses the "total physical cores" minus one. The number of available cores can be checked with the function \code{parallel::detectCores(all.tests = TRUE, logical = FALSE)}
#' @param wl Numeric. Window lenght (in samples)  used in fourier transformation. Default is \code{wl = 512}
#' @param save.file Logic. Do you want to export the results to a csv file automatically?. Default is  \code{save.file = TRUE}
#' @param add.prefix Add a prefix before the date when there is none. The prefix must contain letters for the site and letters for the recorder. By default, if the file name does not contain a prefix, XXXX is added.
#' @details The argument \code{sel.ind} specifies which indexes you want to calculate. The indexes are grouped into the following categories: \code{sel.ind = c("all", "spl", "ent", "bioant", "acousticevents", "tempvar", "efb", "dfb", "trad", "specdist", "clust")}
#'
#' Where: \tabular{lll}{ \bold{sel.ind} \tab \bold{Type} \tab \bold{Indexes that it calculates} \cr "all" \tab Computes all the indices \tab - \cr "spl" \tab Indices based on summary sound intensity \tab \describe{
#' \item{Median background level (mdBGL)}{}
#' \item{Mean sound level (msldB_low, msldB_bio)}{}
#' \item{A-weighted mean sound level  (msldBA_low, msldBA_bio)}{}
#' \item{Average signal amplitude (avgAMP)}{}
#' \item{L10 exceedance level (L10AMP)}{}
#' \item{Median sound level (Mamp)}{}
#' \item{Difference in exceedance levels (dif_L10L90)}{}} \cr "aci" \tab Acoustic complexity index \tab Acoustic complexity index (ACIout) \cr "ent" \tab Entropy indices \tab \describe{
#' \item{Spectral entropy (Hf)}{}
#' \item{Temporal entropy (Ht)}{}
#' \item{Total entropy (EI)}{}
#' \item{Entropy of spectral maxima (Hm)}{}
#' \item{Entropy of spectral variance (Sound Pressure level) (HvPres)}{}
#' \item{Entropy of spectral variance (dB) (HvSPL)}{}} \cr "bioant" \tab Biophony/anthrophony \tab \describe{
#' \item{Normalized difference soundscape index (NDSI)}{}
#' \item{Biophony (BioPh)}{}
#' \item{Anthrophony (AthPh)}{}
#' \item{Ratio of biophony to anthrophony (Bio_anth)}{}} \cr "acousticevents" \tab Difference with background noise \tab \describe{
#' \item{Acoustic activity (AA, AAanth)}{}
#' \item{Count of acoustic events (AAc, AAcanth)}{}
#' \item{Duration of acoustic events (AAdur, AAduranth)}{}} \cr "tempvar" \tab Temporal variation of acoustic energy \tab \describe{
#' \item{Roughness (Rough)}{}
#' \item{Acoustic complexity index (ACI)}{}} \cr "specdist"  \tab Energy distribution in the frequency spectrum  \tab \describe{
#' \item{Acoustic diversity index (ADI)}{}
#' \item{Acoustic evenness index (AEI)}{}
#' \item{Peak frequency (pk)}{}
#' \item{Spectral Kurtosis (pkd)}{}
#' \item{Spectral Skewness (pks)}{}} \cr "clust" \tab Indexes based on cluster analysis  \tab \describe{
#' \item{Spectral persistence (SP2)}{}
#' \item{Spectral diversity (NumCL)}{}} \cr "trad"  \tab Traditional indices (commonly used) \tab \describe{
#' \item{Acoustic diversity index (ADI_raw)}{}
#' \item{Acoustic evenness index (AEI_raw)}{}
#' \item{Acoustic complexity index (ACI_raw)}{}
#' \item{Bioacoustic Index (BIO_raw)}{}
#' \item{Normalized Difference Soundscape Index (NDSI_raw)}{}
#' \item{Median amplitude envelope (MAE_raw)}{}
#' \item{Spectral entropy (Hf_raw)}{}
#' \item{Temporal entropy (Ht_raw)}{}
#' \item{Total entropy (TE_raw)}{}
#' \item{Roughness (rough_raw)}{}
#' \item{Number of peaks (NP_raw)}{}} \cr "spec"     \tab Statistical spectral indices \tab \describe{
#' \item{Centroid of the spectrum (SpecCent)}{}
#' \item{Interquantile of the spectrum (SpecIQR)}{}
#' \item{Third quantile of the spectrum (SpecQ3)}{}
#' \item{Skewness of the spectrum (SpecSkew)}{}
#' \item{Kurtosis of the spectrum (SpecKurt)}{}} \cr "efb" \tab Energy in frequency bands \tab \describe{
#' \item{Average dominant frequency (dfreq)}{}
#' \item{Energy in frequency bands (1 kHz band resolution)}{}} \cr "dfb" \tab Dominance in frequency bands \tab \describe{
#' \item{Average dominant frequency (dfreq)}{}
#' \item{Dominance in frequency bands (1 kHz band resolution)}{}} \cr "noise" \tab Indexes for the classification of rain and cicadas sounds \tab \describe{
#' \item{Signal to noise cicada (S2N_chi)}{}
#' \item{It also includes the indexes included in spec and efb}{}} }
#' @return Generates a dataset with calculated index values defined in \code{ind.sel} for 1 min intervals in each audio file. The resulting dataset also includes an identification of each minute with the columns: Name, Code, Site, Recorder, Date, Year, Month, Day, Hour and Minute.

#' @references \itemize{
#' \item{Buxton, R., McKenna, M. F., Clapp, M., Meyer, E., Stabenau, E., Angeloni, L. M., ... & Wittemyer, G. 2018. Efficacy of extracting indices from large‐scale acoustic recordings to monitor biodiversity. Conservation Biology.}{}
#'
#'\item{Merchant, N. D., Fristrup, K. M., Johnson, M. P., Tyack, P. L., Witt, M. J., Blondel, P., & Parks, S. E. (2015). Measuring acoustic habitats. Methods in Ecology and Evolution, 6(3), 257-265.}{}
#'
#' \item{Sueur, J., Simonis, C., Brown, E., Depraetere, M., Desjonqueres, C., Fabianek, F., Gasc, A., LaZerte, S., Lees, J., Marchal, J., Pavoine, S., Stotz, A., Villanueva-Rivera, L., Ross, Z., Witthoft, C., & Zhivomirov, H. (2016). Paquete para R: Seewave, Sound Analysis and Synthesis. Version 2.0.5}{}
#'
#' \item{Villanueva-Rivera, L. J. & Pijanowski, B. C. 2016. soundecology: Soundscape Ecology. R package version 1.3.2.}{} }
#' @import seewave
#'
#'
#' @export
#'
#' @examples
#'\code{
#'# Load packages
#'
#'library(seewave)
#'
#'# Define the directory of the folder where the audios are located
#'
#'path_folder_audios <- "C:/path/to/your/audio/files"
#'
#'# Set calibration parameters (default for SM4)
#' # c(micsens,gain,vADC)
#' calibration_parameters <- c(-35, 48, 1)
#'
#'# Run the function
#'
#'result <- acoustic_indices(
#'  dir = path_folder_audios,
#'  calibparam = calibration_parameters,
#'  sel.ind = c("all")
#')}

acoustic_indices <- function(dir = NULL, calibparam = "SM4", sel.ind = "all", exclude = "none", noise.ind = TRUE, bioband = c(1000, 11000), noiseband = c(0, 1000), Channel = "both", prefix.format = "SSGG", wl = 512, ncores = NULL, save.file = TRUE, add.prefix = "XXXX"){

#----------------------------------
# Checking requierements
# ---------------------------------

  if(is.null(dir)){
    stop("dir argument not defined. You must specify the path where the audio data is located")
  }

  if(length(calibparam) == 3) {
    micsens = calibparam[1]
    gain = calibparam[2]
    vADC = calibparam[3]
  }
  if(length(calibparam) == 1) {
    if(calibparam == "SM4") {
      micsens = -35
      gain = 48
      vADC = 1
    }
    if(calibparam == "SMmini") {
      micsens = -11
      gain = 18
      vADC = 1.5
    }
    if(calibparam == "SM2+") {
      micsens = -36
      gain = 48
      vADC = 1.41
    }
    if(!calibparam %in% c("SM4", "SMmini", "SM2+")){
      stop("The type of calibration parameters specified is incorrect.")
    }
  }

  if(length(calibparam) != 3 & length(calibparam) != 1 ) {
    stop("Calibration parameters misconfigured, check the function's help page")
  }

if(anyNA(c(micsens, gain, vADC))){
  stop("Some of the calibration parameters have not been defined. Ensure that you have correctly specified the arguments: micsens, gain, and vADC.")
}

  if(!Channel %in% c("left", "right", "both")){
    stop("Incorrect channel definition, it must be one of the following options: 'left', 'right' or 'both' ")
  }

opciones.sel.ind = c("spl", "ent", "bioant", "acousticevents", "tempvar", "clust", "specdist", "efb", "dfb", "trad", "none")

if(sel.ind == "all"){
  sel.ind <- opciones.sel.ind
  if(exclude != "none"){
    sel.ind <- sel.ind[-c(which(sel.ind %in% exclude))]
  }
}



if(!all(sel.ind %in% opciones.sel.ind)){stop("Error in the sel.ind parameter: the indices to be calculated have not been specified or are incorrectly written. Check the documentation using: ?acoustic_indices")}

# ---------------------------------
# Data management
# ---------------------------------

  setwd(dir)
  time1 = proc.time()
  df <- data.frame()

  # Load file names
  files <- list.files(path = dir, pattern = "wav$", ignore.case = T, full.names = F)


  # Working on multiple cores
  noCores <- parallel::detectCores(all.tests = TRUE, logical = FALSE)
  noCores <- noCores - 1
  if(noCores < 1) {
    noCores = 1
  }

  if(!is.null(ncores)){noCores = ncores}
  cl <- parallel::makeCluster(noCores)
  parallel::clusterExport(cl, c("df", "files", "sel.ind", "Channel", "micsens", "gain", "vADC", "noise.ind", "bioband", "noiseband", "prefix.format", "wl", "dir", "save.file", "add.prefix"), envir = environment())
  parallel::clusterEvalQ(cl, {
    library(seewave)
  })
 pbapply::pboptions(type = "txt")
  ap <- pbapply::pblapply(cl = cl, X = files, FUN = function(x){

    name <- as.character(x)
    nameFull <- paste(dir, name, sep = "/")
    readWaveTry <- function(x){tryCatch(tuneR::readWave(x), error=function(e){NA})}
    wavFile <- readWaveTry(name)
    if(!is.na(wavFile)){

      durS <- seewave::duration(wavFile)
      dur <- floor(durS/60)

      minutos <- seq(from = 0, to = dur, by = 1)

      for (i in minutos) {

        if(i == max(minutos)){
          absdur <- durS - (i*60)
          if(absdur < 30){next}
          wav <- tuneR::readWave(nameFull, from = (i*60), to = durS, units = "seconds")

        }else{
          wav <- tuneR::readWave(nameFull, from = (i*60), to = (i*60) + 60, units = "seconds")
          absdur <- 60
        }


        # Mantener asi provisionar hasta solucionar
        canal <- ifelse(Channel == "both", "left", Channel)

        wav <- tuneR::mono(wav, which = canal)

        # Create temporary file
        nombre <- paste0(name,i,".wav")
        tempWav <- paste(dir,nombre, sep = "/")
        tuneR::writeWave(wav, filename = tempWav)

        # Save temporary file name
        ifile <- paste0(strsplit(basename(tempWav), ".wav")[[1]][1], ".wav")

        ## Name format
        filextDay = "%Y%m%d_%H%M%S.wav"
        site = unlist(strsplit( ifile, '_') )[1]
        filename = paste(site, filextDay, sep = "_")

        # Calibration process
        datCalib <- ANECO::PAMCalib(atype = 'PSD', N = wl, timestring = filename, r = 0, outwrite = 1, plottype = "None", calib = 1, envi = "Air", ctype = "TS", Mh = micsens, G = gain, vADC = vADC, tempWav = tempWav, ifile = ifile, welch = "")


        file.remove(tempWav)

        #Compute dBA
        f_dBA <- datCalib[1, 2:ncol(datCalib)]
        # Decibel reference
        aweight <- vector()
        for ( w in 1:length(f_dBA)) {
          aweight[w] <- seewave::dBweight(f_dBA[w])$A
        }
        a <- datCalib[2:nrow(datCalib),2:ncol(datCalib)]
        aA = t(t(a) + aweight)


        # convert to pressure
        press <- rowSums(10^(aA/10))
        dBA = 10*log10(press) #hist(dBA)

        #Data manipulation
        FullMatrixA = as.data.frame(aA)
        rm(aA)
        FullMatrix <- as.data.frame(datCalib)[-1, -1]
        colnames(FullMatrix) <- round(datCalib[1, 2:ncol(datCalib)])
        colnames(FullMatrixA) <- round(datCalib[1, 2:ncol(datCalib)])
        rm(datCalib)
        vFreq <- as.numeric(colnames(FullMatrix))

        # Select "Biological" frequency bands
        PosBio <- which(vFreq >= bioband[1] & vFreq <= bioband[2])
        PosBio_standard <- which(vFreq >= 2000 & vFreq <= 8000)
        BioMatrix <- FullMatrix[, PosBio]
        BioMatrixA <- FullMatrixA[, PosBio]
        BioMatrix_standard <- FullMatrix[, PosBio_standard]


        # Select anthropic frequency bands
        PosNoise <- which(vFreq >= noiseband[1] & vFreq < noiseband[2])
        PosNoise_standard <- which(vFreq >= 0 & vFreq < 2000)
        AntMatrix <- FullMatrix[, PosNoise]
        AntMatrixA <- FullMatrixA[, PosNoise]
        AntMatrix_standard <- FullMatrix[, PosNoise_standard]


        # Number of frequency bands
        nFQ = as.numeric( dim(BioMatrix)[2])
        nFQbk =  as.numeric( dim(AntMatrix)[2])
        # Number of time frames
        fileDur = as.numeric( dim(BioMatrix)[1])

        # Matrix of sound pressure values
        BioMatrixPress <- 10^(BioMatrix/10)
        BioMatrixPress_standard <- 10^(BioMatrix_standard/10)
        AntMatrixPress <- 10^(AntMatrix/10)
        AntMatrixPress_standard <- 10^(AntMatrix_standard/10)

        normdBA <- function (x) { (x-(-10)) / (80 - (-10)) }


        # Function to calculate background noise
        # Modified here as a level for each frequency band

        # recordar que se puede quitar lo de los minimos
        BGL <- function(datos) {
          output <- NULL
          for (ff in 1:dim(datos)[2]) { #loop through each frequency band
            #(1) compute a histogram
            pretmp <- sort(datos[,ff]) #ordenar
            tmp <- pretmp[pretmp <= quantile(pretmp, 0.99)] # excludes the top 1% of the values to avoid outliers
            if(min(tmp) == -Inf){
              output <- 0
              break()
            }
            brks = seq(from = min(tmp), to = max(tmp), length.out = length(tmp)/8)
            histo <- hist(tmp, breaks = brks, plot = FALSE)
            tab <- data.frame(counts = histo[["counts"]], mids = histo[["mids"]])
            tmpval <- tab[,1]
            tmpnam <- tab[,2]
            colmax <- as.numeric(which.max(tmpval))

            #(4) accumulating counts in histogram bins below (3) until 68% of counts below (3) is reached
            cutoff = sum( tmpval[ 1:colmax-1] ) * .68 #value to stop the summation at
            # got an error when the first bin had the max # of values..... so just  use the min value
            if (cutoff == 0) {stploc = 0}
            cntbk <- 0
            for (h in 1:length(tmpval[1:colmax])-1 )
            {
              if (cntbk > cutoff) {
                break }
              loc = colmax - h #find the location to start the summation
              cntbk = cntbk + tmpval[loc]
              stploc = tmpnam[loc-1]
            }
            # (5) final calc: (3) + N(4)
            if (cutoff == 0) {output[ff] = tmpnam[colmax] } else { output[ff] = tmpnam[colmax] + stploc*0.1 }
          }
          return(output)
        }

        # Background noise for each frequency band
        BK_Towsey <- BGL(datos = BioMatrix)
        BK_Towsey_anth <- BGL(datos = AntMatrix)

        if (sum(BK_Towsey) == 0 | sum(BK_Towsey_anth) == 0){
          next()
        }

# --------------------------------------------
# Median background level
# --------------------------------------------
 median_BGL <- NA
        if(any(c("spl") %in% sel.ind) | noise.ind) {
          median_BGL <- median(BK_Towsey)
        }

# --------------------------------------------
# Acoustic complexity index (ACI) calibrado
# --------------------------------------------
 ACIout = NA

        if(any(c("tempvar") %in% sel.ind)){
        t <- array(0, dim(BioMatrix) - c(1, 0))
        for (frq in 1:(dim(t)[2])){
          t[, frq] <- abs(diff(BioMatrixPress[, frq])) / sum(BioMatrixPress[, frq])
        }
        ACIout = sum(t)
        rm(t)
        }

# --------------------------------------------
# Mean sound level
#----------------------------------------------
   msldBA_low = NA
   msldBA_bio = NA
   msldB_low = NA
   msldB_bio = NA

     if(any(c("spl") %in% sel.ind) | noise.ind){
       # dBA
       msldBA_low = 10*log10( sum( 10^(AntMatrixA/10))/ nrow(AntMatrixA))
       msldBA_bio = 10*log10( sum( 10^(BioMatrixA/10))/ nrow(BioMatrixA) )
     }
     if(any(c("spl") %in% sel.ind)){
       # dB
       msldB_low  = 10*log10( sum( AntMatrixPress)/ nrow(AntMatrix))
       msldB_bio = 10*log10( sum(BioMatrixPress)/ nrow(BioMatrix) )
     }

# --------------------------------------------
#  Average signal amplitude
# --------------------------------------------
   avgAMP = NA

#Towsey et al., 2013
# modified as the mean dBA value for the timestep, normalized by the min/max (set at: -10 to 80 dB)
# modified to fit data

        if(any(c("spl") %in% sel.ind)){
        avgAMP = normdBA(mean(dBA))
        }

# --------------------------------------------
# L10 exceedance level
# --------------------------------------------
        L10AMP = NA
        if(any(c("spl") %in% sel.ind)){
        L10AMP = normdBA(quantile (dBA, .1) ) #10 porcentil
        }

# --------------------------------------------
# Entropy indices, code developed by Buxton et al. 2018
# --------------------------------------------

 Hf = NA
 Ht = NA
 EI = NA

        if(any(c("ent") %in% sel.ind)){
        Pres = colMeans(BioMatrixPress)  #Leq for each Fq band over time period (mean of the pressures)
        Pres2 = Pres/sum(Pres)
        Hf = -sum( (Pres2 * log(Pres2))) /log(nFQ)
        Leqt = rowMeans( BioMatrixPress) #Leq for each second over entire band (could also used dBA)
        Leqt2 = Leqt/sum(Leqt)
        Ht = -sum( (Leqt2 * log(Leqt2)) / log(fileDur) )
        EI = Ht * Hf
        }

 # --------------------------------------------
 # Relationship between biophonies and anthropophonies
 # --------------------------------------------
 NDSI = NA
 BioPh = NA
 AthPh = NA
 Bio_anth = NA

      if(any(c("bioant") %in% sel.ind)){

        BioPh = sum( colMeans(BioMatrixPress_standard) )
        AthPh = sum( colMeans(AntMatrixPress_standard) )
        NDSI = (BioPh -AthPh ) / (BioPh + AthPh )
        Bio_anth = (BioPh /AthPh )
      }

 # --------------------------------------------
 # Difference with background noise
 # --------------------------------------------

 AAanth = NA
 AA = NA
 AAcanth = NA
 AAc = NA
 AAduranth = NA
 AAdur = NA

    if(any(c("acousticevents") %in% sel.ind)){
      AA = NULL
      AAc = NULL
      s2n <- NULL
      for (f in 1:length(BioMatrix) ) {
        AA[f] =  length( BioMatrix[, f][ BioMatrix[, f] > BK_Towsey[f] + 3 ] )/ length (BioMatrix[, f]) #proporcion
        AAc[f] = (length( BioMatrix[, f][ BioMatrix[, f] > BK_Towsey[f] + 3 ] ))
        s2n[f] = max((BioMatrix[,f])) - BK_Towsey[f]
      }
      AA= sum(AA)
      AAc = sum(AAc)
      AAdur = NULL
      for (f in 1:length(BioMatrix) ) {
        #logical matrix...
        temp <- BioMatrix[,f] > BK_Towsey[f] + 3

        temp2 <- temp*1 # convertir a 0/1
        #encontrar numeros consecutivos
        rl <- rle(temp2)
        len = rl$lengths
        v =  rl$value
        if (length(v) == 1 )
        {
          if   (v == 0) {AAdur[f] = 0 }
          else {AAdur[f] = len }
          next
        }
        cumsum = NULL
        cntAA = 0
        for ( qq in  seq(from = 2, to = length(v), by =2) )
        {
          cntAA = cntAA + 1
          cumsum[cntAA] = len[qq]
        }
        AAdur[f] = mean(cumsum)
      }
      AAdur = median(AAdur)


      AAanth = NULL
      AAcanth = NULL
      AAduranth = NULL

      for (f in 1:length(AntMatrix) ) {
        AAanth[f] =  length( AntMatrix[, f][ AntMatrix[, f] > BK_Towsey_anth[f]+3 ] )/ length (AntMatrix[, f])
        AAcanth[f] = length( AntMatrix[, f][ AntMatrix[, f] > BK_Towsey_anth[f]+3 ] )
      }

      AAanth= sum(AAanth)
      AAcanth = sum(AAcanth)

      for (f in 1:length(AntMatrix) ) {

        temp <- AntMatrix[, f] > BK_Towsey_anth[f] + 3

        temp2 <- temp*1 # convertir a 0-1
        #encontrar valores consecutivos
        rl <- rle(temp2)
        len = rl$lengths
        v =  rl$value
        if (length(v) == 1 )
        {
          if   (v == 0) {AAduranth[f] = 0 }
          else {AAduranth[f] = len }
          next
        }
        cumsum = NULL
        cntAA = 0
        for ( qq in  seq(from = 2, to = length(v), by = 2) )
        {
          cntAA = cntAA +1
          cumsum[cntAA] = len[qq]
        }
        AAduranth[f] = mean(cumsum)
      }

      AAduranth = median(AAduranth)
    }

 # --------------------------------------------
 # Roughness
 # - ind.Rough = TRUE
 # --------------------------------------------

 Rough = NA

  if(any(c("rough") %in% sel.ind)){
    Rough = NULL
    for (f in 1:length(BioMatrixPress) ){
      x = BioMatrixPress[, f]
      x <- x/max(x)
      deriv2 <- diff(x, 1, 2)
      Rough[f] <- sum(deriv2^2, na.rm = TRUE)
    }
    Rough = median(Rough)
  }


 # --------------------------------------------
 # Acoustic diversity, acoustic eveness
 # - sel.ind = c("adi", "aei")
 # --------------------------------------------
 ADI_step = NA
 ADI_fix = NA
 ADI_1kHz_fix = NA
 AEI_step = NA
 AEI_fix = NA
 AEI_1kHz_fix = NA

 # AEI_step and ADI_step (spectrum divided into 256 frequency bands, background noise threshold adjusted to each band)
 # ADI and AEI Fix ( spectrum divided into 256 freq bands and background noise threshold fixed as the median Background level.
 # ADI and AEI Fix (spectrum divided into 1 kHz freq  bands, and background noise threshold fixed as the median Background level)


if(any(c("specdist") %in% sel.ind)){
  Score = NULL
        for (f in 1:length(BioMatrix) ){
          Score[f] = length( BioMatrix[, f][BioMatrix[, f] > BK_Towsey[f] + 3] )/ length (BioMatrix[, f])
        }

      ADI_step = vegan::diversity(Score, index = "shannon")

      AEI_step = ineq::Gini(Score)

      Score2 = NULL
      for (f in 1:length(BioMatrix) ){
        Score2[f] = length( BioMatrix[, f][BioMatrix[, f] > (median(BK_Towsey) + 3)] )/ length (BioMatrix[, f])
      }

      ADI_fix = vegan::diversity(Score2, index = "shannon")

      AEI_fix = ineq::Gini(Score2)


      rd <- floor(as.numeric(names(BioMatrix))/1000)

      fs <- unique(rd)
      f_perc <- NA
      c1 = 1

      for (frqb in fs) {
      scorefi = NA
      c2 = 1
      fb <- which(rd == frqb)
      for (fi in fb) {
        scorefi[c2] <- length(BioMatrix[, fi][BioMatrix[, fi] > median(BK_Towsey) + 3])
        c2 <- c2 + 1
      }
      f_perc[c1] = sum(scorefi)/(length(fb)*nrow(BioMatrix))

      c1 <- c1 + 1
      }

      f_prop <- f_perc/sum(f_perc)

      ADI_1kHz_fix = vegan::diversity(f_prop, index = "shannon")

      AEI_1kHz_fix = ineq::Gini(f_prop)
}
 # --------------------------------------------
 # Energy distribution in the spectrum
 #  - ind.freqdist = TRUE
 # --------------------------------------------

 pk = NA
 pkd = NA
 pks = NA
 Hm = NA
 HvPres = NA
 HvSPL = NA


if(any(c("ent") %in% sel.ind)){
peakf = NULL
        for (j in 1:dim(BioMatrix)[1] )
        {   peakf[j] = (which.max( BioMatrix[j, ] ) )  }
        pk2 = matrix(0, 1, dim(BioMatrix)[2])
        for (uu in 1:nFQ)  { pk2[uu] = sum(peakf == uu)  }
        colnames(pk2) = colnames(BioMatrix)
        pk2nor = pk2/fileDur
        pk = as.numeric(gsub("H", "", colnames(BioMatrix[which.max(pk2)])))

        # kurtosis
        pkd = moments::kurtosis(as.vector(pk2nor))

        # skewness
        pks = moments::skewness(as.vector(pk2nor))

        # Entropy of Spectral Maxima
        pk2[pk2 == 0] <- 1e-07
        pk_prob = pk2/(sum(pk2)) # normalize
        Hm = -sum((pk_prob * log2(pk_prob))) / log2(nFQ)

        # Entropy of spectral variance- pressure and dB
        Press = (BioMatrixPress)
        Pv = NULL
        for (v in 1:dim(Press)[2] ) {Pv[v] = var(Press[ , v]) }
        Pv2 = Pv/sum(Pv)
        HvPres = -sum( (Pv2 * log2(Pv2)) ) / log2(nFQ)

        Pv = NULL
        for (v in 1:dim(BioMatrix)[2] ) {Pv[v] = var(BioMatrix[,v])  }
        Pv2 = Pv/sum(Pv)
        HvSPL = -sum( (Pv2 * log2(Pv2)) ) / log2(nFQ)

        pk = NA
        pkd = NA
        pks = NA
 }


 if(any(c("specdist") %in% sel.ind)){
   peakf = NULL
   for (j in 1:dim(BioMatrix)[1] )
   {   peakf[j] = (which.max( BioMatrix[j, ] ) )  }
   pk2 = matrix(0, 1, dim(BioMatrix)[2])
   for (uu in 1:nFQ)  { pk2[uu] = sum(peakf == uu)  }
   colnames(pk2) = colnames(BioMatrix)
   pk2nor = pk2/fileDur
   pk = as.numeric(gsub("H", "", colnames(BioMatrix[which.max(pk2)])))

   # kurtosis
   pkd = moments::kurtosis(as.vector(pk2nor))

   # skewness
   pks = moments::skewness(as.vector(pk2nor))
 }



 # --------------------------------------------
 # Normalize exceedence levels for dBA
 # --------------------------------------------
 dif_L10L90 = NA

if(any(c("spl") %in% sel.ind)){
        Exceed_norm    = normdBA ( quantile(dBA, c(0.05, 0.1, 0.5, 0.9, 0.95),na.rm = TRUE) )
        L10n = (Exceed_norm[4]) #L10 = Porcentil 90
        L90n = (Exceed_norm[2]) #L90 = Porcentil 10
        dif_L10L90 = L10n - L90n
}

 # --------------------------------------------
 # Median sound level
 # --------------------------------------------
 Mamp = NA

if(any(c("spl") %in% sel.ind)){
  Mamp = quantile( 10*log10(rowMeans(BioMatrixPress)), 0.5)
}


 # --------------------------------------------
 # Cluster indices
 # --------------------------------------------
SP2 = NA
NumCL = NA
if(any(c("clust") %in% sel.ind)){
        x = scale(BioMatrix)  # To standarize the variables

        res = NbClust::NbClust(x, distance = "euclidean", min.nc=2, max.nc=10, method = "kmeans", index = "silhouette")
        # use all algothrums to find optimal cluster... takes way to long!

        x2 =  as.data.frame(t(res$Best.nc))
        x4 = (table(x2$Number_clusters))
        NumCL =  as.numeric(names(which.max(x4))) # "R"

        #......................................
        # (30) Spectral persistence
        BP = as.numeric(res$Best.partition)
        BPrle <- rle(BP)# need to find how many seconds each cluster "persisted" for
        # average duration of the clusters which persist for longer than one frame, so don't count 1s!
        SP2 = mean(BPrle$lengths[BPrle$lengths>1])
}

# --------------------------------------------
# Indices for rain and cicada detection
# --------------------------------------------

ind_noise = NA
if(noise.ind){
  # Prediccion de ruidos
  fr_chi <- which(vFreq >= 3500 & vFreq <= 5500)
  BioMatrix_chi <- FullMatrix[, fr_chi]
  q75_chi <- quantile(as.matrix(10^(BioMatrix_chi/10)), 0.75)
  BG_amp <- 10^(median(BK_Towsey)/10)
  ind_noise <- (q75_chi - BG_amp)/(q75_chi + BG_amp)
}

# Frequence energy
freq.energy <- NA
if(any(c("efb") %in% sel.ind) | noise.ind){
  freq.energy <- list()
  fqmax <- max(vFreq)
  lband <- 1000
  countpos <- 1
  while (lband <= fqmax) {
    nband <- which(vFreq >= (lband - 1000) & vFreq < lband)
    fband <- FullMatrix[, nband]
    energydB = 10 * log10(sum(10^(fband/10))/nrow(fband))
    freq.energy[[countpos]] <- energydB
    names(freq.energy)[countpos] <- paste0("E", countpos - 1, countpos)
    countpos = countpos + 1
    lband = lband + 1000
}
}


#///////////////////////////////////////////////////////////
# The following are indices that do not go through the calibration process
#//////////////////////////////////////////////////////////

# Extract the corresponding audio minute

if(i == max(minutos)){
  absdur <- durS - (i*60)
  if(absdur < 30){next}
  minwav <- tuneR::readWave(name, from = (i*60), to = durS, units = "seconds")
}else{
  minwav <- tuneR::readWave(name, from = (i*60), to = (i*60) + 60, units = "seconds")
}


minwav <- tuneR::channel(minwav, which = Channel)

sampling_rate <- minwav@samp.rate



# Calculate spectrum values
     specm <- seewave::meanspec(minwav, plot = FALSE)

     if(Channel == "both") {
       specm2 <- seewave::meanspec(minwav, channel = 2, plot = FALSE)
       specm <- data.frame(x = specm[ ,1], y = (specm[ ,2] + specm2[ ,2])/2)
       specm <- as.matrix(specm)
     }


 # --------------------------------------------
 # Dominant frequencies
 # --------------------------------------------

Prom.Dom.frec = NA
length.fd.01 <- NA
length.fd.12 <- NA
length.fd.23 <- NA
length.fd.23 <- NA
length.fd.34 <- NA
length.fd.45 <- NA
length.fd.56 <- NA
length.fd.67 <- NA
length.fd.78 <- NA
length.fd.89 <- NA
length.fd.910 <- NA
length.fd.910 <- NA
length.fd.1011 <- NA
length.fd.11 <- NA


fmax_domfreq <- ifelse(sampling_rate/2 >= 11000, 11000, sampling_rate/2)
fmin_domfreq <- 300

dom.freq <- as.data.frame(seewave::dfreq(minwav, channel = 1, f = minwav@samp.rate, wl = 512, threshold = 0.015, plot = FALSE, bandpass = c(fmin_domfreq, fmax_domfreq)))

if(Channel == "both") {
  dom.freq2 <- as.data.frame(seewave::dfreq(minwav, channel = 2, f = minwav@samp.rate, wl = 512, threshold = 0.015, plot = FALSE, bandpass = c(fmin_domfreq, fmax_domfreq)))
  dom.freq <- data.frame(x = dom.freq$x, y = (dom.freq$y + dom.freq2$y)/2)
}

if(any(c("dfb") %in% sel.ind) | noise.ind){

  Prom.Dom.frec <- mean(dom.freq$y, na.rm = TRUE)
  dm.01 <- dom.freq[dom.freq$y <= 1, ]
  dm.12 <- dom.freq[dom.freq$y > 1 & dom.freq$y <= 2, ]
  dm.23 <- dom.freq[dom.freq$y > 2 & dom.freq$y <= 3, ]
  dm.34 <- dom.freq[dom.freq$y > 3 & dom.freq$y <= 4, ]
  dm.45 <- dom.freq[dom.freq$y > 4 & dom.freq$y <= 5, ]
  dm.56 <- dom.freq[dom.freq$y > 5 & dom.freq$y <= 6, ]
  dm.67 <- dom.freq[dom.freq$y > 6 & dom.freq$y <= 7, ]
  dm.78 <- dom.freq[dom.freq$y > 7 & dom.freq$y <= 8, ]
  dm.89 <- dom.freq[dom.freq$y > 8 & dom.freq$y <= 9, ]
  dm.910 <- dom.freq[dom.freq$y > 9 & dom.freq$y <= 10, ]
  dm.1011 <- dom.freq[dom.freq$y > 10 & dom.freq$y <= 11, ]
  dm.11 <- dom.freq[dom.freq$y > 11, ]

  length.fd.01 <- length(dm.01$y)
  length.fd.12 <- length(dm.12$y)
  length.fd.23 <- length(dm.23$y)
  length.fd.23 <- length(dm.23$y)
  length.fd.34 <- length(dm.34$y)
  length.fd.45 <- length(dm.45$y)
  length.fd.56 <- length(dm.56$y)
  length.fd.67 <- length(dm.67$y)
  length.fd.78 <- length(dm.78$y)
  length.fd.89 <- length(dm.89$y)
  length.fd.910 <- length(dm.910$y)
  length.fd.1011 <- length(dm.1011$y)
  length.fd.11 <- length(dm.11$y)
  }


if(any(c("trad") %in% sel.ind) | noise.ind){

  Prom.Dom.frec <- mean(dom.freq$y, na.rm = TRUE)
}

# --------------------------------------------
# Traditional acoustic indices
# - sel.ind = "trad"
# --------------------------------------------

AEIprom = NA
ADIprom = NA
ACIprom = NA
ACIseewave = NA
BIOprom = NA
NDSIprom = NA
TEraw = NA
Ht_raw = NA
Hf_raw = NA
Rough_see = NA
MAE = NA
NP = NA


if(any(c("trad") %in% sel.ind)){
  AEIraw <- soundecology::acoustic_evenness(minwav)
  AEIprom <- AEIraw$aei_left
  ADIraw <- soundecology::acoustic_diversity(minwav)
  ADIprom <- ADIraw$adi_left
  ACIraw <- soundecology::acoustic_complexity(minwav)
  ACIprom <- ACIraw$AciTotAll_left
  ACIseewave <- seewave::ACI(minwav)
  BIOraw <- soundecology::bioacoustic_index(minwav)
  BIOprom <- BIOraw$left_area
  NDSIraw <- soundecology::ndsi(minwav)
  NDSIprom <- NDSIraw$ndsi_left
  TEraw <- seewave::H(minwav)
  envorni <- seewave::env(minwav, f = minwav@samp.rate, plot = FALSE)
  Ht_raw <- seewave::th(envorni)
  speca <- seewave::spec(minwav, f = minwav@samp.rate, plot = FALSE)
  Hf_raw <- seewave::sh(speca)
  Rough_see <- seewave::roughness(speca) #diferente por defecto, usa todo el espectro
  MAE <- seewave::M(minwav)
  peaks <- seewave::fpeaks(specm, plot = FALSE)
  NP <- length(peaks)/2

  if(Channel == "both"){
    AEIprom <- (AEIraw$aei_left + AEIraw$aei_right)/2
    ADIprom <- (ADIraw$adi_left + ADIraw$adi_right)/2
    ACIprom <- (ACIraw$AciTotAll_left + ACIraw$AciTotAll_right)/2
    ACIseewave <- (ACIseewave + seewave::ACI(minwav, channel = 2))/2
    BIOprom <- (BIOraw$left_area + BIOraw$right_area)/2
    NDSIprom <- (NDSIraw$ndsi_left + NDSIraw$ndsi_right)/2
    TEraw <- (TEraw + seewave::H(minwav, channel = 2))/2
    envorni2 <- seewave::env(minwav, f = minwav@samp.rate, channel = 2, plot = FALSE)
    Ht_raw <- (Ht_raw + seewave::th(envorni2))/2
    speca2 <- seewave::spec(minwav, f = minwav@samp.rate, channel = 2, plot = FALSE)
    Hf_raw <- (Hf_raw + seewave::sh(speca2))/2
    Rough_see <- (Rough_see + seewave::roughness(speca2))/2 #diferente por defecto, usa todo el espectro
    MAE <- (MAE + seewave::M(minwav, channel = 2))/2
  }

}

# --------------------------------------------
# Spectral indices
# --------------------------------------------

cent = NA
iqr = NA
tq = NA
skew = NA
kurt = NA

if(noise.ind){

  specStats <- seewave::specprop(specm, f = minwav@samp.rate, flim = c(fmin_domfreq/1000, fmax_domfreq/1000), plot = FALSE)
  cent <- specStats$cent
  iqr <- specStats$IQR
  tq <- specStats$Q75
  skew <- specStats$skewness
  kurt <- specStats$kurtosis
}

# -------------------------------------
# Identification data
#--------------------------------------
   namesplit <- unlist(strsplit(name, "_"))
   if(length(namesplit == 3)){
   conwav = namesplit[3]
   hora <- as.character(unlist(strsplit(conwav, ".wav"))[1])
   codigo <- namesplit[1]
   sitioindex <- unlist(gregexpr('S', prefix.format))
   grabadoraindex <- unlist(gregexpr('G', prefix.format))
   sitio <- substr(codigo, start = min(sitioindex), stop = max(sitioindex))
   grabadora <- substr(codigo, start = min(grabadoraindex), stop = max(grabadoraindex))
   fecha <- namesplit[2]
   ano <- substr(fecha, start = 1, stop = 4)
   mes <- substr(fecha, start = 5, stop = 6)
   dia <- substr(fecha, start = 7, stop = 8)
   }

if(length(namesplit) == 2) {
  name2 <- paste(add.prefix, name, sep = "_")
  namesplit <- unlist(strsplit(name2, "_"))
  conwav = namesplit[3]
  hora <- as.character(unlist(strsplit(conwav, ".wav"))[1])
  codigo <- namesplit[1]
  sitioindex <- unlist(gregexpr('S', prefix.format))
  grabadoraindex <- unlist(gregexpr('G', prefix.format))
  sitio <- substr(codigo, start = min(sitioindex), stop = max(sitioindex))
  grabadora <- substr(codigo, start = min(grabadoraindex), stop = max(grabadoraindex))
  fecha <- namesplit[2]
  ano <- substr(fecha, start = 1, stop = 4)
  mes <- substr(fecha, start = 5, stop = 6)
  dia <- substr(fecha, start = 7, stop = 8)
}

    z <- list(Name = paste0(name,i), Code = codigo, Site = sitio, RecorderID = grabadora, Date = fecha, Year = ano, Month = mes, Day = dia, Hour = hora, Minute = i, Duration = absdur, TE_seewave = TEraw,  Ht_seewave = Ht_raw, Hf_seewave = Hf_raw, MAE_seewave = MAE, NP = NP, ACI_soundecology = ACIprom, ACI_seewave = ACIseewave, rough_seewave = Rough_see, ADI_soundecology = ADIprom, AEI_soundecology = AEIprom, BIO_soundecology = BIOprom, NDSI_soundecology = NDSIprom, msldBA_low = msldBA_low, msldBA_bio = msldBA_bio, msldB_low = msldB_low, msldB_bio = msldB_bio, BioPh = BioPh, AthPh = AthPh, avgAMP = avgAMP, L10AMP = L10AMP, AAanth = AAanth, AA = AA, AAcanth = AAcanth, AAc = AAc, AAduranth = AAduranth, AAdur = AAdur, pk = pk, pkd = pkd, pks = pks, Hf = Hf,Ht = Ht, NDSI = NDSI, Bio_anth = Bio_anth,  rough = Rough, EI = EI, Hm = Hm, HvPres = HvPres, HvSPL = HvSPL, ACI = ACIout, ADI_step = ADI_step, ADI_fix = ADI_fix, ADI_1kHz_fix = ADI_1kHz_fix, AEI_step = AEI_step, AEI_fix = AEI_fix, AEI_1kHz_fix = AEI_1kHz_fix, SP2 = SP2, NumCL = NumCL, Mamp = Mamp, dif_L10L90 = dif_L10L90,  mdBGL = median_BGL, SpecCent = cent, SpecIQR = iqr, SpecQ3 = tq, SpecSkew = skew, SpecKurt = kurt, S2N_chi = ind_noise, dfreq = Prom.Dom.frec, freq.energy, F01 = length.fd.01, F12 = length.fd.12, F23 = length.fd.23, F34 = length.fd.34, F45 = length.fd.45, F56 = length.fd.56, F67 = length.fd.67, F78 = length.fd.78, F89 = length.fd.89, F910 = length.fd.910, F1011 = length.fd.1011, F11 = length.fd.11)

        # Creating the result dataframe
        df <- rbind(df, data.frame(z))
    }
    }else{
      df <- data.frame(Name = name, Code = NA, Site = NA, RecorderID = NA, Date = NA, Year = NA, Month = NA, Day = NA, Hour = NA, Minute = NA, Duration = NA, TE_seewave = NA,  Ht_seewave = NA, Hf_seewave = NA, MAE_seewave = NA, NP = NA, ACI_soundecology = NA, ACI_seewave = NA, rough_seewave = NA, ADI_soundecology = NA, AEI_soundecology = NA, BIO_soundecology = NA, NDSI_soundecology = NA, msldBA_low = NA, msldBA_bio = NA, msldB_low = NA, msldB_bio = NA, BioPh = NA, AthPh = NA, avgAMP = NA, L10AMP = NA, AAanth = NA, AA = NA, AAcanth = NA, AAc = NA, AAduranth = NA, AAdur = NA, pk = NA, pkd = NA, pks = NA, Hf = NA,Ht = NA, NDSI = NA, Bio_anth = NA,  rough = NA, EI = NA, Hm = NA, HvPres = NA, HvSPL = NA, ACI = NA, ADI_step = NA, AEI_step = NA, ADI_fix = NA, AEI_fix = NA, ADI_1kHz_fix = NA, AEI_1kHz_fix = NA, SP2 = NA, NumCL = NA, Mamp = NA, dif_L10L90 = NA,  mdBGL = NA, SpecCent = NA, SpecIQR = NA, SpecQ3 = NA, SpecSkew = NA, SpecKurt = NA, S2N_chi = NA, dfreq = NA, freq.energy = NA, F01 = NA, F12 = NA, F23 = NA, F34 = NA, F45 = NA, F56 = NA, F67 = NA, F78 = NA, F89 = NA, F910 = NA, F1011 = NA, F11 = NA)
    }
    return(df)
  }
  )
  parallel::stopCluster(cl)
  nd <- data.frame()
  for (a in 1:length(ap)) {
    nd <- rbind(nd, ap[[a]])
  }
  cat("!Analysis completed! \n")

  ptimeS <- (proc.time()-time1)[3]
  phours <- floor(ptimeS/60/60)
  pmin <- floor(ptimeS/60) - phours*60
  cat('Processing time: ', phours, "h ", pmin, 'min \n')

  rownames(nd) <- NULL
  NAcol <- apply(X = nd, MARGIN = 2, FUN = function(x){
    all(is.na(x[-1]))})

  filtro1 <- nd[ ,c(which(NAcol == FALSE))]

  NArow <- apply(X = filtro1, MARGIN = 1, FUN = function(x){
    all(is.na(x[-1]))})
  filtro2 <- filtro1[c(which(NArow == FALSE)), ]

  errores <- filtro1[c(which(NArow == TRUE)), ][, 1]
  if(length(errores) > 0){
   warning("The following files produced errors: \n", paste(errores, collapse = ", "))
  }

  if(save.file) {
    codigoarchivo <- filtro2$Code[1]
    write.csv(filtro2, paste0(codigoarchivo, "_", format(Sys.time(), "%Y%m%d_%H%M") ,".csv"),row.names = FALSE)
  }
  return(filtro2)
}
