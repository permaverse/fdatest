#' NASA daily temperatures data set
#'
#' It contains the daily mean temperatures registered from July 1983 to June
#' 2005 and stored in the NASA database Earth Surface Meteorology for Solar
#' Energy of two different geographical locations: the region (45-46 North, 9-10
#' East), including the city of Milan (Italy), and the region (48-49 North, 2-3
#' East), including the city of Paris (France).
#'
#' @format List of 2 elements:
#' 
#' - `milan`: Matrix of dimensions `c(22, 365)` containing the daily mean
#' temperatures of the region (45-46 North, 9-10 East), including the city of
#' Milan (Italy) registered from July 1983 to June 2005 (22 years).
#' - `paris`: Matrix of dimensions `c(22, 365)` containing the daily mean
#' temperatures of the region (48-49 North, 2-3 East), including the city of
#' Paris (France) registered from July 1983 to June 2005 (22 years).
#' 
#' @source These data were obtained from the NASA Langley Research Center
#'   Atmospheric Science Data Center Surface meteorological and Solar Energy
#'   (SSE) web portal supported by the NASA LaRC POWER Project. Data are freely
#'   available at: [NASA Surface Meteorology and Solar Energy, A Renewable
#'   Energy Resource web site (release 6.0)](http://eosweb.larc.nasa.gov).
"NASAtemp"
