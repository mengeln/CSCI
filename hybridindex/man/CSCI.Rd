\name{CSCI}
\alias{CSCI}
\title{Score samples using the CSCI tool}
\usage{
  CSCI(bugs, stations, rand = sample.int(10000, 1))
}
\arguments{
  \item{bugs}{A data frame with BMI data (see details)}

  \item{stations}{A data frame with environmental data (see
  details)}

  \item{rand}{An integer to control the RNG seed for the
  subsampling. By default set to \code{sample.int(10000,
  1)}}
}
\value{
  A list of data frames that serve as reports in varying
  detail: \item{core}{A high level summary of the CSCI
  results} \item{Suppl1_mmi}{A detailed breakdown of the
  results of the MMI component of the CSCI}
  \item{Suppl1_grps}{Probability of biotic group membership
  in a SampleID by Group format} \item{Suppl1_OE}{A
  detailed breakdown of the results of the O/E component of
  the CSCI} \item{Suppl2_mmi}{Similar to Suppl1_mmi, except
  broken down by replicates} \item{Suppl2_OE}{Similar to
  Suppl1_OE, except brown down by replicates}
}
\description{
  A function that aggregates all of the steps involved in
  the CSCI scoring process into a single tool.
}
\details{
  A valid "bugs" data frame consists of the following
  columns: StationCode, SampleID, FinalID (i.e., taxa
  names), LifeStageCode ("A", "L", "P", or "X") BAResult
  (i.e., taxa counts), and Distinct ("1" where the
  taxonomist has indicated distinctiveness, else left
  blank).

  A valid "stations" data frame consists of the following
  columns: StationCode (must match with same column in the
  "bugs" data frame), New_Lat, New_Long, ELEV_RANGE,
  BDH_AVE, PPT_00_09, LPREM_mean, KFCT_AVE, TEMP_00_09,
  P_MEAN, N_MEAN, PRMH_AVE, SITE_ELEV, MgO_Mean, S_Mean,
  SumAve_P, CaO_Mean. Stations must also have at least one
  of the following: AREA_SQKM or LogWSA.

  The data frames are also subject to the following
  constraints: no missing data in either data frame, except
  for the Distinct column; all StationCodes in the "bugs"
  data frame must be represented in the "stations" data
  frame; every SampleID must be associated with only a
  single StationCode; no duplicated data (i.e., every
  concatentation of the SampleID, FinalID, LifeStageCode,
  and Distinct should be unique).

  In order to produce replicable results, the RNG seed can
  be controlled using the \code{rand} argument.
}
\examples{
data(bugs_stations)
results <- CSCI(bugs = bugs_stations[[1]], stations = bugs_stations[[2]])
results$core
}

