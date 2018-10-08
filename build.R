#! /usr/bin/env Rscript 

rsconnect::setAccountInfo(
  name='saladi',
  token='23E9F0422945E3109A0125D66D28CBDB',
  secret=Sys.getenv("SHINYAPPS_SECRET"))

rsconnect::deployApp()

