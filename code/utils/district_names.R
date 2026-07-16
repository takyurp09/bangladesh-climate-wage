## Canonical Bangladesh district name harmonization (GADM / BBS / ERA5)

BBS_COXS_BAZAR <- "Cox\u2019s bazar"

name_recode <- c(
  "bandarban"     = "banderban",
  "bogra"         = "bogura",
  "brahamanbaria" = "brahmmanbaria",
  "chittagong"    = "chattogram",
  "comilla"       = "cumilla",
  "cox's bazar"   = "coxs bazar",
  "coxs bazar"    = "coxs bazar",
  "jessore"       = "jashore",
  "jhalokati"     = "jhalokathi",
  "khagrachhari"  = "khagrachari",
  "maulvibazar"   = "maulavi bazar",
  "nawabganj"     = "chapai nawabganj",
  "netrakona"     = "netrokona",
  "pirojpur"      = "perojpur"
)

normalize_apostrophes <- function(x) {
  x <- gsub("\u2019|\u2018|\u0060", "'", x, perl = TRUE)
  x
}

district_key <- function(x) {
  x <- normalize_apostrophes(trimws(tolower(as.character(x))))
  x <- gsub("['\u2019\u2018]", "", x, perl = TRUE)
  gsub("\\s+", " ", x)
}

clean_name <- function(x) district_key(x)

recode_shp <- function(x) {
  y <- district_key(x)
  mapped <- name_recode[y]
  ifelse(!is.na(mapped), mapped, y)
}

wb_to_bbs <- c(
  "Bandarban"     = "Banderban",
  "Bogra"         = "Bogura",
  "Brahamanbaria" = "Brahmmanbaria",
  "Chittagong"    = "Chattogram",
  "Comilla"       = "Cumilla",
  "Cox's Bazar"   = BBS_COXS_BAZAR,
  "Jessore"       = "Jashore",
  "Jhalokati"     = "Jhalokathi",
  "Khagrachhari"  = "Khagrachari",
  "Maulvibazar"   = "Maulavi Bazar",
  "Nawabganj"     = "Chapai Nawabganj",
  "Netrakona"     = "Netrokona",
  "Pirojpur"      = "Perojpur"
)

map_wb_district <- function(zila, bbs_districts = NULL) {
  if (!is.null(bbs_districts)) {
    cox <- grep("^Cox", bbs_districts, value = TRUE)
    if (length(cox) == 1L) wb_to_bbs["Cox's Bazar"] <- cox
  }
  out <- wb_to_bbs[zila]
  ifelse(is.na(out), zila, unname(out))
}
