library(tidyverse)

psw <- read.csv("data/hlv_psw.csv") %>%
    dplyr::select(group = GSSCode,
                  age = age0, education, weight) %>%
    dplyr::select(group, age, education, weight) %>%
    group_by(group, age, education) %>%
    summarize(w8 = sum(weight))

aux <- read.csv("data/hlv_psw.csv") %>%
    dplyr::select(group = GSSCode,
                  housing, hrsocgrd,
                  weight) %>%
    group_by(group) %>%
    summarize(p.ab = weighted.mean(hrsocgrd == "AB", weight),
              p.owns = weighted.mean(housing == "Owns", weight))


### Scale these
aux <- aux %>%
    mutate_at(vars(p.ab, p.owns), scale)

psw <- merge(psw, aux, by = "group", sort = FALSE)

psw <- psw %>%
    arrange(group)

saveRDS(psw, file = "working/munged_psw.rds")
