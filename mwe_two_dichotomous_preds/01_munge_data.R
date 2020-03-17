library(tidyverse)

psw <- read.csv("data/hlv_psw.csv") %>%
    dplyr::select(group = GSSCode,
                  sex, housing, weight) %>%
    mutate(isFemale = as.numeric(sex == "Female"),
           isOwner = as.numeric(housing == "Owns")) %>%
    dplyr::select(group, isOwner, isFemale, weight) %>%
    group_by(group, isOwner, isFemale) %>%
    summarize(w8 = sum(weight))

aux <- read.csv("data/hlv_psw.csv") %>%
    dplyr::select(group = GSSCode,
                  age0, hrsocgrd,
                  weight) %>%
    group_by(group) %>%
    summarize(p.ab = weighted.mean(hrsocgrd == "AB", weight),
              p.over65 = weighted.mean(age0 %in% c("65-74", "75+"), weight))


### Scale these
aux <- aux %>%
    mutate_at(vars(p.ab, p.over65), scale)

psw <- merge(psw, aux, by = "group", sort = FALSE)

psw <- psw %>%
    arrange(group)

saveRDS(psw, file = "working/munged_psw.rds")
