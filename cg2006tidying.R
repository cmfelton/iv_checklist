#NOTE: This script does not set the working directory. To 
#read in the data, save it as "census.csv" in the same 
#directory as this script and go Session > Set Working 
#Directory > To Source File Location.

#Obtaining the data requires an IPUMS account. Go to
#https://usa.ipums.org/usa-action/variables/group to
#create an extract using the 1990 Census 5% microdata.
#We list the exact variables used on 
#https://github.com/cmfelton/iv_checklist.

#install.packages("tidyverse")
#install.packages("AER")

library(tidyverse)
library(AER)

#Reading in Census data
census <- read_csv("census.csv")

#The select() line removes variables we don't need

#The filter() line subsets to households with
#only one family present and one married, opposite-sex 
#couple present

#The mutate() line creates a binary indicator for foreign-born
#and recodes the education variable

#The second select() line removes variables we no longer need

#The group_by() line groups the data by household

#The second filter() line subsets to households 
#with no stepchildren and at least two (of the 
#parents' own) children present

census <- census %>%
  dplyr::select(-c(YEAR, SAMPLE, HHWT, CLUSTER, 
                   STRATA, PERWT, RACED, HISPAND, 
                   EDUCD, BPLD, RACESINGD, NFAMS)) %>% 
  filter(GQ == 1 | GQ == 2,
         NSUBFAM == 0, 
         NCOUPLES == 1,
         NMOTHERS == 1,
         NFATHERS == 1) %>% 
  mutate(NOTNATIVE = ifelse(BPL < 100, 0, 
                               ifelse(BPL == 999, NA_integer_, 1)),
         EDUC_R = dplyr::recode(EDUC, 
                         "0" = 0L,
                         "1" = 4L,
                         "2" = 8L,
                         "3" = 9L,
                         "4" = 10L,
                         "5" = 11L,
                         "6" = 12L,
                         "7" = 13L,
                         "8" = 14L,
                         "9" = 15L,
                         "10" = 16L,
                         "11" = 17L)) %>%
  dplyr::select(-c(BPL, EDUC, NMOTHERS, NFATHERS, NCOUPLES, NSUBFAM)) %>%
  group_by(SERIAL) %>%
  filter(!any(RELATED == 303),
         any(NCHILD > 1))

#Creating two vectors of var names we'll need for kids

kid_keep <- c("SERIAL", "PERNUM", "MOMLOC", "POPLOC", "NCHILD", "NSIBS", "RELATED", "SEX",
              "AGE", "BIRTHYR", "RACE", "HISPAN", "SCHOOL", "SCHLTYPE", "NOTNATIVE", "EDUC_R")

kid_vars <- c("MOMLOC", "POPLOC", "NCHILD", "NSIBS", "RELATED", "SEX", "AGE", "BIRTHYR", 
              "RACE", "HISPAN", "SCHOOL", "SCHLTYPE", "NOTNATIVE", "EDUC_R")

#The filter() line subsets to biological 
#children of the household head

#The select() line keeps only the vars in kid_keep

#The group_by() line groups by household

#The arrange() line orders the kids by age,
#using birth year to break (some) ties 
#since age is measured in whole years

#The mutate() line creates a variable
#indicating whether the kid is 1st born,
#2nd born, 3rd born, etc.

#ungroup() ungroups the data

#The next filter() line subsets to first-,
#second-, and third-born children only

#The pivot_wider() line creates a dataframe
#with one household per row and child-specific
#variables labeled with _kid1, _kid2, etc.

#The next filter() line subsets to households
#where the 2nd-born is school-age

#The next filter() line removes households
#where two kids (either 1st- and 2nd-born or 
#2nd- and 3rd-born) are either twins (or part 
#of a set of triplets) or we can't determine 
#the birth order

#The next filter() line subsets to white
#households where the second-born is a boy

#The next mutate() line creates a var
#indicating age difference b/t first-
#and second-born

#The resulting dataframe is saved as a
#new object, kids

kids <- census %>%
  filter(RELATED == 301) %>%
  dplyr::select(all_of(kid_keep)) %>%
  group_by(SERIAL) %>%
  arrange(desc(AGE), BIRTHYR) %>%
  mutate(number_id = row_number(),
         kid_id = paste0("kid", number_id)) %>%
  ungroup() %>%
  filter(number_id < 4) %>%
  pivot_wider(id_cols = SERIAL, names_from = kid_id, values_from = all_of(kid_vars)) %>%
  filter(AGE_kid2 < 18 & AGE_kid2 > 4) %>%
  filter(!((AGE_kid1 == AGE_kid2 & 
              BIRTHYR_kid1 == BIRTHYR_kid2) | 
             (NSIBS_kid1 == 2 & AGE_kid2 == AGE_kid3 & 
                BIRTHYR_kid2 == BIRTHYR_kid3))) %>%
  filter(SEX_kid2 == 1 & RACE_kid2 == 1 & HISPAN_kid2 == 0) %>%
  mutate(age_diff = AGE_kid1 - AGE_kid2)

#Vectors of var names we'll need for parents

parent_keep <- c("SERIAL", "SEX", "AGE", "MARST", "RACE", "HISPAN", "MOMLOC", "POPLOC",
                 "NCHILD", "RELATED", "CHBORN", "NOTNATIVE", "EDUC_R")

parent_vars <- c("AGE", "MARST", "RACE", "HISPAN", "MOMLOC", "POPLOC",
                 "NCHILD", "RELATED", "CHBORN", "NOTNATIVE", "EDUC_R")

#The first filter() line subsets to parents

#The select() line keeps only vars in 
#parent_keep

#The first mutate() line creates a parent_id var

#Next, we group by household and sex

#The next filter() line removes households 
#with non-married parents

#ungroup() ungroups

#pivot_wider creates a dataframe where
#each row represents 1 household
#and vars are labeled _dad or _mom
#mutate creates avg vars needed for 
#2SLS model specification

parents <- census %>%
  filter(RELATED == 101 | RELATED == 201) %>%
  dplyr::select(all_of(parent_keep)) %>%
  mutate(parent_id = ifelse(SEX == 1, "dad", "mom")) %>%
  group_by(SERIAL, SEX) %>%
  filter(!any(MARST > 2)) %>%
  ungroup() %>%
  pivot_wider(id_cols = SERIAL, names_from = parent_id, values_from = all_of(parent_vars)) %>%
  mutate(parent_avg_age = (AGE_mom + AGE_dad) / 2, parent_avg_educ = (EDUC_R_mom + EDUC_R_dad) / 2)

households <- inner_join(kids, parents, by = "SERIAL") 

#should be 1
unique(households$SEX_kid2)

#should be 1 2
unique(households$SEX_kid1)

#should be just TRUE
unique(households$NSIBS_kid2 == households$NSIBS_kid1)

#The mutate() line creates the
#instrument, treatment, and outcome
#vars 

#The filter() line removes households
#where the mother has had fewer children
#than are present in the house
#as well as mixed-race households


households <- households %>%
  mutate(samesex = ifelse(SEX_kid1 == 1, 1, 0),
         priv = ifelse(SCHLTYPE_kid2 == 2, 0, 
                       ifelse(SCHLTYPE_kid2 == 3, 1, NA_integer_)),
         treat = ifelse(NSIBS_kid2 > 1, 1, 0)) %>%
  filter(CHBORN_mom >= NCHILD_mom,
         RACE_mom == 1 & RACE_dad == 1 & HISPAN_mom == 0 & HISPAN_dad == 0)

var_means <- c("priv", "treat", "samesex", "parent_avg_educ", "parent_avg_age",
               "NOTNATIVE_kid2", "AGE_kid2", "age_diff")

#Some descriptive statistics
colMeans(households %>% dplyr::select(all_of(var_means)), na.rm = TRUE)

ivmod <- ivreg(priv ~ treat + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + age_diff +
        factor(AGE_kid2) | samesex + parent_avg_age + parent_avg_educ + NOTNATIVE_kid2 + 
        age_diff + factor(AGE_kid2), data = households)

#Coefficient on treat should be
#-0.057 with estimated SE of 0.0202394
summary(ivmod)

#In Conley and Glauber (2006), 
#estimated effect is -0.046 with
#estimated SE of 0.022 (Table 4
#in Conley and Glauber (2006))

saveRDS(households, "cg2006.rds")
