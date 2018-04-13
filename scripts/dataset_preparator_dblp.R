

# SET THESE PARAMETERS:
# (don't forget to set the working directory)

project_dir <- "results/dblp1/"

files_prefix <- "results_"

# Prepared files will be created in the current working directory

# ---------------------

# uncomment the following three lines and run them to install the packages, if necessary
# install.packages("foreign")
# install.packages("dplyr")
# install.packages("reshape2")


library(foreign)
library(dplyr)
library(reshape2)

# ---------------------

read_stats <- function(directory, filename) {
  filename <- paste0(directory, filename)
  if (!file.exists(filename)) {
    return(NULL)
  }
  stat_lines <- readLines(filename)
  new_stats <- list()
  for (l in stat_lines) {
    colon_pos <- regexpr(':', l)[1]
    key <- substring(l, 1, colon_pos-1)
    value <- substring(l, colon_pos+1)
    new_stats[[key]] <- value
  }
  return(new_stats)
}


# THIS IS SPECIFIC TO DATA:
instance_ids <- 0:199
max_pos_id <- 99
experiments <- 0:19


df_train <- NULL
df_test <- NULL

for (i in seq_along(experiments)) {
  stats <- read_stats(project_dir, paste0(files_prefix, experiments[i], ".txt"))  
  pos_train <- parse_evaluation(stats[["TRAIN_EVALUATION_POSITIVE"]])
  neg_train <- parse_evaluation(stats[["TRAIN_EVALUATION_NEGATIVE"]])
  pos_test <- parse_evaluation(stats[["TEST_EVALUATION_POSITIVE"]])
  neg_test <- parse_evaluation(stats[["TEST_EVALUATION_NEGATIVE"]])
  
  if (is.null(df_test)) {
    df_train <- data.frame(class = rep(c("pos", "neg"), c(length(pos_train), length(neg_train))))
    df_test <- data.frame(class = rep(c("pos", "neg"), c(length(pos_test), length(neg_test))))
  }
  df_train[[paste0("pattern_", i)]] <- c(pos_train, neg_train)
  df_test[[paste0("pattern_", i)]] <- c(pos_test, neg_test)
}


df_train <- select(df_train, starts_with("pattern"), class)
df_test <- select(df_test, starts_with("pattern"), class)


write.arff(df_train, "dblp_ewaldis_train.arff")
write.arff(df_test, "dblp_ewaldis_test.arff")


# BASELINE

baseline_data_train <- read.csv(paste0(project_dir, "baseline_edges_data_train.csv"))
baseline_data_test <- read.csv(paste0(project_dir, "baseline_edges_data_test.csv"))

baseline_data_train <- unique(baseline_data_train)
baseline_data_test <- unique(baseline_data_test)

baseline_data_train <- baseline_data_train[baseline_data_train$timestamp > 0, ]
baseline_data_test <- baseline_data_test[baseline_data_test$timestamp > 0, ]


data_wide_train <- dcast(baseline_data_train, id + class ~ new_label, value.var="new_label",
                         fun.aggregate = length)

data_wide_test <- dcast(baseline_data_test, id + class ~ new_label, value.var="new_label",
                        fun.aggregate = length)

# remove extra columns from test and append missing with zeros (those that are only in train)
data_wide_test <- data_wide_test[, names(data_wide_test) %in% names(data_wide_train)]
missing_columns <- names(data_wide_train)[!(names(data_wide_train) %in% names(data_wide_test))]
data_wide_test[,missing_columns] <- 0

# add missing ids (they are not in data -> they did not have any of those edges -> fill with zeros)
missing_ids_train <- instance_ids[!(instance_ids %in% data_wide_train$id)]
missing_ids_test <- instance_ids[!(instance_ids %in% data_wide_test$id)]

x_train = data.frame(id = missing_ids_train, class = ifelse(missing_ids_train <= max_pos_id, "pos", "neg"))
x_train[, names(data_wide_train)[c(-1, -2)]] <- 0
data_wide_train <- rbind(data_wide_train, x_train)

x_test = data.frame(id = missing_ids_test, class = ifelse(missing_ids_test <= max_pos_id, "pos", "neg"))
x_test[, names(data_wide_test)[c(-1, -2)]] <- 0
data_wide_test <- rbind(data_wide_test, x_test)


# remove id
data_wide_train <- data_wide_train[, -1]
data_wide_test <- data_wide_test[, -1]

data_wide_train <- data_wide_train[ , order(names(data_wide_train))]
data_wide_test <- data_wide_test[ , order(names(data_wide_test))]

data_wide_train <- select(data_wide_train, starts_with("f"), class)
data_wide_test <- select(data_wide_test, starts_with("f"), class)

write.arff(data_wide_train, "dblp_baseline_train.arff")
write.arff(data_wide_test, "dblp_baseline_test.arff")





