---
title: "Metadata Cleaning and Object Creation"
output: html_document
date: "2025-02-03"
---
        
```{r}
library(tidyverse)
```

```{r}
metadata <- read.csv("metadata.csv")
head(metadata)
str(metadata)
```
```{r}
# Check for missing values in each column
colSums(is.na(metadata))
```
```{r}
# Keep relevant columns
# Using column names
metadata_subset <- metadata[, c("Run", "BioProject", "Age", "Sex", "BMI", "Country", "Continent", "Ethnicity", "Study.Group")] 
```

```{r}
head(metadata_subset)
```

```{r}
# save dataframe
write.csv(metadata_subset, "metadata.csv", row.names=FALSE)

# save as rds
saveRDS(metadata_subset, "metadata.rds")
```