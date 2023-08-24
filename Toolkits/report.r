---
title: 'Parameters Report'
author: "Mahe Chen"
date: "June, 2023"
output:
  beamer_presentation:
    theme: AnnArbor
    colortheme: orchid
    fonttheme: structurebold
    slide_level: 2
    includes:
      in_header: preamble.tex
  slidy_presentation: default
classoption: aspectratio=169
fontsize: 12pt
urlcolor: blue
---

```{r setup, include=FALSE}
library(knitr)
opts_chunk$set(echo = FALSE, out.width = "25%")
layout_matrix <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)

include_graphics_row <- function(...) {
  layout(layout_matrix)
  for (img_path in c(...)) {
    knitr::include_graphics(img_path)
  }
}

create_summary <- function(column1_values, column2_values, column3_values) {
  data <- data.frame(Column1 = column1_values, Column2 = column2_values, Column3 = column3_values)
  colnames(data) <- c("min_corr", "min_pnr", "neurons")
  print(data)
}
```

## Row 1: 0.75-0.77, Row 2: 0.78-0.8
```{r, echo=FALSE, out.width="25%"}
include_graphics_row("075img.jpg", "076img.jpg", "077img.jpg", "078img.jpg", "079img.jpg", "08img.jpg")
```

## Row 1: 0.81-0.83, Row 2: 0.84-0.86
```{r, echo=FALSE, out.width="25%"}
include_graphics_row("081img.jpg", "082img.jpg", "083img.jpg", "084img.jpg", "085img.jpg", "086img.jpg")
```

## Row 1: 0.87-0.89, Row 2: 0.9-0.92
```{r, echo=FALSE, out.width="25%"}
include_graphics_row("087img.jpg", "088img.jpg", "089img.jpg", "09img.jpg", "091img.jpg", "092img.jpg")
```

## Row 1: 0.93-0.94, Row 2: -
```{r, echo=FALSE, out.width="25%"}
include_graphics_row("093img.jpg", "094img.jpg")
```

## Row 1: 0.6,6  0.7,5  0.7,6 | Row 2: 0.7,8  0.7,9
```{r, echo=FALSE, out.width="25%"}
include_graphics_row("06_6img.jpg", "07_5img.jpg", "07_6img.jpg", "07_8img.jpg", "07_9img.jpg")
```

## Summary I
```{r, echo=FALSE}
create_summary(c(0.6, 0.7, 0.7, 0.7, 0.7, seq(0.75, 0.81, by = 0.01)),
               c(6, 5, 6, 8, 9, seq(7.5, 8.1, by = 0.1)),
               c(193, 139, 139, 139, 139, 128, 127, 128, 126, 124, 121, 121))
```

# Summary II
```{r, echo=FALSE}
create_summary(seq(0.82, 0.93, by = 0.01),
               seq(8.2, 9.3, by = 0.1),
               c(119, 120, 122, 122, 117, 118, 113, 112, 114, 115, 112, 113))
```

## Summary III
```{r, echo=FALSE}
create_summary(0.94, 9.4, 112)
```