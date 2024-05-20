# data-studios-datasources

## RStudio

### Example 1 - RNASeq

```R
# Adapted from source: https://combine-australia.github.io/RNAseq-R/09-applying-rnaseq-solutions.html

install.packages("BiocManager")
BiocManager::install(c("limma"))
BiocManager::install(c("edgeR"))
BiocManager::install(c("org.Dm.eg.db"))
BiocManager::install(c("gplots"))
BiocManager::install(c("RColorBrewer"))
library(limma)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(org.Dm.eg.db)

counts <- read.delim(file="/workspace/data/datastudios-demo-rstudio/input/2024-05-13/counts_Drosophila.txt")
targets <- read.delim(file="/workspace/data/datastudios-demo-rstudio/input/2024-05-13/SampleInfo_Drosophila.txt")
head(counts)
targets
table(targets$Group)
mycpm <- cpm(counts)
plot(counts[,1],mycpm[,1],xlim=c(0,20),ylim=c(0,5))
abline(v=10,col=2)
abline(h=2,col=4)
thresh <- mycpm > 2
keep <- rowSums(thresh) >= 3
table(keep)
counts.keep <- counts[keep,]
dim(counts.keep)
y <- DGEList(counts.keep)
barplot(y$samples$lib.size)

par(mfrow=c(1,1))
# Get log2 counts per million
logcpm <- cpm(y$counts,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2,outline=FALSE)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs (unnormalised)")

par(mfrow=c(1,2),oma=c(2,0,0,0))
group.col <- c("red","blue")[targets$Group]
boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2,col=group.col,
        pars=list(cex.lab=0.8,cex.axis=0.8))
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs\n(coloured by groups)",cex.main=0.8)

lib.col <- c("light pink","light green")[targets$Library]
boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2, col=lib.col,
        pars=list(cex.lab=0.8,cex.axis=0.8))
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs\n(coloured by library prep)",cex.main=0.8)

par(mfrow=c(1,2))
plotMDS(y,col=group.col)
legend("topright",legend=levels(targets$Group),fill=c("red","blue"))
plotMDS(y,col=lib.col)
legend("topleft",legend=levels(targets$Library),fill=c("light pink","light green"))

logcounts <- cpm(y,log=TRUE)
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]

highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=group.col,scale="row",margins=c(10,5))
```

### Example 2 - Shiny

```R
install.packages("dplyr")
install.packages("gapminder")
install.packages("ggplot2")
install.packages("shiny")

library(shiny)
library(dplyr)
library(ggplot2)
library(gapminder)

# Specify the application port
options(shiny.host = "0.0.0.0")
options(shiny.port = 8180)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      tags$h4("Gapminder Dashboard"),
      tags$hr(),
      selectInput(inputId = "inContinent", label = "Continent", choices = unique(gapminder$continent), selected = "Europe")
    ),
    mainPanel(
      plotOutput(outputId = "outChartLifeExp"),
      plotOutput(outputId = "outChartGDP")
    )
  )
)

server <- function(input, output, session) {
  # Filter data and store as reactive value
  data <- reactive({
    gapminder %>%
      filter(continent == input$inContinent) %>%
      group_by(year) %>%
      summarise(
        AvgLifeExp = round(mean(lifeExp)),
        AvgGdpPercap = round(mean(gdpPercap), digits = 2)
      )
  })

  # Common properties for charts
  chart_theme <- ggplot2::theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

  # Render Life Exp chart
  output$outChartLifeExp <- renderPlot({
    ggplot(data(), aes(x = year, y = AvgLifeExp)) +
      geom_col(fill = "#0099f9") +
      geom_text(aes(label = AvgLifeExp), vjust = 2, size = 6, color = "#ffffff") +
      labs(title = paste("Average life expectancy in", input$inContinent)) +
      theme_classic() +
      chart_theme
  })

  # Render GDP chart
  output$outChartGDP <- renderPlot({
    ggplot(data(), aes(x = year, y = AvgGdpPercap)) +
      geom_line(color = "#f96000", size = 2) +
      geom_point(color = "#f96000", size = 5) +
      geom_label(
        aes(label = AvgGdpPercap),
        nudge_x = 0.25,
        nudge_y = 0.25
      ) +
      labs(title = paste("Average GDP per capita in", input$inContinent)) +
      theme_classic() +
      chart_theme
  })
}

shinyApp(ui = ui, server = server)
```

### Example 3 - Volcano Plot using Shiny

```R
install.packages("shiny")
install.packages("plotly")
install.packages("tidyverse")

library(shiny)
library(plotly)
library(tidyverse)

ui <- fluidPage(
  titlePanel("Volcano Plotly"),
  fluidRow(
    column(
      width = 7,
      plotlyOutput("volcanoPlot", height = "500px")
    ),
    column(
      width = 5,
      dataTableOutput("selectedProbesTable")
    )
  )
)

server <- function(input, output) {

  differentialExpressionResults <-
    read.csv("NKI-DE-results.csv", stringsAsFactors = FALSE) %>%
    mutate(
      probe.type = factor(ifelse(grepl("^Contig", probe), "EST", "mRNA")),
      minusLog10Pvalue = -log10(adj.P.Val),
      tooltip = ifelse(is.na(HUGO.gene.symbol), probe, paste(HUGO.gene.symbol, " (", probe, ")", sep = ""))
    ) %>%
    sample_n(1000)

  output$volcanoPlot <- renderPlotly({

    plot <- differentialExpressionResults %>%
      ggplot(aes(x = logFC,
                 y = minusLog10Pvalue,
                 colour = probe.type,
                 text = tooltip,
                 key = row.names(differentialExpressionResults))) +
      geom_point() +
      xlab("log fold change") +
      ylab("-log10(P-value)")

    plot %>%
      ggplotly(tooltip = "tooltip") %>%
      layout(dragmode = "select")
  })

  output$selectedProbesTable <- renderDataTable({

    eventData <- event_data("plotly_selected")

    selectedData <- differentialExpressionResults %>% slice(0)
    if (!is.null(eventData)) selectedData <- differentialExpressionResults[eventData$key,]

    selectedData %>%
      transmute(
        probe,
        gene = HUGO.gene.symbol,
        `log fold change` = signif(logFC, digits = 2),
        `p-value` = signif(adj.P.Val, digits = 2)
      )
  },
    options = list(dom = "tip", pageLength = 10, searching = FALSE)
  )
}

shinyApp(ui, server, options = list(height = 600))
```

### Example 4 - Load JSON lib, change working directory on startup, read CSV to data frame, export data frame to JSON file

```
install.packages("RJSONIO")
library(RJSONIO)
setwd("/workspace/data")
df <- read.csv("input/2024-01-16/polling_places.csv")
exportJson <- toJSON(df)
write(exportJson, "output/output.json")
```

## JupyterLab

Install all the relevant packages

```
!pip install pandas[pyarrow] jupytext scipy jupyterlab-git qgrid seaborn nb_black
```

Read CSV to data frame, export data frame to JSON file
```
import pandas as pd
df = pd.read_csv('2024-01-16/polling_places.csv', low_memory=False)
df
df.to_json('output.json')
```

## VSCode

TBD
