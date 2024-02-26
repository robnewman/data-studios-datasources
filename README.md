# data-studios-datasources

## RStudio

```
# return list of demos
demos()
# return list of built-in datasets
data()
# access built-in dataset by its name
mtcars
```

Load JSON lib, change working directory on startup, read CSV to data frame, export data frame to JSON file

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
