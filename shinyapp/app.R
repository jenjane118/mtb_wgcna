library(shiny)
library(shinyjs)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(DT)
library(JBrowseR)

# JBrowseR server for hosting local data files (have to have this running before run app)
data_server <- serve_data(here("shinyapp/sequences"))


# # shinyApp function uses ui object and server function to build shiny app object
shinyApp(ui, server, options=c(launch.browser = .rs.invokeShinyPaneViewer))

# jbrowse still showing in browser mode only
data_server$stop_server()
