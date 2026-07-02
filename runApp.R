#!/usr/bin/env Rscript

# Make sure shiny is available
if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("Please install 'shiny' first: install.packages('shiny')")
}

# Always try to launch a browser
options(shiny.launch.browser = TRUE)

# Pick a random free port to avoid 'address already in use'
port <- httpuv::randomPort()

shiny::runApp(
  appDir = ".",
  host = "127.0.0.1",
  port = port,
  launch.browser = TRUE
)
