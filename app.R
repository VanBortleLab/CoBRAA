
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(plotly)
library(dplyr)
library(Matrix)
# library(shinycssloaders)  ##For spinner

source("utils.R")


# ============================================================================
# LOAD DATA AT STARTUP
# ============================================================================

# Load protein data indices
protein_index <<- readRDS("data/tem_chunks/protein_index.rds")
protein_names <<- protein_index$protein
# Load PPI matrix (sparse symmetric)
PPI_sparse <<- readRDS("data/PPI_sparse_sym.rds")
ppi_names <<- rownames(PPI_sparse)
# Load CCI matrix (sparse symmetric)
CCI_sparse <<- readRDS("data/CCI_sparse_sym.rds")
cci_names <<- rownames(CCI_sparse)
# Load BioGRID dbinary matrix 
Biogrid_binary <<- readRDS("data/Biogrid_binary.rds")
# Caches for tem_chunks
z_cache <<- new.env()
rz_cache <<- new.env()




# ============================================================================
# UI
# ============================================================================

ui <- dashboardPage(
  
  dashboardHeader(title = "Protein & Complex Analysis"),
  
  dashboardSidebar(
    
    # Mode selection
    radioButtons(
      "mode",
      "Select Mode:",
      choices = c("Protein" = "protein", "Complex" = "complex"),
      selected = "protein"
    ),
    
    hr(),
    
    # PROTEIN MODE OPTIONS
    conditionalPanel(
      condition = "input.mode == 'protein'",
      
      h4("Protein Analysis"),
      
      radioButtons(
        "protein_submode",
        "Analysis Type:",
        choices = c(
          "Popularity Corrected Interactors" = "popularity",
          "Reciprocal Analysis" = "reciprocal",
          "Dominant Plot" = "dominant"
        ),
        selected = "popularity"
      ),
      
      # Popularity inputs
      conditionalPanel(
        condition = "input.protein_submode == 'popularity'",
        selectizeInput(
          "protein_pop_poi",
          "Select Protein:",
          choices = NULL,
          options = list(maxOptions = 10)
        ),
        sliderInput(
          "protein_pop_labels",
          "Number of Labels:",
          min = 0, max = 200, value = 50, step = 10
        ),
        actionButton("protein_pop_plot", "Generate Plot", class = "btn-primary"),
        br(), br(),
        downloadButton("protein_pop_download", "Download Data (.tsv)")  # ← ADD THIS LINE
      ),
      
      # Reciprocal inputs
      # Reciprocal inputs
      conditionalPanel(
        condition = "input.protein_submode == 'reciprocal'",
        selectizeInput(
          "protein_recip_poi",
          "Select Protein:",
          choices = NULL,
          options = list(maxOptions = 10)
        ),
        sliderInput(
          "protein_recip_threshold",
          "Label Threshold:",
          min = 0, max = 1, value = 0.5, step = 0.01
        ),
        actionButton("protein_recip_plot", "Generate Plot", class = "btn-primary"),
        br(), br(),
        downloadButton("protein_recip_download", "Download Data (.tsv)")  # ← ADD THIS
      ),
      
      # Dominant plot inputs
      conditionalPanel(
        condition = "input.protein_submode == 'dominant'",
        selectizeInput("protein_dom_sel1", "Select 1 (or NULL for triangle):", 
                       choices = NULL, options = list(maxOptions = 10)),
        selectizeInput("protein_dom_sel2", "Select 2:", 
                       choices = NULL, options = list(maxOptions = 10)),
        selectizeInput("protein_dom_sel3", "Select 3:", 
                       choices = NULL, options = list(maxOptions = 10)),
        selectizeInput("protein_dom_sel4", "Select 4:", 
                       choices = NULL, options = list(maxOptions = 10)),
        
        sliderInput("protein_dom_entropy_min", "Entropy Min:", 
                    min = 0, max = 2, value = 0, step = 0.1),
        sliderInput("protein_dom_entropy_max", "Entropy Max:", 
                    min = 0, max = 2, value = 1, step = 0.1),
        sliderInput("protein_dom_maxvalue", "Max Value Min:", 
                    min = 2, max = 50, value = 20, step = 1),
        
        actionButton("protein_dom_plot", "Generate Plot", class = "btn-primary"),
        br(), br(),
        downloadButton("protein_dom_download", "Download Data")
      )
    ),
    
    # COMPLEX MODE OPTIONS
    conditionalPanel(
      condition = "input.mode == 'complex'",
      
      h4("Complex Analysis"),
      
      radioButtons(
        "complex_submode",
        "Analysis Type:",
        choices = c(
          "Network Proximal Complexes" = "network",
          "Dominant Plot" = "dominant"
        ),
        selected = "network"
      ),
      
      # Network inputs
      conditionalPanel(
        condition = "input.complex_submode == 'network'",
        selectizeInput(
          "complex_net_coi",
          "Select Complex:",
          choices = NULL,
          options = list(maxOptions = 10)
        ),
        sliderInput(
          "complex_net_labels",
          "Number of Labels:",
          min = 0, max = 500, value = 10, step = 50
        ),
        actionButton("complex_net_plot", "Generate Plot", class = "btn-primary"),
        br(), br(),
        downloadButton("complex_net_download", "Download Data (.tsv)")
      ),
      
      # Dominant plot inputs
      conditionalPanel(
        condition = "input.complex_submode == 'dominant'",
        selectizeInput("complex_dom_sel1", "Select 1 (or NULL for triangle):", 
                       choices = NULL, options = list(maxOptions = 10)),
        selectizeInput("complex_dom_sel2", "Select 2:", 
                       choices = NULL, options = list(maxOptions = 10)),
        selectizeInput("complex_dom_sel3", "Select 3:", 
                       choices = NULL, options = list(maxOptions = 10)),
        selectizeInput("complex_dom_sel4", "Select 4:", 
                       choices = NULL, options = list(maxOptions = 10)),
        
        sliderInput("complex_dom_entropy_min", "Entropy Min:", 
                    min = 0, max = 2, value = 0, step = 0.1),
        sliderInput("complex_dom_entropy_max", "Entropy Max:", 
                    min = 0, max = 2, value = 1, step = 0.1),
        sliderInput("complex_dom_maxvalue", "Max Value Min:", 
                    min = 0, max = 100, value = 20, step = 5),
        
        actionButton("complex_dom_plot", "Generate Plot", class = "btn-primary"),
        br(), br(),
        downloadButton("complex_dom_download", "Download Data")
      )
    )
  ),
  
  dashboardBody(
    
    fluidRow(
      box(
        width = 12,
        title = textOutput("plot_title"),
        status = "primary",
        solidHeader = TRUE,
        plotlyOutput("main_plot", height = "700px")
      )
    ),
    
    fluidRow(
      box(
        width = 12,
        title = "Data Table",
        status = "info",
        solidHeader = TRUE,
        collapsible = TRUE,
        collapsed = TRUE,
        tableOutput("data_table")
      )
    )
  )
)











# ============================================================================
# SERVER
# ============================================================================


server <- function(input, output, session) {
  
  # Reactive values
  plot_data <- reactiveVal(NULL)
  
  # -------------------------------------------------------------------------
  # UPDATE SELECTIZE INPUTS BASED ON MODE
  # -------------------------------------------------------------------------
  
  observe({
    if (input$mode == "protein") {
      # Update protein mode dropdowns
      updateSelectizeInput(session, "protein_pop_poi", choices = protein_names, server = TRUE)
      updateSelectizeInput(session, "protein_recip_poi", choices = protein_names, server = TRUE)
      
      updateSelectizeInput(session, "protein_dom_sel1", 
                           choices = c("NULL", ppi_names), 
                           selected = ppi_names[1], server = TRUE)
      updateSelectizeInput(session, "protein_dom_sel2", 
                           choices = ppi_names, selected = ppi_names[2], server = TRUE)
      updateSelectizeInput(session, "protein_dom_sel3", 
                           choices = ppi_names, selected = ppi_names[3], server = TRUE)
      updateSelectizeInput(session, "protein_dom_sel4", 
                           choices = ppi_names, selected = ppi_names[4], server = TRUE)
      
    } else {
      # Update complex mode dropdowns
      updateSelectizeInput(session, "complex_net_coi", choices = cci_names, server = TRUE)
      
      updateSelectizeInput(session, "complex_dom_sel1", 
                           choices = c("NULL", cci_names), 
                           selected = cci_names[1], server = TRUE)
      updateSelectizeInput(session, "complex_dom_sel2", 
                           choices = cci_names, selected = cci_names[2], server = TRUE)
      updateSelectizeInput(session, "complex_dom_sel3", 
                           choices = cci_names, selected = cci_names[3], server = TRUE)
      updateSelectizeInput(session, "complex_dom_sel4", 
                           choices = cci_names, selected = cci_names[4], server = TRUE)
    }
  })
  
  # -------------------------------------------------------------------------
  # DYNAMIC PLOT TITLE
  # -------------------------------------------------------------------------
  
  output$plot_title <- renderText({
    if (input$mode == "protein") {
      switch(input$protein_submode,
             "popularity" = "Popularity Corrected Interactors",
             "reciprocal" = "Reciprocal Analysis",
             "dominant" = "Dominant Plot (Tetrahedron/Triangle)")
    } else {
      switch(input$complex_submode,
             "network" = "Network Proximal Complexes",
             "dominant" = "Dominant Plot (Tetrahedron/Triangle)")
    }
  })
  
  # -------------------------------------------------------------------------
  # PROTEIN MODE: POPULARITY CORRECTED INTERACTORS
  # -------------------------------------------------------------------------
  
  # Create reactive value for download data
  protein_pop_data <- reactiveVal(NULL)
  
  observeEvent(input$protein_pop_plot, {
    req(input$protein_pop_poi)
    
    # Generate the plot
    p <- get_popularity_corrected_plot(
      poi = input$protein_pop_poi,
      label_top_n = input$protein_pop_labels
    )
    
    output$main_plot <- renderPlotly(p)
    
    # Prepare table data separately
    tu1 <- sort(PPI_sparse[input$protein_pop_poi, ], decreasing = TRUE)
    biogrid_row <- Biogrid_binary[input$protein_pop_poi, ]
    
    table_data <- data.frame(
      Rank = 1:length(tu1),
      Protein = names(tu1),
      Score = as.numeric(tu1),
      In_BioGRID = as.logical(biogrid_row[names(tu1)]),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        Category = ifelse(
          Score >= 2 & In_BioGRID,
          "Popularity-corrected",
          ifelse(Score >= 2 & !In_BioGRID, "Predicted", "Unassociated")
        )
      ) %>%
      head(100)  # Show top 100 in table
    
    output$data_table <- renderTable(
      table_data,
      striped = TRUE, bordered = TRUE, hover = TRUE
    )
    
    # Store for download
    protein_pop_data(table_data)
  })
  
  # Download handler
  output$protein_pop_download <- downloadHandler(
    filename = function() {
      paste0("PPI_", input$protein_pop_poi, "_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      req(protein_pop_data())
      write.table(protein_pop_data(), file, sep = "\t", 
                  row.names = FALSE, quote = FALSE)
    }
  )
  
  
  # -------------------------------------------------------------------------
  # PROTEIN MODE: RECIPROCAL ANALYSIS
  # -------------------------------------------------------------------------
  # Create reactive value for download data (put this near top of server, before observeEvents)
  protein_recip_data <- reactiveVal(NULL)
  
  observeEvent(input$protein_recip_plot, {
    req(input$protein_recip_poi)
    
    withProgress(message = 'Reciprocal Analysis : ', value = 0, {
      
      setProgress(0.75, detail = "Reciprocating 🔄 .....")
      
      # Get reciprocal plot data
      recip_data <- get_info_from_reciprocal_plot(
        prot = input$protein_recip_poi,
        label_threshold = input$protein_recip_threshold
      )
      
      # Create plot
      p <- get_reciprocal_plot(recip_data, input$protein_recip_poi)
      
      output$main_plot <- renderPlotly(p)
      
      # Prepare download data
      download_data <- recip_data %>% 
        arrange(desc(distances)) %>% 
        transmute(Protein = protein, Category = category, 
                  Z = self, Reciprocal_Z = other, Deviation = distances, 
                  PPI_Score = ppi_score, In_BioGRID = in_biogrid)
      
      # Show table (top 100)
      output$data_table <- renderTable(
        head(download_data, 100),
        striped = TRUE, bordered = TRUE, hover = TRUE
      )
      
      # Store for download
      protein_recip_data(download_data)
      
      setProgress(1, detail = "Complete!")
    })
  })
  
  # Add download handler (put this after the observeEvent above)
  output$protein_recip_download <- downloadHandler(
    filename = function() {
      paste0("Reciprocal_", input$protein_recip_poi, "_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      req(protein_recip_data())
      write.table(protein_recip_data(), file, sep = "\t", 
                  row.names = FALSE, quote = FALSE)
    }
  )
  
  
  
  # -------------------------------------------------------------------------
  # PROTEIN MODE: DOMINANT PLOT
  # -------------------------------------------------------------------------
  
  # Near top of server, with other reactive values
  protein_dom_data <- reactiveVal(NULL)
  
  
  observeEvent(input$protein_dom_plot, {
    req(input$protein_dom_sel1, input$protein_dom_sel2, 
        input$protein_dom_sel3, input$protein_dom_sel4)
    
    withProgress(message = '🔺 Building tetrahedron...', value = 0, {
      
      setProgress(0.3, detail = "Extracting protein interactions...")
      
      # Build list of selections
      selections <- c(input$protein_dom_sel1, input$protein_dom_sel2,
                      input$protein_dom_sel3, input$protein_dom_sel4)
      
      # Extract non-NULL columns first
      non_null_selections <- selections[selections != "NULL"]
      
      if (length(non_null_selections) > 0) {
        selected_cols <- PPI_sparse[, non_null_selections, drop = FALSE]
        
        # Convert to regular matrix
        if (inherits(selected_cols, "Matrix")) {
          selected_cols <- as.matrix(selected_cols)
        }
      } else {
        # All NULL - create empty matrix
        selected_cols <- matrix(0, nrow = nrow(PPI_sparse), ncol = 0)
      }
      
      # Add NULL columns where needed
      final_cols <- matrix(nrow = nrow(selected_cols), ncol = 4)
      colnames(final_cols) <- selections
      
      col_index <- 1
      for (i in 1:4) {
        if (selections[i] == "NULL") {
          final_cols[, i] <- 0  # NULL column = all zeros
        } else {
          final_cols[, i] <- selected_cols[, col_index]
          col_index <- col_index + 1
        }
      }
      
      # Add rownames
      rownames(final_cols) <- rownames(PPI_sparse)
      
      setProgress(0.6, detail = "Computing center of mass...")
      
      # Call your plotTetrahedron function
      result <- plotTetrahedron(
        data = final_cols,
        entropyrange = c(input$protein_dom_entropy_min, input$protein_dom_entropy_max),
        maxvaluerange = c(input$protein_dom_maxvalue, Inf),
        output_table = TRUE
      )
      
      setProgress(0.9, detail = "Rendering 3D plot...")
      
      # Render the plot
      output$main_plot <- renderPlotly(result$plot)
      
      # Show the table
      output$data_table <- renderTable(
        result$table,
        striped = TRUE, bordered = TRUE, hover = TRUE
      )
      
      # Store for download
      protein_dom_data(result$table)
      
      setProgress(1, detail = "Complete!")
    })
  })
  
  # Download handler for dominant plot
  output$protein_dom_download <- downloadHandler(
    filename = function() {
      paste0("Dominant_Plot_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      req(protein_dom_data())
      write.table(protein_dom_data(), file, sep = "\t", 
                  row.names = FALSE, quote = FALSE)
    }
  )


  
  
  
  # -------------------------------------------------------------------------
  # COMPLEX MODE: NETWORK PROXIMAL COMPLEXES
  # -------------------------------------------------------------------------
  complex_net_data <- reactiveVal(NULL)
  
  observeEvent(input$complex_net_plot, {
    req(input$complex_net_coi)
    
    withProgress(message = paste('🕸️ Analyzing', input$complex_net_coi, 'network...'), value = 0, {
      
      setProgress(0.5, detail = "Extracting complexes...")
      
      # Get the row for complex of interest
      coi_scores <- sort(CCI_sparse[input$complex_net_coi, ])
      
      # Create data frame
      plot_data <- data.frame(
        complex = names(coi_scores),
        index = 1:length(coi_scores),
        score = as.numeric(coi_scores),
        log_score = log(as.numeric(coi_scores) + 1),  # log(score + 1)
        stringsAsFactors = FALSE
      )
      
      # Label only top N complexes
      n_complexes <- nrow(plot_data)
      label_top_n <- input$complex_net_labels  # Use slider value
      
      plot_data$label <- ifelse(
        plot_data$index > (n_complexes - label_top_n),
        plot_data$complex,
        ""
      )
      
      setProgress(0.8, detail = "Creating plot...")
      
      # Calculate sparse y-axis ticks
      y_range <- range(plot_data$log_score, na.rm = TRUE)
      y_min <- floor(y_range[1])
      y_max <- ceiling(y_range[2])
      tick_interval <- max(1, ceiling((y_max - y_min) / 10))
      y_ticks <- seq(y_min, y_max, by = tick_interval)
      
      # Create plot
      p <- plot_ly(plot_data) %>%
        add_trace(
          x = ~index,
          y = ~log_score,
          type = "scattergl",
          mode = "markers+text",
          marker = list(
            color = "rgb(100,150,200)",
            size = 10
          ),
          text = ~label,
          textposition = "top center",
          textfont = list(size = 9, color = "black"),
          hovertext = ~paste0(
            "<b>", complex, "</b><br>",
            "Rank: ", index, " / ", n_complexes, "<br>",
            "Score: ", round(score, 4), "<br>",
            "Log(Score+1): ", round(log_score, 4)
          ),
          hoverinfo = "text",
          showlegend = FALSE
        ) %>%
        # Add horizontal threshold line at y = log(2 + 1)
        add_segments(
          x = 1,
          xend = n_complexes,
          y = log(3),
          yend = log(3),
          line = list(color = "red", width = 2, dash = "dash"),
          showlegend = FALSE,
          hoverinfo = "skip"
        ) %>%
        layout(
          title = list(
            text = paste("Network Proximal Complexes for", input$complex_net_coi,
                         sprintf("(Top %d labeled)", label_top_n)),
            font = list(size = 16)
          ),
          xaxis = list(
            title = "",
            showgrid = FALSE,
            showticklabels = FALSE,
            zeroline = FALSE,
            showline = FALSE
          ),
          yaxis = list(
            title = "Log(CCI Score + 1)",
            showgrid = TRUE,
            gridcolor = "rgba(200,200,200,0.2)",
            showticklabels = TRUE,
            tickmode = "array",
            tickvals = y_ticks,
            tickformat = ".0f",
            zeroline = FALSE
          ),
          hovermode = "closest",
          plot_bgcolor = "rgba(255,255,255,1)",
          dragmode = "pan"
        ) %>%
        config(
          scrollZoom = TRUE,
          displayModeBar = TRUE,
          modeBarButtonsToRemove = c("lasso2d", "select2d")
        )
      
      output$main_plot <- renderPlotly(p)
      
      # Prepare download data
      download_data <- plot_data %>%
        arrange(desc(score)) %>%
        select(Rank = index, Complex = complex, Score = score, Log_Score = log_score)
      
      # Show table (top 100)
      output$data_table <- renderTable(
        head(download_data, 100),
        striped = TRUE, bordered = TRUE, hover = TRUE
      )
      
      # Store for download
      complex_net_data(download_data)
      
      setProgress(1, detail = "Complete!")
    })
  })
  
  # Download handler
  output$complex_net_download <- downloadHandler(
    filename = function() {
      paste0("Network_Proximal_", input$complex_net_coi, "_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      req(complex_net_data())
      write.table(complex_net_data(), file, sep = "\t", 
                  row.names = FALSE, quote = FALSE)
    }
  )
  
  # -------------------------------------------------------------------------
  # COMPLEX MODE: DOMINANT PLOT
  # -------------------------------------------------------------------------
  complex_dom_data <- reactiveVal(NULL)
  
  observeEvent(input$complex_dom_plot, {
    req(input$complex_dom_sel1, input$complex_dom_sel2, 
        input$complex_dom_sel3, input$complex_dom_sel4)
    
    withProgress(message = '🔺 Building complex tetrahedron...', value = 0, {
      
      setProgress(0.3, detail = "Extracting complex interactions...")
      
      # Build list of selections
      selections <- c(input$complex_dom_sel1, input$complex_dom_sel2,
                      input$complex_dom_sel3, input$complex_dom_sel4)
      
      # Extract non-NULL columns first
      non_null_selections <- selections[selections != "NULL"]
      
      if (length(non_null_selections) > 0) {
        selected_cols <- CCI_sparse[, non_null_selections, drop = FALSE]
        
        # Convert to regular matrix
        if (inherits(selected_cols, "Matrix")) {
          selected_cols <- as.matrix(selected_cols)
        }
      } else {
        # All NULL - create empty matrix
        selected_cols <- matrix(0, nrow = nrow(CCI_sparse), ncol = 0)
      }
      
      # Add NULL columns where needed
      final_cols <- matrix(nrow = nrow(selected_cols), ncol = 4)
      colnames(final_cols) <- selections
      
      col_index <- 1
      for (i in 1:4) {
        if (selections[i] == "NULL") {
          final_cols[, i] <- 0  # NULL column = all zeros
        } else {
          final_cols[, i] <- selected_cols[, col_index]
          col_index <- col_index + 1
        }
      }
      
      # Add rownames
      rownames(final_cols) <- rownames(CCI_sparse)
      
      setProgress(0.6, detail = "Computing center of mass...")
      
      # Call your plotTetrahedron function
      result <- plotTetrahedron(
        data = final_cols,
        entropyrange = c(input$complex_dom_entropy_min, input$complex_dom_entropy_max),
        maxvaluerange = c(input$complex_dom_maxvalue, Inf),
        output_table = TRUE
      )
      
      setProgress(0.9, detail = "Rendering 3D plot...")
      
      # Render the plot
      output$main_plot <- renderPlotly(result$plot)
      
      # Show the table
      output$data_table <- renderTable(
        result$table,
        striped = TRUE, bordered = TRUE, hover = TRUE
      )
      
      # Store for download
      complex_dom_data(result$table)
      
      setProgress(1, detail = "Complete!")
    })
  })
  
  # Download handler for complex dominant plot
  output$complex_dom_download <- downloadHandler(
    filename = function() {
      paste0("Complex_Dominant_Plot_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      req(complex_dom_data())
      write.table(complex_dom_data(), file, sep = "\t", 
                  row.names = FALSE, quote = FALSE)
    }
  )
  
}

# ============================================================================
# RUN APP
# ============================================================================

shinyApp(ui, server)