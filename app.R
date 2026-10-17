
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(plotly)
library(dplyr)
library(plotly.repel)
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

ABOUT_HTML <- paste0(
  "<h1>Welcome to CoBBRA</h1>\n",
  "<p><strong>CoBBRA (Conductance-Based Reciprocal Association Analysis)</strong> is an interactive platform for exploring and visualizing protein–protein interactions (PPIs) and complex–complex interactions (CCIs). The application includes tools to identify interaction partners and statistically meaningful network relationships between individual protein pairs.</p>\n",
  "<p>CoBBRA offers two analysis modes: <strong>Protein Mode</strong> for protein–protein interaction analysis, and <strong>Complex Mode</strong> for macromolecular interaction analysis and prediction.</p>\n",
  "<h2>1. Protein Mode</h2>\n",
  "<h3>1.1 Popularity-Corrected Interaction Analysis</h3>\n",
  "<p>This tool identifies and visualizes potential interaction partners for a selected protein, ranked according to their PPI scores. Interactions with a PPI score ≥ 2 are highlighted and classified into two categories:</p>\n",
  "<ul><li><strong><span style=\"color:#487A9C\">Blue – Popularity-Corrected Interactions:</span></strong> Previously reported interactions in BioGRID that retain a PPI score ≥ 2 after correcting for protein popularity bias.</li><li><strong><span style=\"color:#FB8451\">Orange – Statistically Meaningful Network Relationships:</span></strong> These indicate greater-than-expected connectivity between protein pairs, despite the absence of previously reported interactions in BioGRID.</li></ul>\n",
  "<p><strong>Interactive parameter:</strong></p>\n",
  "<ul><li><strong>Number of Labels:</strong> Controls how many of the top-ranked interactions are labeled in the plot.</li></ul>\n",
  "<p>This analysis helps distinguish established interactions that remain significant after popularity correction from potentially novel interaction partners.</p>\n",
  "<h3>1.2 Reciprocal Analysis</h3>\n",
  "<p>This tool evaluates reciprocal interaction relationships for a selected protein through an interactive scatter plot. The reciprocal plot visualizes the relationship between interaction scores in both directions, allowing users to identify proteins that deviate from the overall interaction trend for a given protein.</p>\n",
  "<p>Proteins are colored by their association classification:</p>\n",
  "<ul><li><strong><span style=\"color:#487A9C\">Blue – Popularity-Corrected Interactions:</span></strong> Previously reported interactions in BioGRID that retain a PPI score ≥ 2 after correcting for protein popularity bias.</li><li><strong><span style=\"color:#FB8451\">Orange – Statistically Meaningful Network Relationships:</span></strong> These indicate greater-than-expected connectivity between protein pairs, despite the absence of previously reported interactions in BioGRID.</li></ul>\n",
  "<p><em>Note: The reciprocal analysis plot shown here represents only half of a reciprocal association analysis. For any protein pair, significance requires deviation from the expected trend from the perspective of both protein X and protein Y. Blue and orange indicate pairs that meet the two-direction criterion.</em></p>\n",
  "<p><strong>Interactive parameter:</strong></p>\n",
  "<ul><li><strong>Label Threshold:</strong> Controls which interactors are labeled based on their vertical distance above the fitted trendline. Adjusting this threshold helps focus on proteins exhibiting larger deviations.</li></ul>\n",
  "<h3>1.3 Feature Dominance</h3>\n",
  "<p>Feature Dominance analysis explores the relative interaction strengths of proteins with respect to three or four selected proteins of interest. Users can select up to four reference proteins, positioned at the vertices of a tetrahedron. Other proteins are plotted according to the weighted center of mass of their PPI scores with those references.</p>\n",
  "<p>For analyses involving only three reference proteins, set the first protein selection to <strong>NULL</strong> to generate a triangular representation.</p>\n",
  "<p><strong>Interactive parameters:</strong></p>\n",
  "<ul><li><strong>Minimum Entropy:</strong> Sets the lower bound of signal entropy for proteins included in the visualization.</li><li><strong>Maximum Entropy:</strong> Sets the upper bound of signal entropy.</li><li><strong>Max Value Min:</strong> Requires at least this PPI score with one of the selected reference proteins.</li><li><strong>Reference Colors:</strong> Choose a color for each reference vertex using the color selector or a HEX value.</li><li><strong>Repel Labels:</strong> Optionally separates displayed labels in 3D and adds connector lines. <strong>Maximum Labels</strong> controls the number of labels drawn; all qualifying points remain visible.</li></ul>\n",
  "<p>Together, these parameters help investigate distinct or shared interaction preferences and patterns of interaction dominance.</p>\n",
  "<h2>2. Complex Mode</h2>\n",
  "<h3>2.1 Network-Proximal Complexes</h3>\n",
  "<p>This tool extends interaction analysis to protein complexes, enabling identification and visualization of network relationships between complexes. Complex–complex interaction (CCI) scores are calculated for pairs of complexes using a permutation-based statistical framework, as described in the <strong>Methods</strong> section.</p>\n",
  "<p>For a selected complex, the tool ranks potential interacting complexes according to their CCI scores, highlighting associations identified through network proximity analysis.</p>\n",
  "<p><strong>Interactive parameter:</strong></p>\n",
  "<ul><li><strong>Top Hits to Label:</strong> Controls the number of highest-ranking interacting complexes labeled in the visualization.</li></ul>\n",
  "<h3>2.2 Feature Dominance</h3>\n",
  "<p>Complex Feature Dominance extends the protein-level analysis to protein complexes. Users can select three or four reference complexes and visualize relative interaction strengths based on CCI scores. Four references form a tetrahedron; selecting <strong>NULL</strong> for the first reference produces a triangle.</p>\n",
  "<p><strong>Interactive parameters:</strong></p>\n",
  "<ul><li><strong>Minimum Entropy:</strong> Defines the lower bound of signal entropy for complexes included in the visualization.</li><li><strong>Maximum Entropy:</strong> Defines the upper bound of signal entropy.</li><li><strong>Max Value Min:</strong> Requires at least this CCI score with one selected reference complex.</li><li><strong>Reference Colors:</strong> Customize each reference color with the selector or a HEX value.</li><li><strong>Repel Labels / Maximum Labels:</strong> Control 3D label separation and the maximum number of labels, without removing underlying data points.</li></ul>\n",
  "<p>This analysis facilitates identification of complexes with preferential or shared interaction patterns across the selected reference complexes.</p>
"
)


ui <- dashboardPage(
  
  dashboardHeader(
    title = tags$img(
      src = "CoBRAA_logo.jpg",
      height = "100px",width="235px",
      style = "padding: 0px 0px; margin-left: -15px;"
    )
  ),
  
  dashboardSidebar(
    tags$head(
      tags$style(HTML("
    .main-header .logo {
      height: 100px;
    }
    .main-header .navbar {
      height: 100px;
    }
    .content-wrapper, .main-sidebar {
      margin-top: 40px;
    }
    body, .wrapper {
      padding-top: 0px !important;
    }
    .content-wrapper > .content {
      padding-top: 15px;
    }
    a.shiny-download-link {
      color: white !important;
    }
    .btn-primary {
      color: white !important;
    }
    /* Indent h4 headers and download buttons in sidebar */
    .sidebar h4 {
      margin-left: 15px;
    }
    .sidebar .shiny-download-link {
      margin-left: 15px;
    }
  "))
    ),
    
    sidebarMenu(id = "main_tab",
      menuItem("Analysis", tabName = "analysis", icon = icon("chart-line"), selected = TRUE),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    conditionalPanel(
      condition = "input.main_tab == 'analysis' || !input.main_tab",
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
          "Feature Dominance" = "dominant"
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
        downloadButton("protein_pop_download", "Download Data (.tsv)",class = "btn-primary")  # ← ADD THIS LINE
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
        downloadButton("protein_recip_download", "Download Data (.tsv)",class = "btn-primary")  # ← ADD THIS
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
        
        tags$h5("Reference Colors"),
        fluidRow(
          column(6, colourpicker::colourInput("protein_dom_col1", "Reference 1", value = "#E15759", showColour = "both")),
          column(6, colourpicker::colourInput("protein_dom_col2", "Reference 2", value = "#59A14F", showColour = "both")),
          column(6, colourpicker::colourInput("protein_dom_col3", "Reference 3", value = "#4E79A7", showColour = "both")),
          column(6, colourpicker::colourInput("protein_dom_col4", "Reference 4", value = "#F28E2B", showColour = "both"))
        ),
        checkboxInput("protein_dom_repel", "Repel labels in 3D", value = TRUE),
        sliderInput("protein_dom_maxlabels", "Maximum Labels:", min = 0, max = 100, value = 25, step = 5),

        actionButton("protein_dom_plot", "Generate Plot", class = "btn-primary"),
        br(), br(),
        downloadButton("protein_dom_download", "Download Data",class = "btn-primary")
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
          "Feature Dominance" = "dominant"
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
        downloadButton("complex_net_download", "Download Data (.tsv)",class = "btn-primary")
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
        
        tags$h5("Reference Colors"),
        fluidRow(
          column(6, colourpicker::colourInput("complex_dom_col1", "Reference 1", value = "#E15759", showColour = "both")),
          column(6, colourpicker::colourInput("complex_dom_col2", "Reference 2", value = "#59A14F", showColour = "both")),
          column(6, colourpicker::colourInput("complex_dom_col3", "Reference 3", value = "#4E79A7", showColour = "both")),
          column(6, colourpicker::colourInput("complex_dom_col4", "Reference 4", value = "#F28E2B", showColour = "both"))
        ),
        checkboxInput("complex_dom_repel", "Repel labels in 3D", value = TRUE),
        sliderInput("complex_dom_maxlabels", "Maximum Labels:", min = 0, max = 100, value = 25, step = 5),

        actionButton("complex_dom_plot", "Generate Plot", class = "btn-primary"),
        br(), br(),
        downloadButton("complex_dom_download", "Download Data",class = "btn-primary")
      )
    )
    ) # End analysis-only sidebar controls
  ),
  
  dashboardBody(
    tags$head(tags$style(HTML("
      .cobbra-about {max-width: 1050px; margin: 0 auto; padding: 24px 32px; line-height: 1.7; font-size: 15px; background: white; border-radius: 5px;}
      .cobbra-about h1 {font-weight: 700; margin-top: 0; margin-bottom: 22px;}
      .cobbra-about h2 {border-bottom: 1px solid #e5e9ed; padding-bottom: 10px; margin-top: 34px;}
      .cobbra-about h3 {margin-top: 27px; font-size: 20px;}
      .cobbra-about li {margin-bottom: 9px;}
      .cobbra-about .lead {font-size: 16px;}
    "))),
    tabItems(
      tabItem(tabName = "analysis",
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
      ),
      tabItem(tabName = "about",
        div(class = "cobbra-about", HTML(ABOUT_HTML))
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
             "dominant" = "Dominance Plot (Tetrahedron/Triangle)")
    } else {
      switch(input$complex_submode,
             "network" = "Network Proximal Complexes",
             "dominant" = "Dominance Plot (Tetrahedron/Triangle)")
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
      )
    
    output$data_table <- renderTable(
      table_data %>% head(1000),
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
    
    # TRY to generate plot, CATCH errors gracefully ← NEW
    result <- tryCatch({
      
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
          head(download_data, 1000),
          striped = TRUE, bordered = TRUE, hover = TRUE
        )
        
        # Store for download
        protein_recip_data(download_data)
        
        setProgress(1, detail = "Complete!")
      })
      
      return("success")  # If we get here, it worked
      
    }, error = function(e) {
      # HANDLE ERROR - show message instead of crashing
      showModal(modalDialog(
        title = "⚠️ Incomplete Data",
        HTML(paste0(
          "<p><strong>", e$message, "</strong></p>",
          "<p>The app will continue to work for other features.</p>"
        )),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return("error")
    })
    
  })  # ← Close observeEvent
  
  
  
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
        col = c(input$protein_dom_col1, input$protein_dom_col2,
                input$protein_dom_col3, input$protein_dom_col4),
        repel_labels = isTRUE(input$protein_dom_repel),
        max_labels = input$protein_dom_maxlabels,
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
      coi_scores <- CCI_sparse[input$complex_net_coi, ]
      
      # Create data frame and REMOVE the complex of interest itself
      plot_data <- data.frame(
        complex = names(coi_scores),
        score = as.numeric(coi_scores),
        stringsAsFactors = FALSE
      ) %>%
        filter(complex != input$complex_net_coi) %>%  # ← ADD THIS LINE
        arrange(desc(score)) %>%  # Sort by score (highest first)
        mutate(
          rank = row_number(),  # Rank 1 = highest score
          log_score = log(score + 1)
        ) %>%
        arrange(score) %>%  # Re-sort for plotting (visual: low to high)
        mutate(index = row_number())  # Index for x-axis
      
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
      
      label_data <- subset(
        plot_data,
        !is.na(label) & label != "" &
          is.finite(index) & is.finite(log_score)
      )
      
      # Create plot
      p <- plot_ly(plot_data) %>%
        add_trace(
          x = ~index,
          y = ~log_score,
          type = "scattergl",
          mode = "markers",
          marker = list(
            color = "rgb(100,150,200)",
            size = 10
          ),
        #  text = ~label,
         # textposition = "top center",
        #  textfont = list(size = 9, color = "black"),
          hovertext = ~paste0(
            "<b>", complex, "</b><br>",
            "Rank: ", rank, " / ", n_complexes, "<br>",
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
        )
      
      # Add repelled labels
      if (nrow(label_data) > 0) {
        
        p <- p %>%
          plotly.repel::add_text_repel(
            data = label_data,
            x = ~index,
            y = ~log_score,
            text = ~label,
            
            # Prioritize complexes with higher CCI scores
            priority = ~log_score,
            
            # Label repulsion
            force = 3,
            force_pull = 0.8,
            box_padding = 0.7,
            point_padding = 0.3,
            direction = "both",
            
            # Optimization
            max_time_ms = 100,
            max_iter = 1000,
            max_overlaps = Inf,
            seed = 42,
            
            # Connecting lines
            segment = list(
              color = "rgba(100,100,100,0.5)",
              width = 1
            ),
            
            # Label appearance
            font = list(
              size = 10,
              color = "black"
            ),
            
            # Recalculate labels after zooming
            on = c("render", "zoom", "resize")
          )
      }
      
      
      # Plotly configuration
      p <- p %>%
        config(
          scrollZoom = TRUE,
          displayModeBar = TRUE,
          modeBarButtonsToRemove = c("lasso2d", "select2d")
        )
      
      output$main_plot <- renderPlotly(p)
      
      # Prepare download data
      download_data <- plot_data %>%
        arrange(desc(score)) %>%
        select(Rank = rank, Complex = complex, Score = score, Log_Score = log_score)
      
      # Show table (top 100)
      output$data_table <- renderTable(
        head(download_data, 1000),
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
        col = c(input$complex_dom_col1, input$complex_dom_col2,
                input$complex_dom_col3, input$complex_dom_col4),
        repel_labels = isTRUE(input$complex_dom_repel),
        max_labels = input$complex_dom_maxlabels,
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