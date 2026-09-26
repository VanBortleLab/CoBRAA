



ent=function(data)
{
  data=data/rowSums(data)
  eata=log2(data)
  fata=eata %>% as.data.frame() %>% dplyr::mutate(across(.cols = everything(), ~ ifelse(is.infinite(.x), 0, .x)))
  return(-rowSums(data*fata))
}


plotTetrahedron=function(data,
                         entropyrange=c(0,Inf),
                         maxvaluerange=c(0,Inf),
                         col=c("red","green","blue","orange"),
                         labels=NULL,
                         output_table=F,remove_background=T,cex=1,pch=16,
                         title="",push_text=1.2,
                         margin=5)
{
  #require(dplyr)
  # data=y %>% select("POLR3A","POLR3GL","POLR3C","POLR3G") %>% slice_sample(n=500)
  # col=c("red","green","blue","orange")
  # normalized_magnituderange=c(0.01,1)
  # data=hmm;remove_background = T;normalized_magnituderange = c(0.001,1);output_table = T
  
  #Creating a Tetrahedron and finding center of mass
  {  
    data[data<0]=0

    vertices <- matrix(c(   0,           0,  sqrt(2),  # Maps from ( 1,  1,  1)
                            1,           0,        0,  # Maps from ( 1, -1, -1)
                            -1/2,   sqrt(3)/2,        0,  # Maps from (-1,  1, -1)
                            -1/2,  -sqrt(3)/2,        0), # Maps from (-1, -1,  1) This new vertices is used so that the triangle plot remains in 2d xy plane
                       ncol = 3, byrow = TRUE)
    
  com_=data.frame(comx=rowSums(as.matrix(data)%*%diag(vertices[,1])),
                  comy=rowSums(as.matrix(data)%*%diag(vertices[,2])),
                  comz=rowSums(as.matrix(data)%*%diag(vertices[,3])))/rowSums(data)
  vertices11 = as.data.frame(vertices) %>% mutate(labels = colnames(data)) #for plotting vertex
  edges11=as.data.frame(vertices[c(1,2,3,1,4,2,3,4),]) #for plotting edges
  edge_color=rep("black",7)
  if(colnames(data)[1]=="NULL")
  {edge_color[c(1,3,4)]="transparent"
    col[1]="transparent"}
  final=cbind(data,com_)
  }##luka##
  
  #Finding entropy,magnitude,max value and color
  {  
    #Finding entropy and magnitude
    final=final %>% dplyr::mutate(entropy=ent(data))

    #Finding max value and assigning color on the basis of max value
    final=final %>% dplyr::mutate(max=apply(data,1,max)) %>% 
      dplyr::mutate(color=ifelse(.[[1]]==max,col[1],ifelse(.[[2]]==max,col[2],ifelse(.[[3]]==max,col[3],col[4]))))
  }##Luka##
  
  #Filtering on the basis of entropy, magnitude and max value
  {
    #Filtering on the basis of entropy
    final$color[!(final$entropy>=entropyrange[1] & 
                    final$entropy<=entropyrange[2])]="whitesmoke"
    
    #Filtering on the basis of maxvalue
    final$color[!(final$max>=maxvaluerange[1] & 
                    final$max<=maxvaluerange[2])]="whitesmoke"
    #Adding labels and colors
    if(is.null(labels)) 
    {print("No labels provided. Using rownames as labels");labels=rownames(final)}
    final=final %>% mutate(color=factor(color),labels=labels)  #factorizing color coz thats needed in plotly
  }##Luka##
  
  #Plotting
  {
    
    final_good= final %>% filter(color!="whitesmoke")
    if(remove_background==T) {final_bad=final %>% dplyr::slice()} else
    {   final_bad= final %>% filter(color=="whitesmoke")}
    
    p= plot_ly(
      data = vertices11,
      x = ~V1*1.5,y = ~V2*1.5,z = ~V3*1.5,  #multiplying by 1.5 to push away the vertices label
      type = 'scatter3d',mode = 'text',
      text = ~labels,textposition = 'top center',
      textfont = list(size = 20, color = col)
    ) %>%    ##Now adding background
      add_trace(
        data = final_bad ,
        x = ~comx,y = ~comy,z = ~comz,
        type = 'scatter3d',mode = 'markers+text',
        text = ~labels,textposition = 'top center',
        textfont = list(size = 10, color=~color,colors=col),
        marker = list(size = 10,color=~color,colors=col )
      )  %>%   ##Now adding frontground
      add_trace(
        data = final_good ,
        x = ~comx,y = ~comy,z = ~comz,
        type = 'scatter3d',mode = 'markers+text',
        text = ~labels,textposition = 'top center',
        textfont = list(size = 10, color=~color,colors=col),
        marker = list(size = 10,color=~color,colors=col )
      ) 
    
    # %>%   ##Now drawing the edges
    #   add_trace(
    #     data = edges11,
    #     x = ~V1,y = ~V2,z = ~V3,
    #     type = 'scatter3d',mode = 'lines+text',
    #     line = list(color = edge_color, width = 10),
    #     text=NULL
    #   )
    ##Now drawing edges one by one in 7 different layers(traces)
    for (i in 1:(nrow(edges11) - 1)) {
      p <- p %>% add_trace(
        x = edges11$V1[i:(i+1)],
        y = edges11$V2[i:(i+1)],
        z = edges11$V3[i:(i+1)],
        type = 'scatter3d',
        mode = 'lines+text',
        line = list(color = edge_color[i], width = 10),
        text=NULL
      )
    }

        p <- p %>% layout(scene = list(
      camera = list(eye = list(x = 0, y = 0, z = 2.5),up  = list(x = 1, y = 0, z = 0))))     #This line is added to change the viewpoint for new tetrahedorn plot which has a side in 2d xy plane. 

  }##Luka##
  
  #Returning table and plot
  if (output_table) {
    return(list(plot = p, table = final_good))  # `p` is the plotly object
  } else {
    return(p)
  }
  # if(output_table==T) return(final)
}##luka##





get_popularity_corrected_plot <- function(poi, label_top_n = 50) {
  
  # Get the row for protein of interest
  tu1 <- sort(PPI_sparse[poi, ],decreasing = T)
  
  # Get BioGRID info for this protein
  biogrid_row <- Biogrid_binary[poi, ]
  
  # Create data frame (tu1 is already sorted descending, so rank goes 1,2,3...)
  n_proteins <- length(tu1)
  plot_data <- data.frame(
    protein = names(tu1),
    rank = 1:n_proteins,  # Rank: 1 = highest score
    index = n_proteins:1,  # Index reversed for x-axis (visual left-to-right)
    score = as.numeric(tu1),
    in_biogrid = as.logical(biogrid_row[names(tu1)]),
    stringsAsFactors = FALSE
  )
  
  # Add category column
  plot_data$category <- ifelse(
    plot_data$score >= 2 & plot_data$in_biogrid, 
    "Popularity-corrected",
    ifelse(
      plot_data$score >= 2 & !plot_data$in_biogrid,
      "Predicted",
      "Unassociated"
    )
  )
  
  # Assign discrete colors
  plot_data$point_color <- ifelse(
    plot_data$score >= 2 & plot_data$in_biogrid, 
    "#487A9C",      # blue
    ifelse(
      plot_data$score >= 2 & !plot_data$in_biogrid,
      "#FB8451",    # orange
      "#E8E8E8"     # grey
    )
  )
  
  # Label only top N proteins
  n_proteins <- nrow(plot_data)
  plot_data$label <- ifelse(
    plot_data$index > (n_proteins - label_top_n),
    plot_data$protein,
    ""
  )
  
  # Calculate denser y-axis ticks (about 10 ticks)
  score_range <- range(plot_data$score, na.rm = TRUE)
  y_min <- floor(score_range[1])
  y_max <- ceiling(score_range[2])
  tick_interval <- max(1, ceiling((y_max - y_min) / 10))  # About 10 ticks
  y_ticks <- seq(y_min, y_max, by = tick_interval)
  
  # Create interactive Plotly plot
  p <- plot_ly(plot_data) %>%
    add_trace(
      x = ~index,
      y = ~score,
      type = "scattergl",
      mode = "markers+text",
      marker = list(
        color = plot_data$point_color,
        size = 10
      ),
      text = ~label,
      textposition = "top center",
      textfont = list(size = 9, color = "black"),
      hovertext = ~paste0(
        "<b>", protein, "</b><br>",
        "Category: ", category, "<br>",
        "Rank: ", rank, " / ", n_proteins, "<br>",  # ← Changed index to rank
        "Score: ", round(score, 4), "<br>",
        "In BioGRID: ", ifelse(in_biogrid, "Yes", "No")
      ),
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
    # Add horizontal line at y = 2
    add_segments(
      x = 1, 
      xend = n_proteins,
      y = 2, 
      yend = 2,
      line = list(color = "red", width = 2, dash = "dash"),
      showlegend = FALSE,
      hoverinfo = "skip"
    ) %>%
    layout(
      title = list(
        text = paste("PPI Scores for", poi, 
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
        title = "PPI Score",
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
      dragmode = "pan",
      annotations = list(
        list(
          x = 0.02, y = 0.98,
          xref = "paper", yref = "paper",
          text = paste0(
            "<b>Legend:</b><br>",
            "<span style='color:#FB8451'>●</span> Predicted<br>",
            "<span style='color:#487A9C'>●</span> Popularity-corrected<br>",
            "<span style='color:#E8E8E8'>●</span> Unassociated<br>",
            "<span style='color:red'>- - -</span> Threshold (score = 2)"
          ),
          showarrow = FALSE,
          xanchor = "left",
          yanchor = "top",
          bgcolor = "rgba(255,255,255,0.9)",
          bordercolor = "rgba(0,0,0,0.3)",
          borderwidth = 1,
          borderpad = 4,
          font = list(size = 11),
          align = "left"
        )
      )
    ) %>%
    config(
      scrollZoom = TRUE,
      displayModeBar = TRUE,
      modeBarButtonsToRemove = c("lasso2d", "select2d")
    )
  
  return(p)
}


get_info_from_reciprocal_plot <- function(prot, label_threshold = 0.5) {
  
  # Load protein index
  protein_index <- readRDS("data/tem_chunks/protein_index.rds")
  all_proteins <- protein_index$protein
  
  # Find protein index and chunk
  idx <- protein_index$index[protein_index$protein == prot]
  if (length(idx) == 0) stop("Protein not found: ", prot)
  
  chunk_size <- 100
  chunk_num <- ceiling(idx / chunk_size)
  row_in_chunk <- ((idx - 1) %% chunk_size) + 1
  
  # CHECK IF FILES EXIST BEFORE LOADING ← NEW
  z_file <- sprintf("data/tem_chunks/z/%d.Rdata", chunk_num)
  rz_file <- sprintf("data/tem_chunks/rz/%d.Rdata", chunk_num)
  
  if (!file.exists(z_file) || !file.exists(rz_file)) {
    stop(paste0(
      "INCOMPLETE DATA: Missing files for protein '", prot, "'.\n\n",
      "You have an incomplete version of the data folder.\n",
      "Please download the full version from:\n",
      "https://zenodo.org/records/21268626\n\n",  # ← REPLACE WITH YOUR ACTUAL ZENODO LINK
      "and replace the old data folder in your directory."
    ))
  }
  # ← END NEW CODE
  
  # Load chunks (rest of function stays the same)
  z_key <- as.character(chunk_num)
  if (!exists(z_key, envir = z_cache)) {
    z_cache[[z_key]] <- readRDS(z_file)
  }
  if (!exists(z_key, envir = rz_cache)) {
    rz_cache[[z_key]] <- readRDS(rz_file)
  }
  
  # ... rest of your existing code stays exactly the same ...
  # Extract data
  row_data <- z_cache[[z_key]][row_in_chunk, ]
  col_data <- rz_cache[[z_key]][, row_in_chunk]
  
  tub <- data.frame(
    protein = all_proteins,
    self = as.numeric(row_data),
    other = as.numeric(col_data),
    stringsAsFactors = FALSE
  )
  rownames(tub) <- tub$protein
  
  tub <- tub[tub$protein != prot, ]
  tub <- tub[!is.na(tub$self), ]
  tub <- tub %>% arrange(self)
  
  fit <- smooth.spline(tub$self, tub$other, df = 18000, spar = NULL)
  tub$predicted <- predict(fit, x = tub$self)$y
  
  hya <- tub %>% 
    arrange(desc(self)) %>% 
    slice(1:50) %>% 
    mutate(diff = abs(self - lead(self)))
  
  hya <- hya %>% 
    mutate(z_diff = (diff - mean(hya$diff[10:49], na.rm = TRUE)) / sd(hya$diff[10:49], na.rm = TRUE)) %>% 
    filter(z_diff > 50)
  
  tub2 <- tub[!rownames(tub) %in% rownames(hya), ]
  tub2 <- tub2[abs(tub2$other - tub2$predicted) < 0.01, ]
  
  fit2 <- smooth.spline(tub2$self, tub2$other, df = 1000, spar = NULL)
  tub$predicted <- predict(fit2, x = tub$self)$y
  tub <- tub %>% mutate(distances = other - predicted)
  
  ppi_row <- PPI_sparse[prot, ]
  biogrid_row <- Biogrid_binary[prot, ]
  
  tub$ppi_score <- as.numeric(ppi_row[tub$protein])
  tub$in_biogrid <- as.logical(biogrid_row[tub$protein])
  
  tub$ppi_score[is.na(tub$ppi_score)] <- 0
  tub$in_biogrid[is.na(tub$in_biogrid)] <- FALSE
  
  tub$category <- ifelse(
    tub$ppi_score >= 2 & tub$in_biogrid,
    "Popularity-corrected",
    ifelse(
      tub$ppi_score >= 2 & !tub$in_biogrid,
      "Predicted",
      "Unassociated"
    )
  )
  
  tub$point_color <- ifelse(
    tub$ppi_score >= 2 & tub$in_biogrid,
    "#487A9C",
    ifelse(
      tub$ppi_score >= 2 & !tub$in_biogrid,
      "#FB8451",
      "#E8E8E8"
    )
  )
  
  tub$label <- ifelse(
    tub$distances > label_threshold,
    tub$protein,
    ""
  )
  
  tub <- tub %>% arrange(self)
  
  return(tub)
}

get_reciprocal_plot <- function(recip_data, poi) {
  
  # Calculate very sparse y-axis ticks
  y_range <- range(c(recip_data$self, recip_data$other), na.rm = TRUE)
  y_min <- floor(y_range[1])
  y_max <- ceiling(y_range[2])
  y_ticks <- seq(y_min, y_max, by = max(1, ceiling((y_max - y_min) / 5)))
  
  # Calculate very sparse x-axis ticks
  x_range <- range(recip_data$self, na.rm = TRUE)
  x_min <- floor(x_range[1])
  x_max <- ceiling(x_range[2])
  x_ticks <- seq(x_min, x_max, by = max(1, ceiling((x_max - x_min) / 5)))
  
  # Create plot
  p <- plot_ly(recip_data) %>%
    add_trace(
      x = ~self,
      y = ~other,
      type = "scattergl",
      mode = "markers+text",
      marker = list(
        color = recip_data$point_color,
        size = 8
      ),
      text = ~label,
      textposition = "top center",
      textfont = list(size = 9, color = "black"),
      hovertext = ~paste0(
        "<b>", protein, "</b><br>",
        "Category: ", category, "<br>",
        "Z score: ", round(self, 4), "<br>",
        "Reciprocal Z score: ", round(other, 4), "<br>",
        "|Deviation|: ", round(distances, 4), "<br>",
        "PPI Score: ", round(ppi_score, 2)
      ),
      hoverinfo = "text",
      showlegend = FALSE
    ) %>%
    # Add spline fitted line
    add_trace(
      x = ~self,
      y = ~predicted,
      type = "scattergl",
      mode = "lines",
      line = list(color = "black", width = 3),
      hoverinfo = "none",
      showlegend = FALSE
    ) %>%
    layout(
      title = list(
        text = paste("Reciprocal Analysis:", poi),
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Z",
        showgrid = TRUE,
        gridcolor = "rgba(200,200,200,0.2)",
        showticklabels = TRUE,
        tickmode = "array",
        tickvals = x_ticks,
        tickformat = ".0f",
        zeroline = FALSE
      ),
      yaxis = list(
        title = "Reciprocal Z",
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
      dragmode = "pan",
      annotations = list(
        list(
          x = 0.02, y = 0.98,
          xref = "paper", yref = "paper",
          text = paste0(
            "<b>Legend:</b><br>",
            "<span style='color:#FB8451'>●</span> Predicted<br>",
            "<span style='color:#487A9C'>●</span> Popularity-corrected<br>",
            "<span style='color:#E8E8E8'>●</span> Unassociated"
          ),
          showarrow = FALSE,
          xanchor = "left",
          yanchor = "top",
          bgcolor = "rgba(255,255,255,0.9)",
          bordercolor = "rgba(0,0,0,0.3)",
          borderwidth = 1,
          borderpad = 4,
          font = list(size = 11),
          align = "left"
        )
      )
    ) %>%
    config(
      scrollZoom = TRUE,
      displayModeBar = TRUE,
      modeBarButtonsToRemove = c("lasso2d", "select2d")
    )
  
  return(p)
}


#This function worked before the error for downloading tem_chunks from zenodo was setup.
get_info_from_reciprocal_plot_backup <- function(prot, label_threshold = 0.5) {
  
  # Load protein index
  protein_index <- readRDS("data/tem_chunks/protein_index.rds")
  all_proteins <- protein_index$protein
  
  # Find protein index and chunk
  idx <- protein_index$index[protein_index$protein == prot]
  if (length(idx) == 0) stop("Protein not found: ", prot)
  
  chunk_size <- 100
  chunk_num <- ceiling(idx / chunk_size)
  row_in_chunk <- ((idx - 1) %% chunk_size) + 1
  
  # Load chunks
  z_key <- as.character(chunk_num)
  if (!exists(z_key, envir = z_cache)) {
    z_cache[[z_key]] <- readRDS(sprintf("data/tem_chunks/z/%d.Rdata", chunk_num))
  }
  if (!exists(z_key, envir = rz_cache)) {
    rz_cache[[z_key]] <- readRDS(sprintf("data/tem_chunks/rz/%d.Rdata", chunk_num))
  }
  
  # Extract data
  row_data <- z_cache[[z_key]][row_in_chunk, ]
  col_data <- rz_cache[[z_key]][, row_in_chunk]
  
  # Create dataframe
  tub <- data.frame(
    protein = all_proteins,
    self = as.numeric(row_data),
    other = as.numeric(col_data),
    stringsAsFactors = FALSE
  )
  rownames(tub) <- tub$protein
  
  # REMOVE THE POI ITSELF (the self-interaction)
  tub <- tub[tub$protein != prot, ]  # ← Remove by name
  
  # Remove NAs
  tub <- tub[!is.na(tub$self), ]
  
  # Arrange by self
  tub <- tub %>% arrange(self)
  
  # First spline fit
  fit <- smooth.spline(tub$self, tub$other, df = 18000, spar = NULL)
  tub$predicted <- predict(fit, x = tub$self)$y
  
  # Remove outliers on right side
  hya <- tub %>% 
    arrange(desc(self)) %>% 
    slice(1:50) %>% 
    mutate(diff = abs(self - lead(self)))
  
  hya <- hya %>% 
    mutate(z_diff = (diff - mean(hya$diff[10:49], na.rm = TRUE)) / sd(hya$diff[10:49], na.rm = TRUE)) %>% 
    filter(z_diff > 50)
  
  tub2 <- tub[!rownames(tub) %in% rownames(hya), ]
  tub2 <- tub2[abs(tub2$other - tub2$predicted) < 0.01, ]
  
  # Second refined spline fit
  fit2 <- smooth.spline(tub2$self, tub2$other, df = 1000, spar = NULL)
  tub$predicted <- predict(fit2, x = tub$self)$y
  tub <- tub %>% mutate(distances = other - predicted)
  
  
  # Get PPI scores and BioGRID info - SAFE indexing
  ppi_row <- PPI_sparse[prot, ]
  biogrid_row <- Biogrid_binary[prot, ]
  
  tub$ppi_score <- as.numeric(ppi_row[tub$protein])
  tub$in_biogrid <- as.logical(biogrid_row[tub$protein])
  
  # Replace NAs with defaults
  tub$ppi_score[is.na(tub$ppi_score)] <- 0
  tub$in_biogrid[is.na(tub$in_biogrid)] <- FALSE
  
  
  # Assign categories and colors
  tub$category <- ifelse(
    tub$ppi_score >= 2 & tub$in_biogrid,
    "Popularity-corrected",
    ifelse(
      tub$ppi_score >= 2 & !tub$in_biogrid,
      "Predicted",
      "Unassociated"
    )
  )
  
  tub$point_color <- ifelse(
    tub$ppi_score >= 2 & tub$in_biogrid,
    "#487A9C",
    ifelse(
      tub$ppi_score >= 2 & !tub$in_biogrid,
      "#FB8451",
      "#E8E8E8"
    )
  )
  
  # Label based on distance threshold
  tub$label <- ifelse(
    tub$distances > label_threshold,
    tub$protein,
    ""
  )
  
  # Arrange by self for plotting
  tub <- tub %>% arrange(self)
  
  return(tub)
}

