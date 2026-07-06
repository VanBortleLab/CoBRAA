# PPI_CCI
<img width="1031" height="892" alt="Screen Shot 2026-07-06 at 1 05 33 PM" src="https://github.com/user-attachments/assets/d44e01aa-66b3-4fd1-98fb-5b26af749163" />



> **Interactive analysis of Protein-Protein Interactions (PPI) and Complex-Complex Interactions (CCI) using a popularity-bias corrected network**

---

## 🧬 Overview

This Shiny application provides comprehensive tools for analyzing protein and complex interactions, correcting for popularity bias in protein-protein interaction networks. The app includes multiple analysis modes with interactive visualizations.

### Key Features

- **Popularity-Corrected PPI Analysis** - Identify true interactors while accounting for popularity bias
- **Reciprocal Analysis** - Analyse reciprocal plot for your protein of interest
- **Dominant Plot (Tetrahedron)** - 3D visualization of interaction specificity amongst 4(or 3) proteins (or complexes) of interest
- **Interactive Plotly Visualizations** - Zoom, pan, and explore your data
- **Data Export** - Download results as TSV files

---

## 🚀 Quick Start

### Prerequisites

- R (version ≥ 4.0)
- RStudio (recommended)

### Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/VanBortleLab/PPI_CCI_analysis.git
   cd PPI_CCI_analysis
   ```

2. **Download complete data files:**

   ⚠️ Important: This repository contains only a subset of data files needed for reciprocal plots. For full functionality of reciprocal plots, download the complete tem_chunks data:

   📦 Download: tem_chunks_full.zip from Zenodo
   📂 Extract to: data/tem_chunks/

3. **Install required R packages:**

   ```r
   install.packages(c("shiny","shinydashboard","shinyWidgets","plotly","dplyr","Matrix"))
   ```

4. **Run the app:**
   Open app.R and run the app
   ```r
   shiny::runApp(".")
   ```
