library(shiny)
library(ggplot2)
library(Seurat)
library(bslib)
library(ggtree)
library(ggnewscale)
library(ape)
library(phytools)
library(ggtree)
library(ggtreeExtra)
library(DT)
library(cowplot)
library(dplyr)

load("app.RData")

ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "journal", version = 5,  
                          primary = "#2A4765", secondary = "#F1EFE4",
                          success = "#73A790", info = "#A52D27", 
                          warning = "#D7B17C", danger = "#EABAB9"),
  
  br(),
  # App title and logo
  titlePanel(
    fluidRow(
      # Logo
      column(2, 
             img(height = 120, width = 120, src = "vibes.jpg")),
      # Title
      column(10, align = "left",
             div(style = "display: flex; align-items: center; height: 100%;",
                 HTML(
                   "<div style='font-size: 50px; '>
                   <strong>Discovering novel secondary metabolites produced by <i>Streptococcus mutans</i>
                   via genome mining approaches</strong><div>
                   <div style='font-size: 40px; '>
                   CompMicroLab
                   <div>"
                )))
    )), br(),
  
  #Tab: General information
  tabsetPanel(
    tabPanel("Project Overview",
             fluidRow(
               column(12, br(),
                      HTML(
                        "<h3>Project Overview</h3><br>
                        This thesis project focuses on discovering novel secondary metabolites produced by
                        <i>Streptococcus mutans</i> via genome mining approaches.<br>
                        <strong>Background</strong><br>
                        <i>Streptococcus mutans</i> is a gram-positive, facultatively anaerobic bacterium, commonly found in the human oral microbiome.
                        It has been implicated as one of the main contributors to the development of dental caries,as well as other extraoral health issues.
                        <i>S. mutans</i> produces secondary metabolites to interact with its environment and harbors bisynthetic gene clusters (BGCs),
                        some of which are responsible for producing bacteriocins, known as mutacins. These compounds have been under-studied in <i> S. mutans</i>
                        and represent a promising reservoir of novel antimicrobial metabolites.
                        <br>
                        <strong>Aim/Methodology</strong><br>
                        This project aims to provide a large-scale understanding of the biosyntheic potential of <i> S. mutans</i> by:
                        <br>
                        (i) identifying BGCs within all publicly available <i> S. mutans</i> genomes using different genomes detection methods (antiSMASH and GECCO).
                        <br>
                        (ii) partitioning all found <i>S. mutans</i> BGCs as well as thousands of validated BGCs (MIBiG) into gene cluster families (GCFs).
                        <br>
                        (iii) identifying novel GCFs with potential antimicrobial activities.
                        <br>
                        "
                      ),
                      br(), br()
               )
             )
    ),
    
    #Tab: UMAP
    tabPanel("UMAP",
             # Sidebar layout
             fluidRow(
               # Input choices
               column(2, br(),
                      selectInput("clustering", label = "Cluster by:",
                                  choices = c("Seurat cluster"="seurat_clusters",
                                              "Biosynthetic class"="class",
                                              "BGC detection method"="GCF_method", 
                                              "Cluster length"="cluster_length")),
                                  br(),
                                  # Information about method
                                  br(),
                                  htmlOutput("methodInfo")
               ),
               
               # Outputs
               column(8, align = "center", br(),
                      # UMAP
                      plotOutput("umap", width = "750px", height="550px")
               ),
               
               column(2, br(),
                      # Frequency table
                      tableOutput("frequency")
               )
             ), br(),
             
             fluidRow(
               column(2), # empty column
               # table
               column(10,
                      dataTableOutput("table")
               ),
             )
    ),
    
    # Tab: Phylogeny
    tabPanel("Phylogeny",
      fluidRow(
        column(2, br(),
          checkboxGroupInput("phylogenyChoices", "Panels to show:",
            c("Number of BGCs per genome", "Continent", "Isolation source", "MIBiG + antiSMASH/GECCO GCFs")),
          br(),
          tags$style(HTML("
            #submit_btn {
              background-color: #2A4765; /* Custom color */
              color: white;
              border: none;
              padding: 10px 20px;
              text-align: center;
              text-decoration: none;
              display: inline-block;
              font-size: 16px;
              margin: 4px 2px;
              cursor: pointer;
              border-radius: 8px;
            }
            #submit_btn:hover {
              background-color: #1f3450; /* Slightly darker shade for hover effect */
            }
          ")),
          textInput("GCF_phylogeny", "GCF to display:"),
          br(),
          textOutput("GCF_entered")
        ),
        column(7, br(),
               plotOutput("phylogeny", width = "750px", height="550px")
          
        ),
        column(3, br(),
               plotOutput("legends", width = "321px", height = "550px")
        )
      )
    )
  )
)

server <- function(input, output) {
  
  # UMAP plot output
  output$umap <- renderPlot({
    
    # setting conditional variables
    highlight <- NULL
    color <- NULL
    legend <- NULL
    title <- NULL
    
    # conditions for clustering by different methods
    if(input$clustering == "cluster_length"){
      highlight <- which(smutans[[input$clustering]] < 20000)
      color <- "#2A4765"
      legend <- NoLegend()
      title <- "Clusters with a length < 20 Kbp"
    }else if(input$clustering == "seurat_clusters"){
      title <- "Seurat clusters"
    } else if(input$clustering == "GCF_method"){
      title <- "BGC detection method"
    } else if(input$clustering == "class"){
      title <- "Biosynthetic class"
    }
    
    # UMAP plot
    DimPlot(object = smutans, label = FALSE, group.by = input$clustering,
            cells.highlight = highlight, cols.highlight = color, pt.size = 2,
            sizes.highlight = 2) + 
      legend +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      ggtitle(title)
  })
  
  # frequency table
  output$frequency <- renderTable({
    if(input$clustering != "cluster_length") {
      # Table for clustering other than cluster_length
      table_data <- table(smutans[[input$clustering]])
    } else {
      # Table for cluster_length, with modified row and column names
      table_data <- table(smutans[[input$clustering]] < 20000)
      table_data <- as.data.frame(table_data) # Convert to data frame
      colnames(table_data) <- c("cluster_length", "Freq")
      table_data$cluster_length <- c("Cluster Length < 20000 Kbp", "Cluster Length => 20000 Kbp")
    }
    return(table_data)
  })
  
  # conditional text information about the clustering method
  output$methodInfo <- renderUI({
    if(input$clustering == "seurat_clusters"){
      "Seurat uses a graph-based clustering approach to group similar GCFs together.
      The resolution was set to 0.5, resulting in 9 clusters"
    } else if(input$clustering == "cluster_length"){
      "GCFs coloured by cluster length either below or equal to/above 20 Kbp.
      This cutoff was chosen for practical reasons to ensure optimal molecular cloning."
    } else if(input$clustering == "GCF_method"){
      "If any MIBiG BGC is present, the GCF was colored as 'MIBiG'.
      Otherwise it was classified as 'GECCO', 'antiSMASH', or 'mixed'."
    } else if(input$clustering == "class"){
      "GCFs coloured by predicted biosynthetic class of the representative genome"
    }
  })
  
  # table
  output$table <- renderDataTable({
    datatable(
      smutans@meta.data[,c( "GCF","seurat_clusters", "class", "GCF_method", 
                            "cluster_length", "#ofBGCs", "genomes with BGC (%)")],
      rownames = FALSE,
      colnames = c("GCF", "Seurat cluster", "Class", "BGC detection method",
                   "Cluster length", "# of BGCs", "BGC presence (%)"),
      options = list(
        autoWidth = TRUE,
        columnDefs = list(list(width = '100px', targets = "_all"))),
      filter = "top"
    )
  })
  
  #legends
  output$legends <- renderPlot({
    phylogeny_legend <- ggplot() +
      theme_void()
    
    if ("Number of BGCs per genome" %in% input$phylogenyChoices) {
      phylogeny_legend <- phylogeny_legend + 
        annotation_custom(BGC_legend, xmin = 0, xmax = 0.4, ymin = 0.95, ymax = 1)
    }
    
    
    if ("Continent" %in% input$phylogenyChoices) {
      phylogeny_legend <- phylogeny_legend + 
        annotation_custom(continent_legend, xmin = 0, xmax = 0.125, ymin = 0.52, ymax = 1)
    }
    
    if ("Isolation source" %in% input$phylogenyChoices) {
      phylogeny_legend <- phylogeny_legend + 
        annotation_custom(source_legend, xmin = 0, xmax = 1.4, ymin = 0.715, ymax = 1)
    }
    
    if ("MIBiG + antiSMASH/GECCO GCFs" %in% input$phylogenyChoices) {
      phylogeny_legend <- phylogeny_legend + 
        annotation_custom(mixed_MIBiG_legend, xmin = 0, xmax = 0.51, ymin = -0.1, ymax = 1)
    }
    
    print(phylogeny_legend)
  })
  
  # phylogeny output
  output$phylogeny <- renderPlot({
    phylogeny_plot <- ggtree(rerooted_outgroup_droptip) +
      geom_treescale(x = 0.005, y = 0, fontsize = 3) +
      theme(legend.position = "none",
            plot.title = element_text(size = 24, hjust = 0.5)) +
      labs(title = expression(paste(italic("Streptococcus mutans"), " phylogeny")))

    
    if ("Number of BGCs per genome" %in% input$phylogenyChoices) {
      phylogeny_plot <- phylogeny_plot + 
        geom_fruit(data = BGCS_per_genome, geom = geom_col, 
                   mapping = aes(x = GECCO, y = genome), fill = "#73A790", 
                   color = "#73A790") +
        geom_fruit(data = BGCS_per_genome, geom = geom_col, 
                   mapping = aes(x = antiSMASH, y = genome), fill = "#EABAB9", 
                   color = "#EABAB9")
    }
    
    
    if ("Continent" %in% input$phylogenyChoices) {
      phylogeny_plot <- phylogeny_plot + 
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = continent), 
                   width = 0.0005) + 
        scale_fill_viridis_d(option = "D", name="Continent", na.value = "gray") + new_scale_fill() +
        geom_fruit(data = metaphylogeny, geom = geom_tile,  # Add empty geom_fruit layer
                   mapping = aes(y = genome), fill = "white", width = 0.0001)
    }
    
    if ("Isolation source" %in% input$phylogenyChoices) {
      phylogeny_plot <- phylogeny_plot + 
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = `isolation source`), 
                   width = 0.0005) +
        scale_fill_viridis_d(option = "A", name="Isolation source", na.value = "gray") + new_scale_fill() +
        geom_fruit(data = metaphylogeny, geom = geom_tile,  # Add empty geom_fruit layer
                   mapping = aes(y = genome), fill = "white", width = 0.0001)
    }
    
    if ("MIBiG + antiSMASH/GECCO GCFs" %in% input$phylogenyChoices) {
      phylogeny_plot <- phylogeny_plot + 
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = GCF0000070), width = 0.0002, pwidth = 0.05) +
        scale_fill_manual(values = "darkblue", na.value = "white") + new_scale_fill() +
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = GCF0000073), width = 0.0002, pwidth = 0.05) +
        scale_fill_manual(values = "darkred", na.value = "white") + new_scale_fill() +
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = GCF0000075), width = 0.0002, pwidth = 0.05) +
        scale_fill_manual(values = "darkgreen", na.value = "white") + new_scale_fill() +
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = GCF0000082), width = 0.0002, pwidth = 0.05) +
        scale_fill_manual(values = "darkorange", na.value = "white") + new_scale_fill() +
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = GCF0000092), width = 0.0002, pwidth = 0.05) +
        scale_fill_manual(values = "darkorchid", na.value = "white") + new_scale_fill() +
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = GCF0000100), width = 0.0002, pwidth = 0.05) +
        scale_fill_manual(values = "gold3", na.value = "white") + new_scale_fill() +
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = GCF0000923), width = 0.0002, pwidth = 0.05) +
        scale_fill_manual(values = "darkturquoise", na.value = "white") + new_scale_fill()
    }
    
    if(input$GCF_phylogeny %in% smutans$GCF){
      metaphylogeny$GCF_phylogeny <- NA
        
      output$GCF_entered <- renderText(paste(""))
        
      for(i in 1:length(clusters$gcf_id)){
        if(clusters$gcf_id[i] == input$GCF_phylogeny){
          ind <- which(metaphylogeny$genome == clusters$genome[i])
          metaphylogeny$GCF_phylogeny[ind] <- input$GCF_phylogeny
        }
      }
        
      phylogeny_plot <- phylogeny_plot +
        geom_fruit(data = metaphylogeny, geom = geom_tile,  # Add empty geom_fruit layer
                   mapping = aes(y = genome), fill = "white", width = 0.0001) +
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = GCF_phylogeny),
                   width = 0.0005, pwidth = 0.05) +
        scale_fill_manual(values = "black", na.value = "white") + new_scale_fill()
    } 
    
    else if(input$GCF_phylogeny == ""){
      output$GCF_entered <- renderText(paste(" "))
    }
      
    else {
      output$GCF_entered <- renderText(paste("Not a valid GCF"))
    }
    
    print(phylogeny_plot)
  })

}

shinyApp(ui = ui, server = server)
