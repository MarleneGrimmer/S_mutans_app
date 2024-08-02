library(shiny)
library(DT)
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
library(plotly)
library(htmlwidgets)
library(scales)

load("app.RData")

ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "journal", version = 5,        # setting theme and colours 
                          primary = "#2A4765", secondary = "#F1EFE4",
                          success = "#73A790", info = "#A52D27", 
                          warning = "#D7B17C", danger = "#EABAB9"),
  
  br(),
  # App title and logo
  titlePanel(
    fluidRow(
      # Logo
      column(2, align = "center", 
             img(height = 160, width = 160, src = "vibes.jpg")),  # inserting lab logo
      # Title
      column(10, align = "left",
             div(style = "display: flex; align-items: center; height: 100%;",
                 HTML(
                   "<div style='font-size: 45px; '>
                   <strong>Discovering novel secondary metabolites produced by <i>Streptococcus mutans</i>
                   via genome mining approaches</strong><div>
                   <div style='font-size: 35px; '>
                   CompMicroLab
                   <div>"
                )))
    )), br(),
  
  #Tab: General information
  tabsetPanel(
    tabPanel("Project Overview",
             fluidRow(
               column(12, br(),
                      HTML( # project overview text
                        "<h3>Project Overview</h3><br>
                        This thesis project focuses on discovering novel secondary metabolites produced by
                        <i>Streptococcus mutans</i> via genome mining approaches.<br>
                        <strong>Background</strong><br>
                        <i>Streptococcus mutans</i> is a gram-positive, facultatively anaerobic bacterium, commonly found in the human oral microbiome.
                        It has been implicated as one of the main contributors to the development of dental caries, as well as other extraoral health issues.
                        <i>S. mutans</i> produces secondary metabolites to interact with its environment and harbors bisynthetic gene clusters (BGCs),
                        some of which are responsible for producing bacteriocins, known as mutacins. These compounds have been under-studied in <i> S. mutans</i>
                        and represent a promising reservoir of novel antimicrobial metabolites.
                        <br>
                        <strong>Aim/Methodology</strong><br>
                        This project aims to provide a large-scale understanding of the biosynthetic potential of <i> S. mutans</i> by:
                        <br>
                        (i) identifying BGCs within all publicly available <i> S. mutans</i> genomes using different genomes detection methods (antiSMASH and GECCO).
                        <br>
                        (ii) partitioning all found <i>S. mutans</i> BGCs as well as thousands of validated BGCs (MIBiG database) into gene cluster families (GCFs).
                        <br>
                        (iii) identifying novel GCFs with potential antimicrobial activities.
                        <br>
                        "
                      ), br())),
             fluidRow(
               column(12, align = "center",
                      plotlyOutput("UMAP_3D", width = "750px")
               )
             )
    ),
    
    #Tab: UMAP
    tabPanel("UMAP",
             # Sidebar layout
             fluidRow(
               # Input choices
               column(2, br(),
                      # drop-down selection for colouring the GCFs
                      selectInput("clustering", label = "Cluster by:",
                                  choices = c("Seurat cluster"="seurat_clusters",
                                              "Biosynthetic class"="class",
                                              "BGC detection method"="GCF_method", 
                                              "Cluster length"="cluster_length")),
                                  br(),
                                  # Information about method dependent on selection input
                                  br(),
                                  htmlOutput("methodInfo")
               ),
               
               # Outputs
               column(8, align = "center", br(),
                      # UMAP
                      plotOutput("umap", width = "750px", height="550px")
               ),
               
               column(2, br(),
                      # frequency table (how many GCFs per group)
                      tableOutput("frequency")
               )
             ), br(),
             
             fluidRow(
               column(2), # empty column under selection input
               # table underneath the plot listing all GCFs
               column(10,
                      DTOutput("table")
               ),
             )
    ),
    
    # Tab: Phylogeny
    tabPanel("Phylogeny",
      fluidRow(
        column(2, br(), # column with input selections
          checkboxGroupInput("phylogenyChoices", "Panels to show:", # checkbox selection to display panels on phylogeny
            c("Number of BGCs per genome", "Continent", "Isolation source", "MIBiG + antiSMASH/GECCO GCFs")),
          br(),
          selectInput("GCF_select", "GCF to display:", c("",unique(clusters$gcf_id)))
        ),
        column(7, br(),
               plotOutput("phylogeny", width = "750px", height="550px") # phylogeny plot
          
        ),
        column(3, br(),
               plotOutput("legends", width = "321px", height = "550px") # display legends for selected panels 
        )
      )
    )
  )
)

server <- function(input, output) {
  # 3D rotating UMAP
  output$UMAP_3D <- renderPlotly({
    plot_ly(df_3D, x = ~umap1, y = ~umap2, z = ~umap3, color = ~cell,
            colors = hue_pal()(length(levels(df_3D$cell))), marker = list(size = 3)) %>% 
      add_markers() %>%
      layout(scene = list(xaxis = list(title = "umap 1", showticklabels=FALSE), 
                          yaxis = list(title = "umap 2", showticklabels=FALSE),
                          zaxis = list(title = "umap 3", showticklabels=FALSE),
                          camera = list(eye = list(x = 1.25, y = 1.25, z = 1.25),
                                        center = list(x = 0, y = 0, z = 0)
           ))) %>%
      onRender("
      function(el, x){
  var id = el.getAttribute('id');
  var gd = document.getElementById(id);
  Plotly.update(id).then(attach);
  function attach() {
    var cnt = 0;
    
    function run() {
      rotate('scene', Math.PI / 180);
      requestAnimationFrame(run);
    } 
    run();
    
    function rotate(id, angle) {
      var eye0 = gd.layout[id].camera.eye
      var rtz = xyz2rtz(eye0);
      rtz.t += angle;
      
      var eye1 = rtz2xyz(rtz);
      Plotly.relayout(gd, id + '.camera.eye', eye1)
    }
    
    function xyz2rtz(xyz) {
      return {
        r: Math.sqrt(xyz.x * xyz.x + xyz.y * xyz.y),
        t: Math.atan2(xyz.y, xyz.x),
        z: xyz.z
      };
    }
    
    function rtz2xyz(rtz) {
      return {
        x: rtz.r * Math.cos(rtz.t),
        y: rtz.r * Math.sin(rtz.t),
        z: rtz.z
      };
    }
  };
}
    ")
  })
  
  # UMAP plot output
  output$umap <- renderPlot({
    
    # setting conditional variables
    highlight <- NULL
    color <- NULL
    legend <- NULL
    title <- NULL
    
    # conditions for clustering by different methods
    if(input$clustering == "cluster_length"){ # selected cluster_length 
      highlight <- which(smutans[[input$clustering]] < 20000) # highlight all GCFs with length < 20000
      color <- "#2A4765"                                      # highlight colour
      legend <- NoLegend()                                    # remove legend 
      title <- "Clusters with a length < 20 Kbp"              # set title
    } else if(input$clustering == "seurat_clusters"){ # selected seurat_clusters
      title <- "Seurat clusters"                              # set title
    } else if(input$clustering == "GCF_method"){ # selected GCF_method
      title <- "BGC detection method"                         # set title
    } else if(input$clustering == "class"){ # selected biosynthetic class
      title <- "Biosynthetic class"                           # set title
    }
    
    # UMAP plot displaying clustering of GCFs
    DimPlot(object = smutans, label = FALSE, group.by = input$clustering, # grouped by selected input
            cells.highlight = highlight, cols.highlight = color, pt.size = 2, # highlight and color set conditionally
            sizes.highlight = 2) + 
      legend +  # legend set conditionally 
      theme(plot.title = element_text(hjust = 0.5), # title in the middle
            axis.text.x=element_blank(), # x and y axis labels removed 
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
      table_data <- table(smutans[[input$clustering]] < 20000) # frequency of GCFs smaller/longer than 20000 
      table_data <- as.data.frame(table_data) # Convert to data frame
      colnames(table_data) <- c("cluster_length", "Freq")
      table_data$cluster_length <- c("Cluster Length < 20000 Kbp", "Cluster Length => 20000 Kbp")
    }
    return(table_data)
  })
  
  # conditional text information about the clustering method selected
  output$methodInfo <- renderUI({
    if(input$clustering == "seurat_clusters"){
      "Seurat uses a graph-based clustering approach to group similar GCFs together.
      The resolution was set to 0.5, resulting in 9 clusters"
    } else if(input$clustering == "cluster_length"){
      "GCFs with a cluster length below 20 Kbp are highlighted.
      This cutoff was chosen for practical reasons to ensure optimal molecular cloning."
    } else if(input$clustering == "GCF_method"){
      "If any MIBiG BGC is present, the GCF was colored as 'MIBiG'.
      Otherwise it was classified as 'GECCO', 'antiSMASH', or 'mixed'."
    } else if(input$clustering == "class"){
      "GCFs coloured by predicted biosynthetic class of the representative genome"
    }
  })
  
  # GCF table displaying GCFs, assigned cluster, biosynthetic class, detection method, cluster length, number of BGCs in GCF, and presence of GCF in phylogeny 
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
  
  # legends to display alongside phylogeny based on boxes checked
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
  
  # phylogeny output graph
  output$phylogeny <- renderPlot({
    # default phylogeny without added metadata
    phylogeny_plot <- ggtree(rerooted_outgroup_droptip) + 
      geom_treescale(x = 0.005, y = 0, fontsize = 3) +
      theme(legend.position = "none",
            plot.title = element_text(size = 24, hjust = 0.5)) +
      labs(title = expression(paste(italic("Streptococcus mutans"), " phylogeny")))

    # adding panels based on boxes checked 
    if ("Number of BGCs per genome" %in% input$phylogenyChoices) {
      phylogeny_plot <- phylogeny_plot + 
        geom_fruit(data = BGCS_per_genome, geom = geom_col, 
                   mapping = aes(x = GECCO, y = genome), fill = "#7CAE00", 
                   color = "#7CAE00") +
        geom_fruit(data = BGCS_per_genome, geom = geom_col, 
                   mapping = aes(x = antiSMASH, y = genome), fill = "#F8766D", 
                   color = "#F8766D")
    }
    
    
    if ("Continent" %in% input$phylogenyChoices) {
      phylogeny_plot <- phylogeny_plot + 
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = continent), 
                   width = 0.0005) + 
        scale_fill_viridis_d(option = "D", name="Continent", na.value = "gray") + new_scale_fill() +
        geom_fruit(data = metaphylogeny, geom = geom_tile,  # Add empty geom_fruit layer for spacing
                   mapping = aes(y = genome), fill = "white", width = 0.0001)
    }
    
    if ("Isolation source" %in% input$phylogenyChoices) {
      phylogeny_plot <- phylogeny_plot + 
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = `isolation source`), 
                   width = 0.0005) +
        scale_fill_viridis_d(option = "A", name="Isolation source", na.value = "gray") + new_scale_fill() +
        geom_fruit(data = metaphylogeny, geom = geom_tile,  # Add empty geom_fruit layer for spacing
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
    
    # displaying presence of BGCs in entered GCF
    if(input$GCF_select != ""){ # valid GCF entered
      metaphylogeny$GCF_select <- NA # empty column
        
      for(i in 1:length(clusters$gcf_id)){ # going through all BGCs and ther GCF IDs
        if(clusters$gcf_id[i] == input$GCF_select){ # if the current BGC's GCF matches the input GCF
          ind <- which(metaphylogeny$genome == clusters$genome[i]) # which genome is the BGC in
          metaphylogeny$GCF_select[ind] <- input$GCF_select # add the GCF to the empty column at the right genome
        }
      }
        
      phylogeny_plot <- phylogeny_plot + # add the GCF presence layer to the end of the phylogeny graph
        geom_fruit(data = metaphylogeny, geom = geom_tile,  # Add empty geom_fruit layer for spacing
                   mapping = aes(y = genome), fill = "white", width = 0.0001) +
        geom_fruit(data = metaphylogeny, geom = geom_tile,
                   mapping = aes(y = genome, fill = GCF_phylogeny),
                   width = 0.0005, pwidth = 0.05) +
        scale_fill_manual(values = "black", na.value = "white") + new_scale_fill()
    } 
    
    print(phylogeny_plot)
  })

}

shinyApp(ui = ui, server = server)
