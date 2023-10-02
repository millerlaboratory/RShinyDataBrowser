library(shiny)
library(readr)
library(ggplot2)
library(scales)
library(dplyr)
library(DescTools)# added for Closest()
library(data.table)
library(DT)
library(shinythemes)
library(patchwork)




data <- read.delim("Data/tr_1000g.tsv")


data$copy_number <- as.numeric(data$copy_number)
data$size <- as.numeric(data$size)

data <- data.frame(data)

chr <- c(1:22, "X")
chr <- tibble(chr=paste0("chr", chr))


methyl_data <- read.delim("Data/chrX_methylationperIsland_subset.tsv")

cytoBandColor <- read.delim("Data/cytoband_color.tsv") %>%
  dplyr::filter(chr=="chrX")


merged_BWS <- read_delim("Data/hg38_1000g_BWS_all.tsv", 
                         delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

merged_PWA <- read_delim("Data/hg38_1000g_PWA_all.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

ui <- tagList(
  navbarPage(theme = shinytheme("cerulean"), "1000G ONT Data Explorer (Beta)",
             h4("**Nothing is final, adjustments ongoing**"),
  tabPanel("Tandem Repeat Copy Numbers (Beta)",
           h2("Tandem Repeat Copy Numbers"),
    sidebarPanel(
                 p("This is the Tandem Repeat Copy Number App. There is the Copy Number Distribution tab
                 and the Two Alleles tab.
                 Currently over 40 disease-causing tandem repeat 
                 loci can be viewed. Soon it will be over 500,000 simple repeat loci.
                  On the Copy Number Distribution tab, select the desired repeat according to chromosome 
                  and adjust the binwidth to your liking.
                   Click on a histogram bin to view the samples included.
                   The Two Alleles tab shows the copy numbers for both alleleles fro samples with 2.
                   If you navigate away from the Copy Number Distribution Tab and return, please re-select your locus of 
                   interest to re-configure the app.
                   Copy numbers are estimated with straglr (Chiu et al., 2021)."),
      
      selectInput("chr", h3("Select Chromosome:"), choices = unique(chr$chr)),
      selectInput("gene", h3("Select Known Tandem Repeat:"), choices = NULL),
      selectInput("repeat_motif", h3("Select Repeat Motif"), choices= NULL),
      sliderInput("binwidth", h3("Binwidth:"), min = 0.1, max = 25, value = 1, step = 0.2)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Copy Number Distribution", 
              plotOutput("hist", click = 'plot_click'),  # added plot_click
              DTOutput("selected_rows"), ),
        
        tabPanel("Two Alleles", 
                 plotOutput("twoAlleles"))
        )

    )
  ),
  tabPanel("Methylation Plots (Beta)",
           h2("Methylation Plots"),
      sidebarPanel(
                    h5("Panel 1: Whole X Chromosome Methylation"),
                    p("This is a viewer for the average fraction of reads methylated at a selection of informative CpG Islands on the X chromosome. 
                      If there is little data on the X chromosome plot, most of it unmethylated, the sample is XY. 
                      Improvements to labeling are ongoing"),
                    h5("Panel 2: Sample Autosomal DMRs"),
                    p("This viewer plots the average methylation of several DMRs at loci correlating to imprinting disorders.
                    Each haplotype is represented as a box at the average fraction of reads methylated spanning the DMR locus. 
                    With location relative to the chromosomal region."),
      selectInput("sample", h3("Select Sample:"), choices = unique(data$sample)
               )
             ),
      mainPanel(
        tabsetPanel(
          tabPanel("Whole X Chromosome",
                   plotOutput("methylPlot")
          ),
          tabPanel("Autosomal DMR sample",
                   plotOutput("BWSPlot"),
                   plotOutput("PWAPlot")
          ),
        ),
       
        
           ),
  ),
  #tabPanel("Known Autosomal DMR Browser", "Coming soon!"),
  #tabPanel("Structural Variant Tool (Coming soon)", "Coming soon!"),
  tabPanel("Acknowledgements",
           p("These data were made available by the 1000G ONT Consortium and the tools were developed in the Miller Lab at 
             The University of Washington."),
           p("Questions? Email Sophia Gibson at sophiabg@uw.edu"))
 )
)







server <- function(input, output, session) { 
  
  toListen <- reactive({list( input$chr # user updates the range slider
    , input$gene , input$repeat_motif   # user updates the number input
  )
  })
  
  observe({
    selected_chr <- input$chr
    available_genes <- data %>%
      filter(chr == selected_chr) %>%
      dplyr::select(gene) %>%
      distinct() %>%
      arrange(gene) %>%
      pull()
    updateSelectInput(session, "gene", choices = available_genes)
  })
  
  observe({
    selected_gene <- input$gene
    available_repeats <- data %>%
      filter(gene == selected_gene) %>%
      dplyr::select(repeat_unit) %>%
      distinct() %>%
      arrange(repeat_unit) %>%
      pull()
    updateSelectInput(session, "repeat_motif", choices = available_repeats)
  })
  
  
  
  
  observeEvent(toListen(), {
    
    filtered_data <- data %>%
      filter(chr == input$chr, gene == input$gene, repeat_unit== input$repeat_motif) %>%
      dplyr::arrange(copy_number)
    
    
    output$hist <- renderPlot(
      
      ggplot(filtered_data, aes(x=copy_number, stat="identity")) +
        geom_histogram(binwidth = input$binwidth, fill="#FF9999") +
        xlab("Copy number") +
        ylab("Sample count") +
        theme(axis.text = element_text(size=18),
              axis.title = element_text(size=18),
              panel.background = element_rect(fill = "white"),
              panel.border = element_rect(colour = "black", fill=NA),
              panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                              colour = "grey"))
      
    )
      
    
  })
  
  
  two <- reactive({
    data %>%
      filter(chr == input$chr, gene == input$gene, repeat_unit== input$repeat_motif) %>%
      filter(allele_count == 2)
  })
  
  output$twoAlleles <- renderPlot({
    ggplot(two(), aes(x=sample, y=copy_number, fill=allele)) +
      geom_bar(stat="identity", position=position_dodge()) +
      xlab("Sample")+
      ylab("Copy number")+
      scale_fill_manual(values=c("#CC6666", "#9999CC"))+
      theme(axis.text = element_text(size=18),
            axis.title = element_text(size=18),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            panel.background = element_rect(fill = "white"),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                            colour = "grey"))
  })
  
  
 output$selected_rows <- renderDT({
    if (is.null(input$plot_click$x)) return()
    
#find the closest bin center to show where the user clicked on the histogram
   bd2 <- ggplot_build(ggplot2::last_plot())$data[[1]]
    
   cBin <- which(bd2$x == Closest(bd2$x, input$plot_click$x))
    
    chr <- input$chr
    
    bin <- data.table(bin=1:nrow(bd2), chr=chr, start=bd2$xmin+0.01, end=bd2$xmax)
    data.table::setkey(bin, chr, start, end)
    
    filtered_data2 <- data %>%
      filter(chr == input$chr, gene == input$gene, repeat_unit== input$repeat_motif) %>%
      dplyr::arrange(copy_number) %>%
      tidyr::drop_na()
    
    cn <- data.table(chr=chr, start=filtered_data2$copy_number, 
                     end=filtered_data2$copy_number+.001, 
                     sample=filtered_data2$sample, gene=filtered_data2$gene)
    
    intersect <- foverlaps(cn, bin, by.x=c("chr", "start", "end")) %>%
      dplyr::select(i.start, sample, bin, gene) %>%
      dplyr::rename(copy_number = i.start)
    
    bin_data <- left_join(filtered_data2, intersect)
    
   table  <- bin_data %>% 
     filter(bin == cBin) %>% 
     dplyr::select(sample, gene, size, copy_number, repeat_unit, allele) %>%
     dplyr::rename(repeat_length=size, "repeat"=repeat_unit, copy_number=copy_number)
    table
  })
 
 
 
 
 
 
 filtered_BWS_data <- reactive({
   merged_BWS %>%
     filter(sample %in% input$sample) %>%
     dplyr::filter(DMR=="H19" | DMR =="KvDMR1")
 })
 
 
output$BWSPlot <- renderPlot({
  
  
  b <- ggplot(filtered_BWS_data(), aes(x=(start+stop)/2,y=mean_per_mod, fill=Haplotype)) +
    geom_point(shape=22,color="black", size=8)+
    scale_fill_manual(values=c("#56B4E9", "#D81B60"))+
    scale_color_manual(values=c("#56B4E9", "#D81B60"))+
    geom_text(data=filtered_BWS_data(), aes(x=(start+stop)/2, y=(mean_per_mod+4),  
                             label = ifelse(Haplotype == "hp2", DMR, "")), vjust = -0.5)+
    xlab("Position (bp)")+
    ylab("Average Fraction Methylated")+
    ylim(0, 100)+
    ggtitle("Beckwith-Wiedemann Syndrome DMRs")+
    coord_cartesian(xlim=c(0,2800000))+
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          strip.text = element_text(size=15))
  
  
  cytoBandColor <- read.delim("Data/cytoband_color.tsv") %>%
    dplyr::filter(label=="p15.5")
  
  
  b_rect <- ggplot(cytoBandColor) +
    geom_rect(aes(xmin = start, xmax = stop, ymin = 0, ymax = 0.25, fill=color), color = "black", size = 0.1, alpha = 0.5) +
    theme_void() +
    scale_fill_identity()+
    annotate("text", x = (cytoBandColor$start + cytoBandColor$stop) / 2, y = -.1,
             label = cytoBandColor$label, vjust = 0.2, hjust = 0, size = 2.5) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(10, 10, 30, 10)) 
  
  combined_plot_BWS <- (b / b_rect) +
    plot_layout(ncol = 1, heights = c(9, 2))
  
  combined_plot_BWS
  
  
  #ggplot(filtered_BWS_data()) +
    #geom_rect(aes(xmin = start, xmax = stop, ymin = mean_per_mod, ymax = mean_per_mod, fill=Haplotype, color=Haplotype), size = 10, alpha = 0.5) +
    #geom_errorbar(aes(ymin = mean_per_mod - stderr, ymax = mean_per_mod + stderr), width=0.2) +
    #scale_fill_manual(values=c("#56B4E9", "#D81B60"))+
    #scale_color_manual(values=c("#56B4E9", "#D81B60"))+
    #annotate("text", x = (filtered_BWS_data()$start + filtered_BWS_data()$stop) / 2, y = filtered_BWS_data()$mean_per_mod + 6,
             #label = filtered_BWS_data()$DMR) +
    #xlab("Position (bp)")+
    #ylab("Average Fraction Methylated")+
    #ylim(0, 100)+
    #theme(panel.background = element_blank(),
          #panel.border = element_rect(color = "black", fill = NA),
          #strip.text = element_text(size=15))
  
})
 
filtered_PWA_data <- reactive({
  merged_PWA %>%
    filter(sample %in% input$sample)
})
 
 
output$PWAPlot <- renderPlot({
  
  
  p <- ggplot(filtered_PWA_data(), aes(x=(start+stop)/2,y=mean_per_mod, fill=Haplotype)) +
    geom_point(shape=22,color="black", size=8)+
    scale_fill_manual(values=c("#56B4E9", "#D81B60"))+
    scale_color_manual(values=c("#56B4E9", "#D81B60"))+
    geom_text(data=filtered_PWA_data(), aes(x=(start+stop)/2, y=(mean_per_mod+4),  
                                            label = ifelse(Haplotype == "hp2", DMR, "")), vjust = -0.5)+
    xlab("Position (bp)")+
    ylab("Average Fraction Methylated")+
    ylim(0, 100)+
    ggtitle("Prader-Willi/Angleman Syndrome DMRs")+
    coord_cartesian(xlim=c(20500000, 25500000))+
    theme(panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          strip.text = element_text(size=15))
  
  
  cytoBandColor <- read.delim("Data/cytoband_color.tsv") %>%
    dplyr::filter(chr=="chr15") %>%
    dplyr::filter(label=="q11.2")
  
  
  p_rect <- ggplot(cytoBandColor) +
    geom_rect(aes(xmin = start, xmax = stop, ymin = 0, ymax = 0.25, fill=color), color = "black", size = 0.1, alpha = 0.5) +
    theme_void() +
    scale_fill_identity()+
    annotate("text", x = (cytoBandColor$start + cytoBandColor$stop) / 2, y = -.1,
             label = cytoBandColor$label, vjust = 0.2, hjust = 0, size = 2.5) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(10, 10, 30, 10)) 
  
  combined_plot_PWA <- (p / p_rect) +
    plot_layout(ncol = 1, heights = c(9, 2))
  
  combined_plot_PWA
  
  
})
 
 


 filtered_methyl_data <- reactive({
   methyl_data %>%
     filter(sample %in% input$sample)
 })
 
 output$methylPlot <- renderPlot({
   
   rect <- ggplot(cytoBandColor) +
     geom_rect(aes(xmin = start, xmax = stop, ymin = 0, ymax = 0.25, fill=color), color = "black", size = 0.1, alpha = 0.5) +
     theme_void() +
     scale_fill_identity()+
     annotate("text", x = (cytoBandColor$start + cytoBandColor$stop) / 2, y = -1,
              label = cytoBandColor$label, vjust = 0.2, hjust = 0, size = 2.5, angle = 90) +
     theme(axis.title.x = element_blank(),
           axis.text.x = element_blank(),
           axis.ticks.x = element_blank(),
           plot.margin = margin(10, 10, 30, 10)) 
   
   
   plot <- ggplot(filtered_methyl_data(), aes(x=start,y=mean_per_mod, fill=Haplotype)) +
     geom_errorbar(aes(ymin = mean_per_mod - stderr, ymax = mean_per_mod + stderr), width=0.2) +
     geom_point(shape=21,color="black", size=3) +
     scale_fill_manual(values=c("#56B4E9", "#D81B60"))+
     xlab("CpG Island Start Position (bp)")+
     ylab("Average Fraction Methylated")+
     ggtitle("Whole X chromosome CpG island methylation")+
     ylim(0, 100)+
     xlim(0,156040895)+
     theme(panel.background = element_blank(),
           panel.border = element_rect(color = "black", fill = NA),
           strip.text = element_text(size=15))
   
   # Arrange the two plots using patchwork library
   combined_plot <- (plot / rect) +
     plot_layout(ncol = 1, heights = c(9, 2))
   # Print the combined plot
   combined_plot
   
  

 })
 
 
 
 
}

shinyApp(ui = ui, server = server)