library(shiny)
library(ggvis)

  ui <- fluidPage(
    
    titlePanel("Genewise Summary Data"),
      
    fluidRow(
      column(4,
             selectInput("infilter",
                         "Filter:",
                         c('none', 'phyloCSF', 'ncs', 'all'))),
      column(6,
             selectInput("ethnicity",
                         "Ancestry Group:",
                         c('All', 'Africa', 'Americas', 'ASJ', 'East Asia', 'Finland', 'Europe', 'Other')))
      ),
    
    tableOutput("data"),
    mainPanel(
      tabsetPanel(
        id = 'dataset',
        tabPanel("genes", DT::dataTableOutput("mytable1")),
        tabPanel("snps", DT::dataTableOutput("mytable2")),
        tabPanel('plot', ggvisOutput("plot2"))
      )
    ),
    fluidRow(
      column(9, DT::dataTableOutput('x1')),
      column(3, verbatimTextOutput('x4')),
      
      column(10,tags$em(paste0(
        "Note: PhyloCSF filters out variants that don't exhibit a pattern of conservation typical of protein-coding exons",
        " and those with reading frames that are likely offset.",
        " NCS filters splice variants that fall in a non-canonical splice site.",
        " All applies both filters."
      )))
    )
  )



# Define server logic required to draw a histogram
server <- function(input, output) {
    library('data.table')
    library('tidyverse')
    library('DT')
    library('ggvis')
  
  totalsnpdf<-fread('joined_snpdf')
  totaldf<-fread('totaldf')
  totaldf<-totaldf%>%ungroup%>%mutate(rowid=seq(n()))
  totalsnpdf<-totalsnpdf%>%ungroup%>%mutate(rowid=seq(n()))
  output$mytable1 <- DT::renderDataTable(
    DT::datatable({totaldf%>%filter(filter %in% input$infilter)%>%select(chr,gene_name,mips.ind,acmg.ind,variants,frameshift_variants,stopgain_variants,splice_variants,nmd_variants,n_heterozygotes,n_homozygotes)},  
                  rownames=FALSE, 
                  colnames = c('gene'='gene_name', 'MIPS indicator'='mips.ind', 'ACMG indicator'='acmg.ind', 'frameshift variants'='frameshift_variants',
                               'stopgain variants'='stopgain_variants','splice DA variants'='splice_variants', 'NMD variants'='nmd_variants',
                               'n heterozygotes'='n_heterozygotes',
                               'n homozygotes'='n_homozygotes'),
                  #extensions='Buttons',
                  options = list(
                    #dom = 'Bfrtip',
                    #buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
                    pageLength=100,
                  order = list(3, 'desc'),
                  lengthMenu = list(c(100, 500, 1000, -1), list('100', '500', '1000','All')),
                  columnDefs = list(list(visible=TRUE, targets=2))))%>%formatStyle('MIPS indicator',
                                                                                    target = 'row',
                                                                                    color = styleEqual(c(0, 1), c('black', 'red')),
                                                                                    fontWeight = styleEqual(c(0, 1), c('normal','bold')))
  )
  
  output$mytable2 <- DT::renderDataTable({
    s = input$mytable1_rows_selected
    gene = totaldf[totaldf$filter==input$infilter,]$gene_name[s]
    DT::datatable({totalsnpdf%>%filter(filter %in% input$infilter & gene_name %in% gene)%>%ungroup%>%select(chr, pos, type, ref, alt, gene_name,het,homozygotes,gnomad.indicator,
                                                                                                            gnomad_all_maf, gnomad_afr_maf, gnomad_amr_maf,gnomad_asj_maf,gnomad_eas_maf,gnomad_fin_maf,gnomad_nfe_maf,gnomad_oth_maf)}, 
                  rownames=FALSE, 
                  colnames = c('gene'='gene_name', 
                               'heterozygotes'='het'),
                            options(pageLength=100,
                            order = list(3, 'desc'),
                            lengthMenu = list(c(100, 500, 1000, -1), list('100', '500', '1000','All'))))
  })
  
  output$x4 = renderPrint({
    s = input$mytable1_rows_selected
    gene = totaldf[totaldf$filter==input$infilter,]$gene_name[s]
    if (length(s)) {
      cat('Genes: ')
      cat(gene, sep = ', ')
      cat('\n\n')
      cat('Filter: ')
      cat(input$infilter)
    }
  })
  
  gene_tooltip <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.null(x$rowid)) return(NULL)
    
    genedf <- totalsnpdf[totalsnpdf$rowid == x$rowid, ]
    
    paste0("<b>", genedf$chr, ":", genedf$pos, genedf$ref, '-', genedf$alt,"</b><br>",
           genedf$type,"</b><br>",
           'MAFs:',"</b><br>",
           "hunt: ", round(genedf$raf,4),"</b><br>",
           "gnomad:</b><br>", "all: ", round(genedf$gnomad_all_maf,4),"</b><br>",
           "afr: ", round(genedf$gnomad_afr_maf,4),"</b><br>",
           "asj: ", round(genedf$gnomad_asj_maf,4),"</b><br>",
           "eas: ", round(genedf$gnomad_eas_maf,4),"</b><br>",
           "fin: ", round(genedf$gnomad_fin_maf,4),"</b><br>",
           "nfe: ", round(genedf$gnomad_nfe_maf,4),"</b><br>",
           "oth: ", round(genedf$gnomad_oth_maf,4)
           #"$", format(movie$BoxOffice, big.mark = ",", scientific = FALSE)
    )
  }
  
  vis <- reactive({
    # Lables for axes
    s = input$mytable1_rows_selected
    gene = totaldf[totaldf$filter==input$infilter,]$gene_name[s]
    n<-nrow ( totalsnpdf%>% filter(filter %in% input$infilter))
    eth <- input$ethnicity
    #'All', 'Africa', 'Americas', 'ASJ', 'East Asia', 'Finland', 'Europe', 'Other'
    totalsnpdf %>% ungroup %>%mutate(plotgroup=case_when(eth=='All' ~ gnomad_all_maf,
                                             eth=='Africa' ~ gnomad_afr_maf,
                                             eth=='Americas' ~ gnomad_amr_maf,
                                             eth=='ASJ' ~ gnomad_asj_maf,
                                             eth=='East Asia' ~ gnomad_eas_maf,
                                             eth=='Finland' ~ gnomad_fin_maf,
                                             eth=='Europe' ~ gnomad_nfe_maf,
                                             eth=='Oth' ~ gnomad_oth_maf,
                                             TRUE ~ 0)) %>% filter(filter %in% input$infilter & !is.na(plotgroup)) %>% mutate(mips.ind=as.factor(mips.ind)
                                                                                                                              #gene.ind=as.factor(ifelse(gene_name %in% gene, 1, 0))
                                                                                                                              ) %>%
      ggvis(x = ~raf, y = ~plotgroup) %>%
      layer_points(size := 50, size.hover := 200,
                   fillOpacity := 0.4, 
                   fillOpacity.hover := 0.5,
                   key := ~rowid) %>%
      add_tooltip(gene_tooltip, "hover") %>%
      layer_paths(x=~raf,y=~raf,stroke:='lightgrey') %>%
      add_axis("x", title = 'hunt MAF') %>%
      add_axis("y", title = 'gnomad MAF') %>%
            set_options(width = 800, height = 500)
  })
  
  vis %>% bind_shiny("plot2")
  
  # output$plot2<-renderPlot({
  # ggplot(totaldf%>%filter(filter %in% input$infilter),aes(variants))+geom_density()+facet_wrap(~chr)},height = 400,width = 600)
}
  

# Run the application 
shinyApp(ui = ui, server = server)

