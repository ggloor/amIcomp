#
# This is a Shiny web application. 
#

library(shiny)
library(ALDEx2)
library(aIc)

ui <- fluidPage(

  # App title ----
  titlePanel("Compositional test"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(


     
      fileInput("upload", NULL, accept=c(".csv", ".tsv", ".txt"), buttonLabel = "Upload .tsv ..."), 
# hard coded to 2 lines
      numericInput("z.rem", "max 0 proportion", value = 0.95, min = 0, max=1, step = 0.05),
      tableOutput("files") ,
      numericInput("group_1_size", "group 1 size", value=7, min=3, step=1),
      numericInput("group_2_size", "group 2 size", value=7, min=3, step=1),
#      tableOutput('groups'),
      
      # Input: Selector for test ----
      selectInput("test", "Compositional test:",
                  c("distance dominance" = "dom",
                    "perturbation invariance" = "pert",
                    "scale invariance" = "scale",
                    "correlation coherence" = "cohere",
                    "data singularity" = "sing")),

      # Input: Selector for normalization ----
      selectInput("norm", "Data normalization:",
                  c("proportion" = "prop",
                    "centred log ratio" = "clr",
                    "intraquartile log ratio" = "iqlr",
                    "low V, high abund log ratio" = "lvha",
                    "edgeR TMM" = "TMM",
                    "edgeR TMMwsp" = "TMMwsp",
                    "DESeq RLE" = "RLE",
                    "none" = "none")),

      # Input: Selector for zero treatment ----
      selectInput("zero", "Zero replacement:",
                  c("prior" = "prior",
                    "Geometric Baysian" = "GBM",
                    "Count zero mult" = "CZM",
                   "none" = "none")),

      # Input: Checkbox for whether outliers should be included ----
      checkboxInput("log", "Take log of transform (prop, TMM, RLE)", F)

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Formatted text for caption ----
      h4(textOutput("caption")),

      # Output: Plot of the requested variable against mpg ----
      plotOutput("testPlot"),
      
      tableOutput("groups"),
      tableOutput("head")
    )
  )
)

# Define server logic to plot various variables and load data ----
server <- function(input, output) {

  
   # Generate a plot of the requested variable against mpg ----
  # and only exclude outliers if requested
  output$testPlot <- renderPlot({

 up.data <- reactive({
    req(input$upload)
    
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
      csv = vroom::vroom(input$upload$datapath, delim = ","),
      tsv = read.table(input$upload$datapath,header=T, row.names=1, sep = "\t"),
      validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
  output$head <- renderTable({
    head(up.data(), 2)
  })
  
#  # set groups 
    group <- c(rep('A', input$group_1_size), rep("B", input$group_2_size))
    zero.method = input$zero
    
    if(input$zero=="none") {zero.method = NULL }
     
# perturbation 
    if(input$test=='pert'){
      x <- aIc.perturb(up.data(), norm.method=input$norm, zero.method=zero.method, zero.remove=input$z.rem,  log=input$log, group=group)
      plot(x$plot, main=x$main, xlab=x$xlab, ylab=x$ylab)
      abline(v=0, lty=2, lwd=3, col='red')
      if(x$is.perturb == 'Yes'){
        output$caption <- renderText({paste('The data are approximately perturbation invariant with transform ', input$norm,". The maximum observed absolute perturbation is: ", x$ol, "%.", sep="")})
      } else if(x$is.perturb == 'No') {
        output$caption <- renderText({paste('The data are not perturbation invariant with transform ', input$norm,'. The maximum observed absolute perturbation is: ', x$ol, '%. Please try the clr transform on this dataset.', sep="")})
      }

# dominance
    } else if (input$test=='dom'){
      x <- aIc.dominant(up.data(), norm.method=input$norm,zero.method=zero.method, zero.remove=input$z.rem,  log=input$log, group=group)
      plot(x$plot, main=x$main, xlab=x$xlab, ylab=x$ylab)
      abline(v=0, lty=2, lwd=3, col='red')
      if(x$is.dominant == 'Yes'){
        output$caption <- renderText({paste('The data are distance dominant with transform ', input$norm,". The proportion of non-dominant distances in the sub-compositon is: ", round(x$ol,2), "%.", sep="")})
      } else if(x$is.dominant == 'No') {
        output$caption <- renderText({paste('The data are not distance dominant with transform ', input$norm,'. The proportion of non-dominant distances in the sub-compositon is: ', round(x$ol,2), '%. Please try the clr transform on this dataset.', sep="")})
      }

# scale
    } else if (input$test=='scale'){
      x <- aIc.scale(up.data(), norm.method=input$norm, zero.method=zero.method, zero.remove=input$z.rem,  log=input$log, group=group)
      plot(x$plot, main=x$main, xlab=x$xlab, ylab=x$ylab)
      abline(v=0, lty=2, lwd=3, col='red')
      if(x$is.scale == 'Yes'){
        output$caption <- renderText({paste('The data are scale invariant with transform ', input$norm,". The proportion of non-scale invariand distances in the sub-compositon is: ", round(x$ol,2), "%.", sep="")})
      } else if(x$is.scale == 'No') {
        output$caption <- renderText({paste('The data are not distance dominant with transform ', input$norm,'. The proportion of non-scale consistent distances is: ', round(x$ol,2), '%. Please try the clr transform on this dataset.', sep="")})
      }

# coherence      
    } else if (input$test=='cohere'){
      x <- aIc.coherent(up.data(), norm.method=input$norm, zero.method=zero.method, zero.remove=input$z.rem,  log=input$log, group=group)
      plot(x$plot[,1], x$plot[,2], main=x$main, xlab=x$xlab, ylab=x$ylab)
      abline(0,1, lty=2, lwd=3, col='red')
      if(x$is.coherent == 'Yes' & input$norm != 'none'){
        output$caption <- renderText({paste('The data are sub-compositionally coherent with ', input$norm,". Network analysis and correlation inference may be OK.", sep="")})
      } else if(x$is.coherent == 'No' & input$norm != 'none') {
        output$caption <- renderText({paste('The data are not sub-compositionally coherent with transform ', input$norm,'. You should use some form of compositional association test and should check out the propr R package. In addition, when doing dimension reduction you need to be aware of compositional effects; see Greenacre and Aitchison 2002 "Biplots of compositional data" JRSA 51:375', sep="")})
      } else if(x$is.coherent == 'Yes' & input$norm == 'none') {
        output$caption <- renderText({paste('The data must be transformed in some way prior to use. You should use some form of compositional association test and should check out the propr R package. In addition, when doing dimension reduction you need to be aware of compositional effects; see Greenacre and Aitchison 2002 "Biplots of compositional data" JRSA 51:375', sep="")})
      }

    # singularity      
    } else if (input$test=='sing'){
      x <- aIc.singular(up.data(), norm.method=input$norm, zero.method=zero.method, zero.remove=input$z.rem,  log=input$log, group=group)
      if(x$is.singular == 'Yes'){
       output$caption <- renderText({paste('The data are singular with transform ', input$norm,'. When doing dimension reduction you need to be aware of compositional effects; see Greenacre and Aitchison 2002 "Biplots of compositional data" JRSA 51:375. You likey will be best served using a compositional analysis approach.', sep="")})
      } else if (x$is.singular == 'No'){
       output$caption <- renderText({paste('The data are not singular with transform ', input$norm,'. Go about your daily business.', sep="")})
      }
    }
  })
  
 
}

# Create Shiny app ----
shinyApp(ui, server)