library(shiny)
library(bslib)


ui <-page_sidebar( 
      # App title ----
      title = "Hello Shiny",
      #Sidebar panel for for inputs
      sidebar = sidebar(
        # Input: Slider for the number of bins  ---
        sliderInput(
          inputId = "bins",
          label = "number of bins:",
          min = 1,
          max = 50,
          value = 30
        ) 
      ),
      #output: histogram ---
      plotOutput( outputId = "distPlot")
      )

server <-function(input, output){
  
  output$distPlot <- renderPlot({
    x <- faithful$waiting
    
    bins <- seq(min(x), max(x), length.out = input$bins + 1 )
   
     hist(x, breaks = bins, col = "#007bc2", border = "white", 
         xlab = "waiting time to next eruption (in mins)",
         main = "Histogram of waiting times")
  })
}

#ui <- ...
#server <- ...
shinyApp(ui= ui, server = server)

runApp("my_app")
