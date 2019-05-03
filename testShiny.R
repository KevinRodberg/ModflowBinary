library(shiny)
library(shinyFiles)
# Define UI ----
ui <- fluidPage(
  titlePanel("Basic widgets"),
  
  fluidRow(
    column(3,h3("Buttons"),actionButton("action", "Action"),br(),br(),
      submitButton("Submit")),
    
    column(3,h3("Single checkbox"),checkboxInput("checkbox", "Choice A", value = TRUE)),
    
    column(3,checkboxGroupInput("checkGroup",h3("Checkbox group"),
          choices = list("Choice 1" = 1,"Choice 2" = 2,"Choice 3" = 3),
        selected = 1)),
    
    column(3,
           dateInput("date",
                     h3("Date input"),
                     value = "2014-01-01"))
  ),
  
  fluidRow(
    column(3,
           dateRangeInput("dates", h3("Date range"))),
    
    column(3,
           shinyFilesButton("Btn_GetFile", "Choose a file" ,
                            title = "Please select a file:", multiple = FALSE,
                            buttonType = "default", class = NULL),
           textOutput("txt_file")),
    column(3,
      h3("Help text"),
      helpText(
        "Note: help text isn't a true widget,",
        "but it provides an easy way to add text to",
        "accompany other widgets."
      )
    ),
    
    column(3,
           numericInput("num",
                        h3("Numeric input"),
                        value = 1))
  ),
  
  fluidRow(
    column(3,
           radioButtons(
             "radio",
             h3("Radio buttons"),
             choices = list(
               "Choice 1" = 1,
               "Choice 2" = 2,
               "Choice 3" = 3
             ),
             selected = 1
           )),
    
    column(3,
           selectInput(
             "select",
             h3("Select box"),
             choices = list(
               "Choice 1" = 1,
               "Choice 2" = 2,
               "Choice 3" = 3
             ),
             selected = 1
           )),
    
    column(
      3,
      sliderInput(
        "slider1",
        h3("Sliders"),
        min = 0,
        max = 100,
        value = 50
      ),
      sliderInput(
        "slider2",
        "",
        min = 0,
        max = 100,
        value = c(25, 75)
      )
    ),
    
    column(3,
           textInput("text", h3("Text input"),
                     value = "Enter text..."))
  )
  
)

# Define server logic ----
server <- function(input, output,session) {
  options(shiny.maxRequestSize=3000*1024^2)
  volumes = getVolumes()
  observe({  
    shinyFileChoose(input, "Btn_GetFile", roots = volumes, session = session)
    
    if(!is.null(input$Btn_GetFile)){
      # browser()
      file_selected<-parseFilePaths(volumes, input$Btn_GetFile)
      output$txt_file <- renderText(as.character(file_selected$datapath))
    }
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)