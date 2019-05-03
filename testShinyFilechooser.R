library(shiny)
library(shinyFiles)


ui <- fluidPage(
  shinyFilesButton("Btn_GetFile", "Choose a file" ,
                   title = "Please select a file:", multiple = FALSE,
                   buttonType = "default", class = NULL),
  
  textOutput("txt_file")     
)


server <- function(input,output,session){
  
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
shinyApp(ui = ui, server = server)