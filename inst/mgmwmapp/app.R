library(shiny)
library(ucimudata)
library(classimu)

################
# CONSTANTS
################
const.RENDER_PLOT_WIDTH = 1000
const.RENDER_PLOT_HEIGHT = 600
const.RENDER_PLOT_RES = 100 # default is 72

const.FIGURE_PLOT_HEIGHT = "600px"
const.FIGURE_PLOT_HEIGHT_REDUCED = "400px"
const.FIGURE_PLOT_HEIGHT_LOGO = "100px"

const.nb_of_digits = 7

# convert degrees-per-second to radians-per-second
const.degps_2_radps = 1/360 * 2*pi

# constant default frequency for custom data
const.DEFAULT_FREQ = 1 # [Hz]


################
# FUNCTIONS for COSINE INTEGRAL calculations
################
cos_function <- function(t){
  cos(t)/t
}

Ci <- function(x){
  -integrate(f = cos_function, lower = x, upper = 2e3, subdivisions=10000)$value
}

VCi <- Vectorize(Ci, c("x"))

sigma2_T <- function(T, f0, B){
  2*B*B/pi * ( log(2) - ( (sin(pi*f0*T))^3 ) / (2*(pi*f0*T)^2) * ( sin(pi*f0*T)+4*pi*f0*T*cos(pi*f0*T) ) + VCi(2*pi*f0*T) - VCi(4*pi*f0*T) )
}




# loading the four internal datasets
data("KVH1750imuAcc")

smac_url <- a("https://smac-group.github.io/gui4gmwm/", href="https://smac-group.github.io/gui4gmwm/")
smac_url_description <- "gui4gmwm Overview:"

# increses file limit from default-5MB to 100MB
options(shiny.maxRequestSize=100*1024^2) 



ui <- shinyUI(fluidPage(
  
  shinyjs::useShinyjs(),
  
  tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: red}")),
  tags$style(HTML(".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {background: green}")),
  tags$style(type='text/css', '#summ {background-color: rgba(0,0,200,0.02); color: black; width: 500px; font-size: 14px;}'), 
  
  
  title = "mGMWM GUI",
  tabsetPanel(id = "tabs",
              tabPanel("Model Data", plotOutput(outputId = "plot", height = const.FIGURE_PLOT_HEIGHT)),
              tabPanel("Selected Sensor", plotOutput(outputId = "plot2", height = const.FIGURE_PLOT_HEIGHT)),
              tabPanel("Summary", verbatimTextOutput(outputId = "summ", placeholder = FALSE)),
              tabPanel("Help",
                       # fluidPage("cluster"),
                       h4("Help Tab" ),
                       br(),
                       # actionButton("subClust", label = "Create Subcluster"),
                       # 
                       uiOutput(outputId = "tabhelpurl"),
                       br(),br(),
                       fluidRow(
                         column(5,
                                plotOutput(outputId = "tabhelpplotlogo_pennstate", height = const.FIGURE_PLOT_HEIGHT_LOGO)
                         ),
                         column(5,
                                plotOutput(outputId = "tabhelpplotlogo_epfl", height = const.FIGURE_PLOT_HEIGHT_LOGO)
                         )
                       )
              )
  ),
  
  hr(),
  
  fluidRow(
    column(4,
           radioButtons("data_input_choice", "Select data input:", choices = c("from library" = "library", "custom" = "custom")),
           
           conditionalPanel(
             condition = "input.data_input_choice == 'library'",
             
             selectInput("imu_obj", "Select IMU file:",
                         c("KVH 1750 Acc" = "KVH1750imuAcc"),
                         selected = 1),
             
             selectInput("sensors", "Select sensor", c("1"="1","2"="2", selected = 1))
           ),
           
           conditionalPanel(
             condition = "input.data_input_choice == 'custom'",
             
             fileInput("user_defined_txt_file", "Select INPUT file (max 100MB):",
                       accept = c(
                         "text/txt",
                         "text/comma-separated-values,text/plain",
                         ".txt",
                         ".imu",
                         placeholder = "No file selected")
             ),
             sliderInput("user_defined_txt_file_column", "Select column number:",
                         min=1, max=6, value=1),
             numericInput("user_defined_txt_frequency", label = "Set frequency of dataset", value = const.DEFAULT_FREQ), # frequency defined by the user
             textInput("user_defined_units", "Define units of active dataset", "rad/s")
             
           ),
           
           
           actionButton("fit1", label = "Plot WV"),
           
           br(),
           
           uiOutput("choose_columns")
    )
  )
)
)

server <- function(input, output, session) {
  # library or custom dataset
  v <- reactiveValues(plot = FALSE,
                      fit = FALSE,
                      mgmwm = NULL,
                      form = NULL,
                      freq = 100,
                      first_gmwm = NULL,
                      # how to transorm that in a vector?
                      n = NULL,
                      sensor_name = NULL,
                      sensor_column = NULL,
                      overlap_datasheet = FALSE,
                      y_label_with_dataunits = NA)
  
  # PUSHING ON BUTTON "Plot WV"
  observeEvent(input$fit1, {
    
    withProgress(message = 'Calculating empirical WV...', value = 0, {
      
      v$plot = TRUE
      v$fit = FALSE
      v$overlap_datasheet = "datasheet" %in% input$option_plot
      
      if ("library" %in% input$data_input_choice){ #using library data
        my_data = get(input$imu_obj)
        if(input$sensors == "Accel. X"){
          Xt = input$imu_obj[[1]]
        } else if (input$sensors == "Accel. Y"){
          Xt = input$imu_obj[[2]]
        } else if (input$sensors == "Accel. Z"){
          Xt = input$imu_obj[[3]]
        } 
        
        v$sensor_name = input$imu_obj
        v$sensor_column = input$sensors
        v$freq = attr(my_data, 'freq')
        v$custom_data = FALSE
        if (input$sensors == "Gyro. X" || input$sensors == "Gyro. Y" || input$sensors == "Gyro. Z"){
          v$y_label_with_dataunits = expression(paste("Wavelet Variance ", nu, " [", deg/s, "]"))
        } else if (input$sensors == "Accel. X" || input$sensors == "Accel. Y" || input$sensors == "Accel. Z"){
          v$y_label_with_dataunits = expression(paste("Wavelet Variance ", nu, " [", g, "]"))
        }
      }
      v$n = length(Xt)
      # Transform that in a vector?
      v$form = Xt$variance
      
      updateNavbarPage(session, "tabs", selected = "Selected Sensor")
    })
  })
  
  
}

# calc a specific VW and plot it in the tab "Selected Sensor"
output$plot2 <- renderPlot({
  
  if (v$fit || v$plot){
    # for the real data
    a = v$form
    freq_a = v$freq
    a$scales = a$scales#/freq_a
    duration_a = v$n/(freq_a*60*60)
    
    if (v$plot){ # should i plot just the real data?
      if (v$custom_data){ # is it custom data from a txt file?
        title = paste("Haar Wavelet Variance of TXT-FILE:\n", v$custom_data_name, " (column # ", input$user_defined_txt_file_column, " / ", v$custom_data_tot_colums,
                      ") - Filesize: ", round(v$custom_data_size/1024/1024,2), " [MB] - Duration: ", round(duration_a,1), "(h) @", freq_a, "(Hz)", sep = "")
      }else{ # it is NOT custom data
        title = paste("Haar Wavelet Variance of DATASET:\n", input$imu_obj, " (", input$sensors,
                      ") - Duration: ", round(duration_a,1), "(h) @", freq_a, "(Hz)", sep = "")
      }
      
      if ("datasheet" %in% input$option_plot){
        plot_wv_and_datasheet(a,
                              v$datasheet_noise_model,
                              # v$actual_datasheet_BI_parameter,
                              expression(paste("Scale ", tau, " [s]")),
                              v$y_label_with_dataunits,
                              prov_title = title)
      } else {
        plot(a,
             axis.x.label = expression(paste("Scale ", tau, " [s]")),
             axis.y.label = v$y_label_with_dataunits,
             title = title,
             CI = T, #"ci" %in% input$option_plot,
             title.size = 22, 
             axis.label.size = 20, 
             axis.tick.size = 17, 
             legend.title.size = 19, 
             legend.text.size = 19) + theme(legend.position = c(0, 0))
      }
      
    }else{ # when doing the "gmwm modeling" plot
      
      if (v$custom_data){ # is it custom data from a txt file?
        title = paste("Haar Wavelet Variance of TXT-FILE-DATA:\n", v$custom_data_name, " (column number ", input$user_defined_txt_file_column,
                      ") - Filesize: ", round(v$custom_data_size/1024/1024,2), " [MB] - Duration: ", round(duration_a,1), "(h) @", freq_a, "(Hz)", sep = "")
      }else{ # it is NOT custom data
        title = paste("Haar Wavelet Variance of DATASET:\n", input$imu_obj, " (", input$sensors,
                      ") - Duration: ", round(duration_a,1), "(h) @", freq_a, "(Hz)", sep = "")
      }
      
      if ("datasheet" %in% input$option_plot & !"process_decomp" %in% input$option_plot){
        plot_gmwm_and_datasheet(object = a, 
                                datasheet = v$datasheet_noise_model, 
                                # v$actual_datasheet_BI_parameter,
                                axis.x.label = expression(paste("Scale ", tau, " [s]")),
                                v$y_label_with_dataunits,
                                prov_title = title)
      }else{
        plot(a,
             axis.x.label = expression(paste("Scale ", tau, " [s]")),
             axis.y.label = v$y_label_with_dataunits,
             process.decomp = "process_decomp" %in% input$option_plot,
             CI = T, #"ci" %in% input$option_plot,
             title = title,
             title.size = 22, 
             axis.label.size = 20, 
             axis.tick.size = 17, 
             legend.title.size = 19, 
             legend.text.size = 19) #+ theme(legend.position = c(0.1, 0.1))
        
      }
    }
  }else{
    plot(NA)
  }
}, height = const.RENDER_PLOT_HEIGHT, width = const.RENDER_PLOT_WIDTH, res = const.RENDER_PLOT_RES)

shinyApp(ui, server)