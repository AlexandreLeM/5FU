## app.R ##
library(shiny)
library(shinydashboard)
library(shiny)
library(mapbayr)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plyr)


adapt_poso <- function(AUC,est, amt, Tk0, DV, AUC_cible){
  stopifnot(inherits(est, "mapbayests"))
  CL <- get_param(est, "CL")
  paste0("Patient's clearance: ", round(CL), " L/h.", "\n",
         "AUC estime par le modele: ", round(max(AUC$AUC)/1000, digit = 2), "mg/L.", "\n",
         "Nouvelle dose ", round(AUC_cible * CL / 50) *50, " mg pour une AUC a",AUC_cible,".","\n")}


code <- {"$PARAM @annotated
TVCL : 150 : Cl 
TVVC : 90 : VC

ETA1 : 0 : CL
ETA2 : 0 : VC

$OMEGA 
0.565
0.221

$Sigma
0.01 // proportional
0.01 // additive

$CMT @annotated
CENT  : Central compartment (mg/L) [ADM][OBS]
AUC_cpt : Area under the curve

$TABLE
double DV = (CENT/VC) *(1 + EPS(1)) + EPS(2);
double AUC = AUC_cpt;

$MAIN
double CL = TVVC * exp(ETA1 + ETA(1)) ;
double VC = TVCL * exp(ETA2 + ETA(2)) ;
double KE = CL/VC ;

$ODE @annoted
dxdt_CENT = - CENT * KE ;
dxdt_AUC_cpt = CENT/VC ;

$CAPTURE DV KE CL AUC"}

my_model <- mcode("5FU_model", code)

ui <- dashboardPage(skin = "green",
                    dashboardHeader(title = "Starter-BFC_5FU"),
                    
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Estimation of 5FU IV exposure", tabName = "5FU"),
                        menuItem("Info", tabName = "Info"))),
                    
                    dashboardBody(
                      tabItems(
                        tabItem(tabName = "5FU",
                                fluidRow(
                                  box(height = 600,
                                      title = "Molecule :",solidHeader = T,status = "success",
                                      numericInput("amt", "Dose (mg)", 4800),
                                      numericInput("Tk0", "temps de perfusion", 46),
                                      numericInput("DV", "Concentration mesuree en ug/L", 550),
                                      numericInput("time", "temps entre le debut de la perfusion et le dosage", 8),
                                      numericInput("AUC_cible", "AUC cible", 25),
                                      numericInput("TimeSimu", "Temps de la simulation", 48),
                                      h4("Estimation"),
                                      actionButton("GO", "Let's GO !"),
                                  ),
                                  
                                  box(
                                    title = "Resultats :",solidHeader = T,status = "success",
                                    plotOutput("CV"),verbatimTextOutput("results"),height = 600)),
                                
                                fluidRow(
                                  box(title = "Simu",solidHeader = T,status = "success",width = 12,
                                      numericInput("Newamt", "Nouvelle Dose", 5500),
                                      numericInput("NewTk0", "Temps de perfusion", 46),
                                      actionButton("SimuGO", "Let's GO !"),
                                      verbatimTextOutput("newdose"))),
                        
                        fluidRow(
                          box(title = "Dose versus AUC",solidHeader = T,status = "success",width = 12,
                              sliderInput("range", "Range:",
                                          min = 0, max = 10000,step = 50,
                                          value = c(4800,5800)),
                              dataTableOutput('Dose_vs_AUC')))
                      ),
                      
                        tabItem(tabName = "Info",
                                box(title ="Model",solidHeader = T,status = "success",width = 12,   
                                    paste0("Parameter"),
                                    textOutput("text"),
                                    img(src = "model.png",height="100%", width="100%"))
                                
                        ))))


server <- function(input, output) {
  my_data <- reactive({
    my_model %>%
      adm_lines(time = 0, amt = (input$amt*1000), rate = (input$amt*1000)/(input$Tk0),cmt = 1) %>%
      obs_lines(time = input$time,DV = input$DV) %>%
      get_data() %>% 
      dplyr::filter(!((mdv==0)&(is.na(time)|is.na(DV))))
  })
  
  my_est <- eventReactive(input$GO, {
    mapbayest(my_model, my_data(), verbose = F)
    
  })
  
  AUC <- eventReactive(input$GO, {
    Calcul <- my_est()%>%  use_posterior() %>% 
      data_set(expand.ev(amt = input$amt*1000, dur = input$Tk0) %>% mutate(rate = amt/dur)) %>% 
      carry_out() %>% 
      mrgsim(end = input$Tk0 * 2)
  })
  
  Max.Data <- eventReactive(input$GO,{
    Max.Data <- my_est() %>% 
      use_posterior() %>% 
      data_set(expand.ev(amt = c(seq(min(input$range),max(input$range), by = 50)), dur = c(46)) %>% mutate(rate = amt/dur)) %>% 

      mrgsim(end = input$Tk0 * 2)
})

  
  DF <- eventReactive(input$GO,{
  Max.Data <- Max.Data()
  Dose <- seq(min(input$range),max(input$range), by = 50)
  df.all <- data.frame(Max.Data$ID,round(Max.Data$AUC,digits = 1),Dose)
  colnames(df.all) <- c("ID", "AUC","dose (en mg)")
  DF <- as.data.frame(df.all %>%
    group_by(ID) %>%
    slice(which.max(AUC))) 
  })
  
  output$Dose_vs_AUC <- renderDataTable({
    DF()
    })
    
  NEWAUC <- eventReactive(input$SimuGO, {
    Calcul <- my_est()%>%  use_posterior() %>% 
      data_set(expand.ev(amt = input$Newamt*1000, dur = input$NewTk0) %>% mutate(rate = amt/dur)) %>% 
      carry_out() %>% 
      mrgsim(end = input$Tk0 * 2)
  })
  
  output$CV <- renderPlot({
    my_est()
    My_simu <- my_est()%>%
      augment(start = 0, delta = 1, end = input$TimeSimu,ci =TRUE,ci_method ="delta",ci_width = 95)
    
    simu <- My_simu$aug_tab
    Conc_Pat <- My_simu$aug_tab$value[My_simu$aug_tab$type == "IPRED"]
    time <- My_simu$aug_tab$time[My_simu$aug_tab$type == "IPRED"]
    Ci_low <- My_simu$aug_tab$value_low[My_simu$aug_tab$type == "IPRED"]
    Ci_up <- My_simu$aug_tab$value_up[My_simu$aug_tab$type == "IPRED"]
    mycol <- rgb(0,86,27, max = 255, alpha = 50, names = "green")
    
    plot <-  ggplot(as.data.frame(Conc_Pat), aes(x=time,y=Conc_Pat), height =600) + 
      geom_ribbon(aes(ymin = Conc_Pat, ymax = Ci_up), fill = mycol) +
      geom_ribbon(aes(ymin = Ci_low, ymax = Conc_Pat), fill = mycol) +
      geom_point(simu, mapping =  aes(x=input$time, y=input$DV), size=4,color = "black") +
      geom_line(aes(y = Conc_Pat),size = 1,color="blue")+
      scale_x_continuous("Time (h)", limits = c(0, input$TimeSimu)) +
      scale_y_continuous("Concentrations (µg/L)") + 
      ggtitle(paste0("Pour une perfusion de ",input$amt," mg pendant ",input$Tk0," h."))
    
    print(plot)
  })
  
  output$newdose <- renderText({
    
    A<-NEWAUC()
    paste0("AUC estime par le modele: ", round(max(A$AUC)/1000, digit = 2), "mg/L.")
  })
  output$results <- renderText({
    adapt_poso(
      est = my_est(),
      amt = isolate(input$amt),
      Tk0 = isolate(input$Tk0),
      DV = isolate(input$DV),
      AUC_cible = isolate(input$AUC_cible),
      AUC = AUC()
    )
  })
  
  output$text <- renderText({ 
    HTML(paste("Les parametres du modele :Cl : 80 L.mn")) })
  

}



shinyApp(ui= ui, server = server)
