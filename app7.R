library(shiny)
install.packages("mclust")
library(mclust)


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("EM-Algorithm (Gaussian Mixture Model)"),
  
  # Sidebar with inputs
  sidebarLayout(
    sidebarPanel(
      # Input CSV
      shinyUI(fluidPage(
        fileInput('data', 'Choose CSV file',
                  accept=c('text/csv', 'text/comma-separated-values','.csv')),
      )
      ),
      # CSV Header
      checkboxInput("hasHeader", "Header", value = FALSE),
      # CSV Separator 
      radioButtons("separator", "Separator", 
                   choices = c("Comma", "Semicolon", "Tab"),
                   selected = "Comma"
      ),
      # CSV Column to look at
      numericInput("column", "Pick Column:", value = 1, min = 1, max = 10 ),
      # Divider
      div(style = "border-top: 1px solid #ddd; margin-top: 20px; margin-bottom: 20px;"),
      # Number of Modes
      numericInput("modes", "# of Modes:", value = 1, min = 1, max = 10),
      # Divider
      div(style = "border-top: 1px solid #ddd; margin-top: 20px; margin-bottom: 20px;"),
      # EM-Step:
      numericInput("emSteps", "EM-step:", value = 2, min = 1),
      # Divider
      div(style = "border-top: 1px solid #ddd; margin-top: 20px; margin-bottom: 20px;"),
      # Show options
      h5("Show"),
      fluidRow(
        column(4, checkboxInput("showInitialEst", "Initial Est.", value = FALSE)),
        column(4, checkboxInput("showFinalEst", "Final Est.", value = FALSE)),
        column(4, checkboxInput("showLegend", "Legend", value = FALSE))
      ),
      # Number of bins
      sliderInput("bins",
                  "Number of bins:",
                  min = 1,
                  max = 50,
                  value = 30),
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Data", tableOutput("datatable")),  #plots the dataframe of data
                  tabPanel("Model Selection", plotOutput("aicplot")),
                  tabPanel("Plot", plotOutput("distPlot")),  # Show a plot of the generated distribution
                  tabPanel("Summary", tableOutput("line"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  x = NULL
  # Function to compute the probability of each data point belonging to each component
  compute_probabilities <- function(data, mu, sigma) {
    likelihood <- rep(0, length(mu))
    likelihood <- matrix(rep(likelihood, length(data)),ncol=length(mu))
    probability_component <- rep(0, length(mu))
    probability_component <- matrix(
      rep(probability_component, length(data)),ncol=length(mu))
    
    
    for(i in 1:length(mu)){
      likelihood[,i] <- dnorm(data, mean = mu[i], sd = sigma[i])
    }
    total_likelihood <- rowSums(likelihood)
    
    for(i in 1:length(mu)){
      probability_component[,i] <- likelihood[,i] / total_likelihood
    }
    return(probability_component)
  }
  
  # Function to update the parameters using the computed probabilities
  update_parameters <- function(data, probabilities) {
    sum_prob<-rep(0, dim(probabilities)[2])
    mu<-rep(0, dim(probabilities)[2])
    sigma<-rep(0, dim(probabilities)[2])
    
    for(i in 1:dim(probabilities)[2]){
      sum_prob[i] <- sum(probabilities[,i])
      mu[i] <- sum(probabilities[,i] * data) / sum_prob[i]
      sigma[i] <- sqrt(sum(probabilities[,i] * (data - mu[i])^2) / sum_prob[i])
    }
    return(list(mu = mu, sigma = sigma))
  }
  
  
  em_algorithm <- function(data, num_iterations, modes) {
    epsilon <- 0.000001
    if (is.null(modes)){
      modes<-input$modes
    }
    
    # Initialize parameters
    mem <- kmeans(data, modes)$cluster
    mu <- c()
    sigma <- c()
    
    for(i in 1:modes)
    {
      mu <- c(mu,mean(data[mem==i]))
      sigma <- c(sigma,sd(data[mem==i]))
    }
    
    var1 <- c()
    
    for (iteration in 1:num_iterations) {
      # E-step: Compute probabilities
      probabilities <- compute_probabilities(data, mu, sigma)
      
      # M-step: Update parameters
      parameters <- update_parameters(data, probabilities)
      
      
      var2 <- ((norm(t(mu),"2")+norm(t(sigma),"2")) - (norm(t(parameters$mu),"2")+norm(t(parameters$sigma),"2")))**2
      
      if(var2 <epsilon**2){
        globalValues$convergence <-iteration
        #updateNumericInput(session, "emSteps", value = iteration)
        break
      }
      
      # Update parameters for the next iteration
      mu <- parameters$mu
      sigma <- parameters$sigma
      
      var1 <- rbind(var1,c(mu=mu,sigma=sigma))
      
    }
    globalValues$result <- var1
    return(list(var1))
  }
  
  globalValues <- reactiveValues(my_modes = 1,
                                 result = NULL,
                                 emsteps=2,
                                 convergence = NULL)
  observeEvent(input$modes, {
    globalValues$my_modes <- input$modes
    x <- data()
    x <- x[,input$column]
    em_algorithm(x, 1000, globalValues$my_modes)
    
  })
  observeEvent(input$emSteps, {
    globalValues$emsteps <- input$emSteps
  })
  
  #generates or takes data that is uploaded by user
  data <- reactive({
    inFile <- input$data
    if (is.null(inFile) && is.null(x)){
      
      # Generate some data from a mixture of two normal distributions
      set.seed(1)
      random_data <- c(rnorm(500, mean = 5, sd = 2), rnorm(500, mean = 15, sd = 3))
      random_data2 <- c(rnorm(300, mean = 5, sd = 2), rnorm(300, mean = 15, sd = 3), rnorm(300, mean = 30, sd = 4))
      random_data3 <- c(rnorm(300, mean = 5, sd = 2), rnorm(300, mean = 15, sd = 3), rnorm(300, mean = 30, sd = 4), rnorm(300, mean = 10, sd = 1))
      data <- as.data.frame(cbind(random_data, random_data2, random_data3))
    }
    else {
      #based off what the separator is and if there is a header
      #it reads in file that way
      seperator <- switch(input$separator,
                          "Comma" = ",",
                          "Semicolon" = ";",
                          "Tab" = "\t"
      )
      read.csv(inFile$datapath, sep = seperator, header = input$hasHeader)
    }  
    
  })
  
  #creates a histogram plot
  output$distPlot <- renderPlot({
    chosen <- globalValues$emsteps
    
    #makes sure chosen doesn't go larger than convergence
    if(globalValues$convergence <= chosen)
    {
      chosen <- globalValues$convergence -1
    }
    
    x <- data()
    x <- x[,input$column]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, prob = TRUE, col = 'darkgray', border = 'white',
         xlab = 'Data',
         ylab = 'Density',
         main = 'Histogram and GMM fits')
    
    t<- seq(min(x), max(x), 0.01)
    
    ## NEED TO KNOW THE NUMBER OF STEP PLOTTING
    density_init <- t*0
    density_final <- t*0
    density_chosen <- t*0
    
    #CHANGE: IF NULL THEN NEED TO RUN EM
    result <- globalValues$result
    
    
    #print(result)
    
    for(i in 1:globalValues$my_modes)
    {
      #NEED TO CHANGE 1 TO THE VALUE WE WANT TO PLOT
      density_init<- density_init + dnorm(t, result[1,i], result[1,i+globalValues$my_modes])
      density_final<- density_final + dnorm(t, result[length(result[,1]),i], result[length(result[,1]),i+globalValues$my_modes])
      #print(length(result[,1]))
      density_chosen<- density_chosen + dnorm(t, result[chosen,i], result[chosen,i+globalValues$my_modes])
      
    }
    
    # chosen estimate
    lines(t,density_chosen/globalValues$my_modes, lwd = 2, col = 'black')
    if (input$showInitialEst){
      # initial estimate
      lines(t,density_init/globalValues$my_modes, lwd = 2, lty = 2, col = 'blue')
    }
    if (input$showFinalEst){
      # final estimate
      lines(t,density_final/globalValues$my_modes, lwd = 2, lty = 2, col = 'red')
    }
    
    
    # CHANGE: SO IT SAYS STEP #
    if (input$showLegend){
      #legend
      if (input$showInitialEst & input$showFinalEst){ legend('topright', legend = c('Initial', 'Final', paste('Step', chosen)), col = c('blue', 'red', 'black'), lty = c(2, 2, 1), title = "Legend")}
      else if (input$showInitialEst){ legend('topright', legend = c('Initial', paste('Step', chosen)), col = c('blue', 'black'), lty = c(2, 1), title = "Legend")}
      else if (input$showFinalEst){ legend('topright', legend = c('Final', paste('Step', chosen)), col = c('red', 'black'), lty = c(2, 1), title = "Legend")}
      else { legend('topright', legend = c(paste('Step', chosen)), col = c('black'), lty = 1, title = "Legend")}
    }
    
  })
  
  # creates a table of the data
  output$datatable <- renderTable({
    data()  # Use the reactive data instead of reading from inFile$datapath
  }, rownames = FALSE, digits = 4)
  
  # creates plot of AIC & BIC
  output$aicplot <- renderPlot({
    x <- data()
    x <- x[,input$column]
    
    em<- c()
    
    n_components <- 1:10
    clfs <- lapply(n_components, function(n) Mclust(x, G = n, modelNames = "V"))
    bics <- sapply(clfs, function(clf) BIC(clf))
    aics <- sapply(clfs, function(clf) AIC(clf))
    
    #only return AIC because... ask us
    globalValues$my_modes <- which.min(aics)
    updateNumericInput(session, "modes", value = globalValues$my_modes)
    #print(globalValues$my_modes)
    
    plot(n_components, aics, type = 'b', lwd = 2, col = 'blue', ylab = 'BIC/AIC', xlab = '# of Modes', main = 'Information Criterion', ylim = c(min(aics), max(bics)))
    lines(n_components, bics, type = 'b', lwd = 2, col = 'red')
    legend('topright', legend = c('AIC', 'BIC'), col = c('blue', 'red'), lty = 1)
    
    
  })
  
  
  # EM - ALgorithm
  output$line <- renderTable({
    # input pi as vector of probabilites,
    # modes as integer
    x <- data()
    x <- x[,input$column]
    
    
    # Run the EM algorithm
    num_iterations <- 1000
    
    
    #print(globalValues$my_modes)
    result <- em_algorithm(x, num_iterations, globalValues$my_modes)
    
    # Print the results
    #print("Estimated Parameters:")
    #print(result)
    
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)


