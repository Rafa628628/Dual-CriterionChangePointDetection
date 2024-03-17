rm(list=ls())
library(TSA)
library(MTS)
library(energy)
library(plotly)
library(dplyr)
library(isotree)
library(shiny)
library(plotly)
library(ggplot2)
library(gridExtra)
library(ggrepel)
# Euclidean distance function
euclideanDistance <- function(x, y) {
  sqrt(sum((x - y)^2))
}

# Custom energy distance function (simplified version)
customEnergyDistance <- function(A, B) {
  dist_AB <- as.matrix(dist(rbind(A, B)))
  nA <- nrow(A)
  nB <- nrow(B)
  dist_AA <- dist_AB[1:nA, 1:nA]
  dist_BB <- dist_AB[(nA+1):(nA+nB), (nA+1):(nA+nB)]
  dist_AB <- dist_AB[1:nA, (nA+1):(nA+nB)]
  
  d = mean(dist_AB) - 0.5 * mean(dist_AA) - 0.5 * mean(dist_BB)
  return(d)
}

# Past Divergence function
pastDivergence <- function(timeSeries, b_1, b_2,d_c, distance) {
  n <- nrow(timeSeries)
  d_vec <- numeric(n)
  
  for (i in (2+b_1):(n-b_2)) {
    if (distance == "Euclidean") {
      distances <- sapply((i-b_1):(i-2), function(j) {
        if(j > 0) euclideanDistance(timeSeries[i, , drop = FALSE], timeSeries[j, , drop = FALSE])
        else 0
      }, simplify = TRUE, USE.NAMES = FALSE)
      d_vec[i] <- min(c(mean(distances, na.rm = TRUE), euclideanDistance(timeSeries[i, , drop = FALSE], timeSeries[i-1, , drop = FALSE])), na.rm = TRUE)
    } else if (distance == "Energy") {
      #for (i in 3:n) {  # Start from 3 to ensure we have A_{i-2} available
      # Define sets A_i, A_{i-1} based on the conditions
      indices_i <- which(sapply(1:n, function(j) abs(j-i) < b_1 && euclideanDistance(timeSeries[i, , drop = FALSE], timeSeries[j, , drop = FALSE]) < d_c))
      indices_i_1 <- which(sapply(1:n, function(j) abs(j-(i-1)) < b_1 && euclideanDistance(timeSeries[i-1, , drop = FALSE], timeSeries[j, , drop = FALSE]) < d_c))
      indices_i_2 <- which(sapply(1:n, function(j) abs(j-(i-2)) < b_1 && euclideanDistance(timeSeries[i-2, , drop = FALSE], timeSeries[j, , drop = FALSE]) < d_c))
      
      # Combine A_i and A_{i-1} for edist
      combined1 <- rbind(timeSeries[indices_i, , drop = FALSE], timeSeries[indices_i_1, , drop = FALSE])
      sizes1 <- c(length(indices_i), length(indices_i_1))  # Sizes of A_i and A_{i-1}
      
      # Combine A_i and A_{i-1} for edist
      combined2 <- rbind(timeSeries[indices_i, , drop = FALSE], timeSeries[indices_i_2, , drop = FALSE])
      sizes2 <- c(length(indices_i), length(indices_i_2))  # Sizes of A_i and A_{i-1}
      
      if (length(indices_i) > 0 && length(indices_i_1) > 0) {
        # Calculate energy distance using edist
        d_vec[i] <- min(edist(combined1, sizes = sizes1, distance = FALSE),edist(combined2, sizes = sizes2, distance = FALSE))
      }
      #}
    }
  }
  
  return(d_vec)
}

# Future Density function
futureDensity <- function(timeSeries,b_1, b_2, d_c, distance) {
  n <- nrow(timeSeries)
  p_vec <- numeric(n)
  
  for (i in (2+b_1):(n-b_2)) {
    futureIndices <- (i+1):min(n, i+b_2)
    p_vec[i] <- sum(sapply(futureIndices, function(j) {
      dist <- if (distance == "Euclidean") euclideanDistance(timeSeries[i, , drop = FALSE], timeSeries[j, , drop = FALSE])
      dist <= d_c
    }, simplify = TRUE, USE.NAMES = FALSE))
  }
  
  return(p_vec)
}
pastDivergenceBiggerThanIt <- function(dat) {
  # Initialize the result vector with default values equal to the number of rows in dat
  df_i <- rep(0, dim(dat)[1])
  
  # Iterate through each row in the dataframe
  for (i in 1:nrow(dat)) {
    # Extract the dis value for the current point
    current_dis <- dat$dis[i]
    
    # Find indices where dis is greater than the current dis
    higher_dis_indices <- which(dat$dis > current_dis)
    
    # Calculate the distance to all points with a higher dis value
    if (length(higher_dis_indices) > 0) {
      # Calculate absolute differences
      distances <- abs(i - higher_dis_indices)
      
      # Find the minimum distance and update df_i for the current point
      df_i[i] <- min(distances)
    }else{
      df_i[i]=max(i,dim(dat)[1]-i)
    }
  }
  
  return(df_i)
}
normedZ <- function(plotdata) {
  # Ensure d and p are treated as numeric
  plotdata$dis <- as.numeric(plotdata$dis)
  plotdata$rho <- as.numeric(plotdata$rho)
  
  # Compute d' using the provided formula
  min_d <- min(plotdata$dis[plotdata$dis > 0])
  max_d <- max(plotdata$dis)
  plotdata$d_prime <- (plotdata$dis - min_d) / (max_d - min_d)
  
  # Compute p' using the provided formula
  min_p <- min(plotdata$rho[plotdata$rho > 0])
  max_p <- max(plotdata$rho)
  plotdata$p_prime <- (plotdata$rho - min_p) / (max_p - min_p)
  
  # Calculate z_i = d'_i * p'_i
  plotdata$distimerho <- ifelse(plotdata$d_prime * plotdata$p_prime< 0, 0, plotdata$d_prime * plotdata$p_prime)*plotdata$disBigger/dim(plotdata)[1]
  plotdata=plotdata[order(plotdata$distimerho,decreasing = T),]
  plotdata$order=seq(1,dim(plotdata)[1],1)
  return(plotdata)
}


ui <- fluidPage(
  titlePanel("Interactive Plot with Click Events"),
  plotlyOutput("interactivePlot"),
  verbatimTextOutput("infoText")
)

server <- function(input, output) {
  output$interactivePlot <- renderPlotly({
    p_plotly
  })
  
  output$infoText <- renderPrint({
    eventdata <- event_data("plotly_click")
    if (!is.null(eventdata)) {
      clickedPoint <- eventdata
      topRightPoints <- plotdata[plotdata$rho > clickedPoint$x & plotdata$dis > clickedPoint$y, ]
      print(topRightPoints)
      
      # Update the plot to highlight top-right points
      p_highlight <- p + geom_point(data = topRightPoints, aes(x = rho, y = dis), color = "blue", size = 3)
      ggplotly(p_highlight)
    }
  })
}

# Example usage:
# Generate sample data

m=1000
n1=1000
ourlies=30
block_before=20
block_after=20
sigma1=0.5
sigma2=2
sigma3=0.5
mu=4
muo=0
thres=0.9
ts.x<-sample(c(rnorm(m,0,sigma1),rnorm(ourlies,muo,sigma2)),m+ourlies)
ts.y<-sample(c(rnorm(n1,0,sigma1),rnorm(ourlies,muo,sigma2)),n1+ourlies)+mu 
ts.x1<-sample(c(rnorm(m,0,sigma1),rnorm(ourlies,muo,sigma2)),m+ourlies)
dat=matrix(c(ts.x,ts.y,ts.x1),ncol=1)


# Initialize plotdata frame
plotdata <- data.frame(point = seq(1, nrow(dat), 1), dis = 0, rho = 0)

# Calculate distance threshold
distance <- as.matrix(dist(dat))
dc <- mean(distance) / 3

# Compute past divergence and future density
plotdata$dis <- pastDivergence(dat, b_1=block_before,block_after,d_c=1,distance="Euclidean")
plotdata$rho <- futureDensity(dat, block_before,block_after, dc,distance="Euclidean")
plotdata$disBigger <- pastDivergenceBiggerThanIt(plotdata)


plotdata1=normedZ(plotdata)

iso <- isolation.forest(data.frame(plotdata1[,c("distimerho")]), ntrees = 100, nthreads = 1)

### print all the change points.
plotdata1$pred <- predict(iso, data.frame(plotdata1[,c("distimerho")]))
print(plotdata1[(plotdata1$pred>thres),"point"])

dat_line=data.frame(t=seq(1,dim(dat)[1],1),x=dat[,1])
line=ggplot(dat_line, aes(x = t, y = x)) +
  geom_line() + # Connect points with lines
  #geom_point() + # Optionally add points
  labs(x = "t", y = "X") +
  theme_bw() +
  theme(panel.grid = element_blank())+
  labs(title="(a) Sequence")
  #theme_minimal()
line
# dim(dat)

# Plot the decision graph using ggplot2
dual_criterion_decision_graph=ggplot(plotdata, aes(x = rho, y = dis)) +
  geom_point(color = "black") +
  geom_text(aes(label = point), hjust = 0, vjust = 0)+
  labs(x = "Future density: p",y = "Past divergence: d")+
  theme_bw() +
  theme(panel.grid = element_blank())+
  labs(title="(b) Dual criterion decision graph")
dual_criterion_decision_graph



z_score <- ggplot(data = plotdata1, aes(x = distimerho, y = 0)) +
  geom_point(size = 3, alpha = 0.5) + # Increase size and adjust alpha
  geom_text_repel(aes(label = ifelse(distimerho > 0.01, as.character(point), NA)), 
                  na.rm = TRUE) +
  labs(x = "Z score", y = "") + # Setting y label to empty
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(), # Remove y-axis labels
        axis.title.y = element_blank()) + # Remove y-axis title
  labs(title="(c) Z score")
z_score

setwd("D:/Work/Work/afterFord/newCP/code/2024Revised")
pdf(paste0("res/twoD_r1.pdf"),width=8,height=4)
bottom_row <- arrangeGrob(line, z_score, ncol = 1)
grid.arrange(bottom_row, dual_criterion_decision_graph, ncol = 2)
dev.off()

pdf(paste0("res/line.pdf"),width=3,height=3)
line
dev.off()
pdf(paste0("res/dual_criterion_decision_graph.pdf"),width=3,height=3)
dual_criterion_decision_graph
dev.off()
pdf(paste0("res/z_score.pdf"),width=3,height=3)
z_score
dev.off()

# # Convert the ggplot object 'p' to a plotly object
# p_plotly <- ggplotly(p)
# # Run the app
# shinyApp(ui = ui, server = server)
