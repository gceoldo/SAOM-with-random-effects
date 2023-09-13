library(igraph)
library(purrr)

igw1 <- graph_from_adjacency_matrix(adj.matrix.1, "directed")
V(igw1)$name <- 1:39
igw2 <- graph_from_adjacency_matrix(adj.matrix.2, "directed")
V(igw2)$name <- 1:39

oserved_out_degrees_w1 <- map_int(1:39, ~ sum(adj.matrix.1[.x,]))
oserved_out_degrees_w2 <- map_int(1:39, ~ sum(adj.matrix.2[.x,]))

N <- 39
scale_degrees <- 0.03
lim_graph <- 2.2
bound_right <- lim_graph - .1
bound_left <- -lim_graph + .1
bound_top <- .9
space_between_segments <- 0.05
width_segments <- 4
how_above_text <- .2
how_below_text <- .1
degrees_labelled <- c(0,10,20)
labels_distribution <- c("out-degree","out-degree")
pos_x_left_text <- -1
pos_x_right_text <- 1
graphs_to_plot <- list(igw1, igw2)
degrees_to_plot <- list(oserved_out_degrees_w1, oserved_out_degrees_w2)
text_graph_description <- list('wave 1','wave 2')
next_line <- .1

function_plot <- function() {
  #l <- layout.mds(graphs_to_plot[[1]])
  #l <- layout.kamada.kawai(graphs_to_plot[[1]])
  l <- layout.fruchterman.reingold(graphs_to_plot[[2]])
  #l <- layout.circle(graphs_to_plot[[1]])
  
  for(k in 1:2) {
    if(k%%2 == 1) { # left plot
      
      plot(graphs_to_plot[[k]],
           edge.arrow.size = .4, 
           vertex.size = 6,
           vertex.label.cex = 0.7,
           vertex.label = "",
           vertex.color = 1 + high.status, 
           edge.color = "black",
           vertex.shape = c("circle"),#, "square")[2-high.status],
           xlim=c(-1,lim_graph), layout=l)
      
      abline(h=1.3)
      abline(h=-1.3)
      
      text(bound_right-scale_degrees*degrees_labelled[2],
           bound_top+how_above_text,
           labels_distribution[k])
      
      for(j in 1:length(degrees_labelled)) {
        segments(bound_right-scale_degrees*degrees_labelled[j], 
                 bound_top+space_between_segments,
                 bound_right-scale_degrees*degrees_labelled[j], 
                 bound_top-space_between_segments*N, 
                 lty=c(1,3,3)[j])
        
        text(bound_right-scale_degrees*degrees_labelled[j], 
             bound_top-space_between_segments*N-how_below_text, 
             as.character(degrees_labelled[j]))
        
        for(i in 1:N) {
          #if(i %in% c(1,)) text(bound_right,bound_top - (i-1)*space_between_segments,i)
          if(degrees_to_plot[[k]][i]>0) {
            segments(bound_right-scale_degrees*degrees_to_plot[[k]][i], 
                     bound_top - (i-1)*space_between_segments,
                     bound_right, 
                     bound_top - (i-1)*space_between_segments,
                     lwd=width_segments, col = categorical_pal(3)[1+high.status[i]])
          }
        }
      }
      text(pos_x_left_text,
           bound_top+how_above_text,
           text_graph_description[[k]][1])
      text(pos_x_left_text,
           bound_top+how_above_text-next_line,
           text_graph_description[[k]][2])
      text(pos_x_left_text,
           bound_top+how_above_text-2*next_line,
           text_graph_description[[k]][3])
    }
    else {
      plot(graphs_to_plot[[k]],
           edge.arrow.size = .4, 
           vertex.size = 6,
           vertex.label.cex = 0.7,
           vertex.label="",
           vertex.color = 1 + high.status, 
           edge.color = "black",
           vertex.shape = c("circle"),#, "square")[2-high.status],
           xlim=c(-lim_graph,1), layout=l)
      
      abline(h=1.3)
      abline(h=-1.3)
      
      text(bound_left+scale_degrees*degrees_labelled[2],
           bound_top+how_above_text,
           labels_distribution[k])
      
      for(j in 1:length(degrees_labelled)) {
        text(bound_left+scale_degrees*degrees_labelled[j], 
             bound_top-space_between_segments*N-how_below_text, 
             as.character(degrees_labelled[j]))
        
        segments(bound_left+scale_degrees*degrees_labelled[j], 
                 bound_top+space_between_segments,
                 bound_left+scale_degrees*degrees_labelled[j], 
                 bound_top-space_between_segments*N, 
                 lty=c(1,3,3)[j])
        
        for(i in 1:N) {
          if(degrees_to_plot[[k]][i]>0) {
            segments(bound_left, 
                     bound_top - (i-1)*space_between_segments,
                     bound_left+scale_degrees*degrees_to_plot[[k]][i], 
                     bound_top - (i-1)*space_between_segments,
                     lwd=width_segments, col = categorical_pal(3)[1+high.status[i]])
          }  
        } 
      } 
      text(pos_x_right_text,
           bound_top+how_above_text,
           text_graph_description[[k]][1])
      text(pos_x_right_text,
           bound_top+how_above_text-1*next_line,
           text_graph_description[[k]][2])
      text(pos_x_right_text,
           bound_top+how_above_text-2*next_line,
           text_graph_description[[k]][3])
      
      
    }
  }
}


scalepdf <- 3
pdf("data.pdf",width=3 * scalepdf,height=3 * scalepdf)

par(mfrow=c(1,2))
par(mar=c(1,1,1,1))

function_plot()

dev.off()

