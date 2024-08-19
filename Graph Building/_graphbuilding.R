library(igraph)
library(dplyr)
library(ggplot2)
library(tidyverse)

### HELPER FUNCTIONS ###

find_nearest_node <- function(nodeID, d_mat) {
  # Return node ID that is closest to the given node in procedure matrix
  # p_mat = procedure matrix, d_mat = distance matrix
  
  # identify shortest distance in the distance matrix
  min_dist <- min(slice(d_mat, nodeID),na.rm=TRUE)
  
  # obtain node ID with the shortest distance
  return(as.numeric(which(d_mat[nodeID] == min_dist)[1]))
}

find_connected_nodes <- function(nodeID, adj_mat){
  return(which(slice(adj_mat, nodeID) == 1))
}

# Check if newly formed node would form a small angle (i.e. <45 deg)
new_edge_forms_small_angle <- function(node0, node1, adj_mat, p_mat) {
  # identify connected nodeIDs
  connected_nodes <- find_connected_nodes(node0, adj_mat)
  # construct vector to compare
  uv <- c(p_mat$X[node1], p_mat$Y[node1]) - c(p_mat$X[node0], p_mat$Y[node0])
  dist_uv <- sqrt(sum(uv*uv))
  
  # compare each vector
  for(n in connected_nodes){
    # compare vector
    uw <- c(p_mat$X[n], p_mat$Y[n]) - c(p_mat$X[node0], p_mat$Y[node0])
    dist_uw <- sqrt(sum(uw*uw))
    uv_dot_uw <- sum(uv*uw)
    if(acos(uv_dot_uw/(dist_uv*dist_uw)) <= pi/3){
      return(TRUE)
    }
  }
  return(FALSE)
}

# Check if two nodes are connectable
is_connectable <- function(node1, node2, d_mat, p_mat, adj_mat){
  # Priority 1: check if nodes can't accommodate any more connections
  if(p_mat$d[node1] >= 3 || p_mat$d[node2] >= 3){
    return(FALSE)
  }
  # Priority 2: check if newly formed node would form a small angle (i.e. <45 deg)
  else if(new_edge_forms_small_angle(node1, node2, adj_mat, p_mat)){ 
    return(FALSE) 
  }
  # Priority 3: check if nodes are already connected and are close to each other
  else if(!is.na(d_mat[node1, node2])){
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

# infers a honeycomb lattice graph given a file of coordinates (ID, X, Y)
build_graph <- function(fname, n_max_iter=10000, start_node=1){
#fname is a csv with absolute X and Y coordinates
  
  coord_list <- read.csv(fname, header=TRUE, sep=",") %>% 
    rename("ID" = "X.1")
  # plot(coord_list$X, coord_list$Y)
  
  nnodes <- as.numeric(nrow(coord_list))
  
  # calculate distance between each single node
  dist_matrix <- as.data.frame(as.matrix(dist(select(coord_list,X,Y))))
  dist_matrix <- replace(dist_matrix, dist_matrix == 0, NA)
  
  # create adjacency matrix
  adj_matrix <- data.frame(matrix(0,ncol = nnodes, nrow = nnodes))
  # create weighted adjacency matrix
  weighted_adj_matrix <- adj_matrix
  colnames(adj_matrix) <- c(1:nnodes)
  
  # initiate the procedure matrix
  # ID, X, Y, d=degree
  proc_matrix <- data.frame(coord_list,
                            d = 0)
  
  ## INIT ##
  # max iterations
  
  # set starting node
  curr_node <- start_node
  # initiate the list of visited nodes
  visited <- list()
  
  ## LOOP ##
  
  for(ii in 1:n_max_iter){
    # find the nearest node
    nearby_node <- find_nearest_node(curr_node, dist_matrix)
    
    # if nearby node is connectable, do the following:
    if(is_connectable(node1 = curr_node,
                      node2 = nearby_node,
                      d_mat = dist_matrix,
                      p_mat = proc_matrix,
                      adj_mat = adj_matrix)){
      # update the adjacency matrix
      adj_matrix[curr_node, nearby_node] <- 1
      adj_matrix[nearby_node, curr_node] <- 1
      
      # update the weighted matrix
      d <- dist_matrix[curr_node, nearby_node]
      weighted_adj_matrix[curr_node, nearby_node] <- d
      weighted_adj_matrix[nearby_node, curr_node] <- d
      
      # increase # of connections by 1
      proc_matrix$d[curr_node] <- proc_matrix$d[curr_node]+1
      proc_matrix$d[nearby_node] <- proc_matrix$d[nearby_node]+1
      
      # update the distance matrix so now it's NA
      # we do this so the min can find the next smallest distance
      dist_matrix[curr_node, nearby_node] <- NA
      dist_matrix[nearby_node, curr_node] <- NA
      
      # we are done visiting this node
      visited <- c(curr_node, visited)
      
      # update the next node
      next_node <- nearby_node
      
    } else {
      
      # in the case that this is the last node we are visiting
      if(length(visited) == 1){
        # we are done!
        # break
        cat("Finished constructing graph in",ii,"iterations \n")
        break
      }
      
      # set the new node to be the previous node
      next_node <- as.numeric(visited[2])
      
      # truncate our visited list at the head
      visited <- visited[-1]
    }
    
    # update the current node to the next node
    curr_node <- next_node
    
  }
  
  g = graph_from_adjacency_matrix(adjmatrix = as.matrix(adj_matrix), 
                                     mode="undirected")
  V(g)$name <- proc_matrix$ID
  V(g)$X <- proc_matrix$X
  V(g)$Y <- proc_matrix$Y
  
  return(g)
}

# new ROI box object
new_boxROI <- function(X.abs=0, Y.abs=0, X.local=6, Y.local=6, width=12, height=12){
  return(list(X = X.abs,
              Y = Y.abs,
              X.local = X.local,
              Y.local = Y.local,
              width = width,
              height = height
              ))
}

# convert local position to absolute position
# (pos, box) -> vector
local_box_pos_to_absolute_pos <- function(local_pos, boxROI=new_boxROI()){
  offset <- c(boxROI$X.local, boxROI$Y.local) - local_pos
  return(c(boxROI$X, boxROI$Y) - offset)
}

# determines whether point is in the boxROI
# (pos, boxROI) -> bool
is_in_boxROI <- function(abs_pos, boxROI=new_boxROI()){
  pos_x <- abs_pos[1]
  pos_y <- abs_pos[2]
  x_lower <- local_box_pos_to_absolute_pos(c(1,1),boxROI)[1]
  y_lower <- local_box_pos_to_absolute_pos(c(1,1),boxROI)[2]
  x_higher <- local_box_pos_to_absolute_pos(c(boxROI$width,boxROI$height),boxROI)[1]
  y_higher <- local_box_pos_to_absolute_pos(c(boxROI$width,boxROI$height),boxROI)[2]
  return(pos_x >= x_lower && pos_x <= x_higher &&
           pos_y >= y_lower && pos_y <= y_higher)
}

# We want to map a list of box ROIs to a coordinate on the graph
# then remove all the unmapped coordinates
build_graph_from_boxROI_coords <- function(referenceGraph, box_ROI_coords_file){
  g <- referenceGraph
  # read the boxROI coord file
  box_ROI_coords <- read.csv(box_ROI_coords_file, header=TRUE, sep=",") %>% 
    rename("ID" = "X.1")
  
  # make a list of boxes, then assign each box
  box_ROI_list <- list()
  for(ID in box_ROI_coords$ID){
    b <- new_boxROI(X.abs = box_ROI_coords$X[ID],
                    Y.abs = box_ROI_coords$Y[ID])
    box_ROI_list <- append(box_ROI_list, list(b))
  }
  
  # matching algorithm
  
  for(v in V(g)){
    # get coords from graph
    v_coords <- c(V(g)$X[V(g) == v], V(g)$Y[V(g) == v])
    
    for(boxROI in box_ROI_list){
      # get coords from box_ROI_list
      if (is_in_boxROI(v_coords,boxROI)){
        # replace the coords of the graph
        V(g)$boxROI[V(g) == v] <- list(boxROI)
      }
    }
  }
  
  g_new <- g
  for(ii in 1:length(V(g))){
    if(is.null(V(g)$boxROI[[ii]])){
      g_new <- delete_vertices(g_new, V(g_new)[V(g_new)$name == V(g)[ii]])
    }
  }
  
  return(g_new)
  
}

