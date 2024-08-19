source("_graphbuilding.R")

## HELPER FUNCTIONS ##

is_positively_oriented <- function(v1, v2){
  cross_product <- v1[1]*v2[2] - v1[2]*v2[1] 
  if(cross_product > 0){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

draw_ring_path <- function(starting_tail_node, starting_head_node, g){
  
  # initialize adjacency matrix
  adjm <- as.data.frame(as_adjacency_matrix(g, sparse=FALSE))
  
  # check if starting nodes are connected
  if(adjm[starting_tail_node, starting_head_node] == 0){
    return(NULL)
  }
  
  tail <- starting_tail_node
  head <- starting_head_node
  
  # initialize path
  p <- c(starting_tail_node, starting_head_node)
  
  max_iter <- 15
  
  for(ii in 1:max_iter){
    
    # check if head is the starting tail
    if(head == starting_tail_node){
      return(p)
    }
    
    # initialize base vector
    v1 <- c(V(g)$X[head], V(g)$Y[head]) - c(V(g)$X[tail], V(g)$Y[tail])
    
    # find neighboring nodes
    neighbors <- find_connected_nodes(head, adjm)
  
    # find the node whose path is positively oriented and is not the previous one
    for(neighbor in neighbors){
      v2 <- c(V(g)$X[neighbor], V(g)$Y[neighbor]) - c(V(g)$X[head], V(g)$Y[head])
      if(is_positively_oriented(v1, v2)){
        # update path
        p <- c(p, neighbor)
        
        # update players
        tail <- head
        head <- neighbor
        break
      }
    }
  }
  return(NULL)
}

calculate_ring_perimeter <- function(path, g){
  # check if path is circular
  if((path[1] != path[length(path)]) || is.null(path)){
    return(NULL)
  }
  
  head <- path[2]
  tail <- path[1]
  d <- 0
  
  # LOOP #
  for(ii in 1:(length(path)-1)){
    # create the vector
    v <- c(V(g)$X[head], V(g)$Y[head]) - c(V(g)$X[tail], V(g)$Y[tail])
    
    # calculate its distance
    d <- d + sqrt(sum(v*v))
    
    # update players
    head <- path[ii+2]
    tail <- path[ii+1]
  }
  
  return(d)
}

calculate_ring_area <- function(path, g){
  # check if path is circular
  if((path[1] != path[length(path)]) || is.null(path)){
    return(NULL)
  }
  
  head <- path[2]
  tail <- path[1]
  a <- 0
  
  # LOOP #
  for(ii in 1:(length(path)-1)){
    # create the vector
    coord1 <- c(V(g)$X[tail], V(g)$Y[tail])
    coord2 <- c(V(g)$X[head], V(g)$Y[head])

    # calculate its determinant
    a <- a + coord1[1]*coord2[2] - coord1[2]*coord2[1]
    
    # update players
    head <- path[ii+2]
    tail <- path[ii+1]
  }
  
  return(a/2)
}

calculate_shape_factor <- function(path,g){
  P <- calculate_ring_perimeter(path, g)
  A <- calculate_ring_area(path, g)
  if(is.null(path)){
    return(NULL)
  }
  return(4*pi*A/(P*P))
}

calculate_ring_centroid_from_path <- function(path, g){
  # check if path is circular
  if((path[1] != path[length(path)]) || is.null(path)){
    return(NULL)
  }
  v <- c(0,0)
  # LOOP #
  for(ii in 1:(length(path)-1)){
    # add to the vector
    v <- v + c(V(g)$X[path[ii]], V(g)$Y[path[ii]])
  }
  
  # averages all the coordinates
  return(v/(length(path)-1))
}
