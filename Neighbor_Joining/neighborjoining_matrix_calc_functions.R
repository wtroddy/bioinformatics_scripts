# calculate the prime matrix
calc_d_prime <- function(input_matrix, d_star_matrix, min_pos_col, min_pos_row) {
    # get stats from input matrix
    n <- ncol(input_matrix)
    min_distance <- input_matrix[min_pos_row, min_pos_col]
    total_d_col <- sum(input_matrix[, min_pos_col])
    total_d_row <- sum(input_matrix[min_pos_row,])
    delta_min <- ((total_d_col - total_d_row) / (n - 2))
    min_distance_x <- ((min_distance + delta_min) / 2)
    min_distance_y <- ((min_distance - delta_min) / 2)
    
    # set default values
    cur_col <- 1
    cur_row <- 1
    d_prim_vals <- c()
    
    # loops 
    for (i in 1:ncol(input_matrix)) {
      for (j in 1:nrow(input_matrix)) {
        # logic to see if this is a collapsed neighbor
        if (cur_col == min_pos_col) {
          collapsed_value <- ((input_matrix[cur_row, min_pos_col] + input_matrix[cur_row, min_pos_row] - min_distance) / 2)
          d_prim_vals <- append(d_prim_vals, collapsed_value, after = length(d_prim_vals))
          # increment and repeat
          cur_row <- cur_row + 1
        } else {
          # increment ncol for loop 
          cur_row <- cur_row + 1
        }
      }
      cur_row <- 1
      cur_col <- cur_col + 1
    }
    
    ### create new matrix
    d_prime_matrix <- matrix(nrow = (nrow(input_matrix)-1), ncol = (ncol(input_matrix)-1))
    repl_row_col <- min(min_pos_col, min_pos_row)
    orig_input_matrix <- input_matrix[-min_pos_row:-min_pos_col,-min_pos_row:-min_pos_col]
  
    #fill in loop
    for (i in 1:ncol(d_prime_matrix)) {
      for (j in 1:nrow(d_prime_matrix)) {
        if (i == repl_row_col || j == repl_row_col) {
          d_prime_matrix[,repl_row_col] <- d_prim_vals[-min_pos_col]
          d_prime_matrix[repl_row_col,] <- d_prim_vals[-min_pos_col]
        } else
          if (i > repl_row_col || j > repl_row_col) {
            if (i > repl_row_col && j < repl_row_col) {
              d_prime_matrix[j,i] <- input_matrix[j,i+1]
            } else if (i < repl_row_col && j > repl_row_col ) {
              d_prime_matrix[j,i] <- input_matrix[j+1,i]
            } else {
              d_prime_matrix[j,i] <- input_matrix[j+1,i+1]
            }
            
          } else if (i < repl_row_col || j < repl_row_col) {
            d_prime_matrix[j,i] <- input_matrix[j,i]
          }
      }
    }
    
    # create ouput list
    calc_d_prime_list <- list(
      min_distance = min_distance,
      total_d_col = total_d_col,
      total_d_row = total_d_row,
      delta_min = delta_min,
      min_distance_x = min_distance_x,
      min_distance_y = min_distance_y,
      d_prime_collapsed_values = d_prim_vals,
      d_prime_matrix = d_prime_matrix
    )
    
    # return the list
    return(calc_d_prime_list)
}

# calculate d star 
calc_d_star <- function(input_matrix) {
  # set default values 
  cur_col <- 1
  cur_row <- 1
  n <- ncol(input_matrix)
  d_star_vals <- c()
  
  # loop through each position
  for(i in 1:ncol(input_matrix)){
    for (j in 1:nrow(input_matrix)) {
      d_star_pos <- (n-2)*input_matrix[cur_col,cur_row]-(sum(input_matrix[cur_col,]))-sum(input_matrix[cur_row,])
      d_star_vals <- append(d_star_vals, d_star_pos, after=length(d_star_vals))
      cur_row <- cur_row+1
    }
    cur_row <- 1
    cur_col <- cur_col+1
  }
  
  # create matrix
  d_star_matrix <- matrix(d_star_vals, nrow =  nrow(input_matrix), ncol =  ncol(input_matrix), byrow = TRUE)
  
  # remove identity positions 
  for(i in 1:ncol(input_matrix)) {
    d_star_matrix[i,i] <- 0L
  }
  
  # return output matrix
  return(d_star_matrix)
}
