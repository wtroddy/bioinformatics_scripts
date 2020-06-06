# ### neighbor_joining_driver.R
# r script used to drive use of functions for matrix calculations 
# to run this script, 
# 	1. put this file and the "neighborjoining_matrix_calc_functions.R" in the same directory
# 	2. define input values and settings in code below
#
# output:
#	this will write the matrix at each iteration to the cur wd as a csv
#
# to do's:
#	automate the iteration process with loops
#	add customization of what should be output and saved to directory 
#   wrap everything in a single function with input 

# load functions
source("neighborjoining_matrix_calc_functions.R")

# Input Values and Settings 
  input_vals <- c(0,6,9,11,15,6,0,11,13,17,9,11,0,12,16,11,13,12,0,10,15,17,16,10,0)
  input_rows <- 5
  input_cols <- 5
  input_edge_names <- c("A", "B", "C", "D", "E")

# create matrix
inmat <- matrix(input_vals, nrow=input_rows, ncol = input_cols, byrow = TRUE)

### first iteration ---
  d1_star <- calc_d_star(inmat)
  # d prime 1, use position [4,5] 
  d1_prime <- calc_d_prime(inmat, d1_star, 4, 5)
  
  # add edge names
  d1_edge_names <- c("A", "B", "C", "M")
  d1_prime_matrix_edge_names <- d1_prime$d_prime_matrix
  write.csv(d1_prime_matrix_edge_names, "d1_prime_matrix_edge_names.csv", row.names = d1_edge_names)

  
### second iteration ----
  d2_star <- calc_d_star(d1_prime$d_prime_matrix)
  # d prime 2, use 3,4
  d2_prime <- calc_d_prime(d1_prime$d_prime_matrix, d2_star, 3, 4)
  
  # add edge names
  d2_edge_names <- c("A", "B", "M_2")
  d2_prime_matrix_edge_names <- d2_prime$d_prime_matrix
  row.names(d2_prime_matrix_edge_names) <- d2_edge_names
  write.csv(d2_prime_matrix_edge_names, "d2_prime_matrix_rownames.csv")
  

### third iteration ----
  d3_star <- calc_d_star(d2_prime$d_prime_matrix)
  # d prime 3, use 1,3
  d3_prime <- calc_d_prime(d2_prime$d_prime_matrix, d3_star, 2, 3)

  # add edge names
  d3_edge_names <- c("A", "M_3")
  d3_prime_matrix_edge_names <- d3_prime$d_prime_matrix
  row.names(d3_prime_matrix_edge_names) <- d3_edge_names
  write.csv(d3_prime_matrix_edge_names, "d3_prime_matrix_edge_names.csv")
  