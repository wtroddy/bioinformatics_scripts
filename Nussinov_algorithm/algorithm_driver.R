# ### algorithm_driver.R
# r script used to drive use of Nussinov algorithms from file "Nussinov_algorithm.R"
# to run this script, put this file and the "Nussinov_algorithm.R" in the same directory
# input sequences are defined by input_seq
# this will store algorithm outputs in the R Environment and print the solutions

### input ----
source("Nussinov_algorithm.R")
input_seq <- "GAAGUGCUUC"

### init the matrix ----
init_list <- Nussinov_algorithm_init_matrix(input_seq)

### run Nussinov algorithm function ----
Nussinov_alg <- Nussinov_algorithm_solver(init_list)

### output solved DP matrix ----
nussinov_matrix <- Nussinov_alg$nussinov_matrix
rownames(nussinov_matrix) <- Nussinov_alg$sequence
colnames(nussinov_matrix) <- Nussinov_alg$sequence
print(nussinov_matrix)

### output traceback  ----
traceback_solution <- Nussinov_traceback(Nussinov_alg$traceback_tracker)
print(traceback_solution)