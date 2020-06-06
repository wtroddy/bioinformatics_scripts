# ### Nussinov_algorithm.R 
# This file includes three functions to execute and 'solve' a Nussinov algorithm
#
# The functions are expected to be run in this order, with the noted input:
#    Order | Function Name                   | Purpose                           | Input            | Output                  |
#   --------------------------------------------------------------------------------------------------------------------------|
#          | Nussinov_algorithm_init_matrix  | produce a blank matrix with zeros | Sequence string  | List including          |
#      1.  |                                 | and handle formatting of          | to be predicted  | matrix and 'clean'      |
#          |                                 | the input sequence string         |                  | string.                 |
#   --------------------------------------------------------------------------------------------------------------------------|
#          | Nussinov_algorithm_solver       | complete the dynamic programming  | List from 1.     | List including complete |
#      2.  |                                 | table with the maximum score      | _init_matrix     | DP table matrix and     |
#          |                                 |                                   |                  | traceback shadow matrix |
#   --------------------------------------------------------------------------------------------------------------------------|
#          | Nussinov_traceback              | produce a table with each step    | trackback shadow | data.frame of each step |
#      3.  |                                 | taken and "traceback" from which  | matrix from 2.   | and the associated calcs|
#          |                                 | score lead to the final pathway   | _algorithm_solver| from the top to bottom  |
#          |                                 | this is a recurisve function      |                  |                         |
#   --------------------------------------------------------------------------------------------------------------------------|


### function to initialize matrix ----
Nussinov_algorithm_init_matrix <-
  function(input_seq) {
    ### check sequence format ----
    if (length(input_seq) == 1 && nchar(input_seq) > 1) {
      seq_split <- strsplit(input_seq, split = "")[[1]]
    } else {
      seq_split <- input_seq;
    }
    
    ### setup matrix ----
    nussinov_matrix <-
      matrix(NA, nrow = length(seq_split), ncol = length(seq_split))
    
    #### initialize the matrix with zeros ----
    for (i in 1:nrow(nussinov_matrix)) {
      for (j in 1:ncol(nussinov_matrix)) {
        if (i == j) {
          nussinov_matrix[i, j] <- 0
        } else if ((i - 1) == (j)) {
          nussinov_matrix[i, j] <- 0
        }
      }
    }
    return(list(sequence = seq_split, nussinov_matrix = nussinov_matrix))
  }

### function to implement recursive algorithm ----
Nussinov_algorithm_solver <-
  function(nussinov_list, randomize_direction = TRUE) {
    ### setup ----
    ref_basepairs <-
      list(c("A", "U"), c("U", "A"), c("C", "G"), c("G", "C"))
    seq_split <- nussinov_list$sequence
    nussinov_matrix <- nussinov_list$nussinov_matrix
    
    ### define traceback tracker df to fill ----
    traceback_tracker <- data.frame(
      cur_diag = character(),
      i_posit = character(),
      j_posit = character(),
      ij_posit = character(),
      i_value = character(),
      j_value = character(),
      ij_value = character(),
      mi_score = character(),
      mj_score = character(),
      md_score = character(),
      mk_score = character(),
      max_dir = character(),
      ref_posit = character(),
      max_val = character(),
      mat_dir = character(),
      stringsAsFactors = FALSE
    )
    
    ### fill in nussinov_matrix ----
    # set cur diag tracker
    cur_diag = 1
    # loop and fill
    for (L in 1:length(seq_split)) {
      for (j in 1:ncol(nussinov_matrix)) {
        for (i in 1:nrow(nussinov_matrix)) {
          if (is.na(nussinov_matrix[i, j])) {
            if ((cur_diag - abs(i - j)) == 1) {
              if (i < j) {
                #calc if d is a match
                if (list(c(seq_split[i], seq_split[j])) %in% ref_basepairs) {
                  d <- 1
                } else {
                  d <- 0
                }
                
                ### calculate max direction ----
                # mi and md
                if (i < nrow(nussinov_matrix)) {
                  mi <- nussinov_matrix[i + 1, j]
                  md <- (nussinov_matrix[i + 1, j - 1]) + d
                }
                # mj
                mj <- nussinov_matrix[i, j - 1]
                
                ### bifurcation
                if (i < nrow(nussinov_matrix) && j < ncol(nussinov_matrix)) {
                  # init values
                  k_values <- seq((i+1),(j-1))
                  mk_values <- data.frame(k = c(),
                                          val = c())
                  # get all values of k
                  for (k in k_values){
                    #mk_values <- append(mk_values, nussinov_matrix[i,k]+nussinov_matrix[k+1,j])
                    mk_values <- rbind(mk_values, data.frame(k = k,
                                                             val = nussinov_matrix[i,k]+nussinov_matrix[k+1,j]
                                                             )
                                       )

                  }
                  # select maximum value and k
                  mk <- max(mk_values$val, na.rm = TRUE)
                  kmax <- min(mk_values$k[mk_values$val == mk], na.rm = TRUE)
                } else {
                  mk <- NA
                  kmax <- NA
                }

                
                ### record maximum direction chosen for traceback ----
                if (is.na(mk)) {
                  max_dir_df <- data.frame(dir = c("md", "mi", "mj"),
                                           vals = c(md, mi, mj), 
                                           stringsAsFactors = FALSE)
                } else {
                  max_dir_df <- data.frame(dir = c("md", "mi", "mj", "mk"),
                                           vals = c(md, mi, mj, mk), 
                                           stringsAsFactors = FALSE)
                }
                
                
                ### randomly select direction ---
                if (randomize_direction==TRUE) {
                  max_dirs <- max_dir_df[which(max_dir_df$vals == max(max_dir_df$vals, na.rm = TRUE)),]
                  if (nrow(max_dirs)>1) {
                    max_dir_name <- sample(max_dirs$dir, size = 1)
                  } else {
                    max_dir_name <- max_dir_df[which.max(max_dir_df$vals),"dir"]
                  }
                 } else {
                  max_dir_name <- max_dir_df[which.max(max_dir_df$vals),"dir"]
                 }
                
                # set the value selected 
                max_dir_val <- max_dir_df$vals[max_dir_df$dir == max_dir_name]
                
                # update tracker
                if (!(is.null(max_dir_name))) {
                  if (max_dir_name == "mi") {
                    ref_posit <- paste(i + 1, j, sep = ",")
                    mat_dir <- "down"
                  } else if (max_dir_name == "mj") {
                    ref_posit <- paste(i, j - 1, sep = ",")
                    mat_dir <- "left"
                  } else if (max_dir_name == "md") {
                    ref_posit <- paste(i + 1, j - 1, sep = ",")
                    mat_dir <- "diagonal"
                  } else if (max_dir_name == "mk") {
                    #ref_posit <- paste(i, j, kmax, sep = ",")
                    ref_posit <- paste(paste(i, kmax, sep = ","),
                                       paste(kmax+1, j, sep = ",")
                                       , sep = "|")
                    mat_dir <- "bifurcation"
                  }
                  
                  # fill df
                  traceback_tracker[nrow(traceback_tracker) + 1,] <-
                    c(
                      cur_diag,
                      i,
                      j,
                      ij_posit = paste(i, j, sep = ","),
                      seq_split[i],
                      seq_split[j],
                      paste0(seq_split[i], seq_split[j]),
                      mi,
                      mj,
                      md,
                      mk,
                      max_dir_name,
                      ref_posit,
                      max_dir_val,
                      mat_dir
                    )
                  
                }
                
                ### update matrix ----
                nussinov_matrix[i, j] <- max_dir_val
                
              }
            }
          }
          
        }
      }
      cur_diag = cur_diag + 1
    }
    
    ### return output ----
    return(
      list(
        sequence = seq_split,
        nussinov_matrix = nussinov_matrix,
        traceback_tracker = traceback_tracker
      )
    )
  }

### function to traceback
Nussinov_traceback <-
  function(nussinov_traceback,
           nussinov_traceback_map = NULL) {
    ### find next row in traceback ----
    # logic to check if there's an existing map already
    # if not, create a map 
    if (is.null(nussinov_traceback_map)) {
      max_posit_j <- max(as.numeric(nussinov_traceback$j_posit))
      min_posit_i <- min(as.numeric(nussinov_traceback$i_posit))
      nussinov_traceback_map <-
        subset(
          nussinov_traceback,
          nussinov_traceback$ij_posit == paste(min_posit_i, max_posit_j, sep = ",")
        )
      last_ref_posit <-
        nussinov_traceback_map$ref_posit[nrow(nussinov_traceback_map)]
    # if there is a map, then select the last reference position
    } else {
      last_ref_posit <-
        nussinov_traceback_map$ref_posit[nrow(nussinov_traceback_map)]
      
      # logic to handle bifurcation 
      if (nussinov_traceback_map$mat_dir[nrow(nussinov_traceback_map)]=="bifurcation") {
        posits <- strsplit(last_ref_posit, "\\|")[[1]]
        
        next_ref_row <- data.frame()
        
        for (p in posits) {
          if (p %in% nussinov_traceback$ij_posit) {
            sub_ref_row <- subset(nussinov_traceback,
                                  nussinov_traceback$ij_posit == p)
            
            sub_ref_row <- Nussinov_traceback(nussinov_traceback, sub_ref_row)
            
            next_ref_row <- rbind(next_ref_row, sub_ref_row)
          }
        }
      } ### if not a bifurcated row, select the next
        else {
        next_ref_row <- subset(nussinov_traceback, nussinov_traceback$ij_posit == last_ref_posit)
      }
      
      # add the next ref to the map 
      nussinov_traceback_map <-
        rbind(nussinov_traceback_map, next_ref_row)
    }
    
    
    ### return output ----
    if (last_ref_posit %in% nussinov_traceback$ij_posit) {
      return(Nussinov_traceback(nussinov_traceback, nussinov_traceback_map))
    } else {
      return(nussinov_traceback_map)
    }
    
  }
