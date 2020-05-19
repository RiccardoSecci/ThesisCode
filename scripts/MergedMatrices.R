finalmatrixcreate <- function(input_matrixup,input_matrixdown){
	
	if ( all(colnames(input_matrixup) ==colnames( input_matrixdown)) &
		 all(colnames(input_matrixup) ==colnames( input_matrixdown))){


	final_matrix <- matrix(nrow = nrow(input_matrixup),ncol = ncol(input_matrixup), dimnames = dimnames(input_matrixup) )

	for(rowindex in 1:nrow(final_matrix)){
		print(rowindex)
		row_up = input_matrixup[rowindex,]
		row_down = input_matrixdown[rowindex,]
		index_sign_down = which(row_down != row_up & row_down == 0.5)
		#index_non_sign_down = which(row_up == 2 & row_down == 1)
		result_addition= row_up + row_down
		#result_addition[index_non_sign_down ] = -result_addition[index_non_sign_down ]
		if (length(index_sign_down)!=0 ){
			result_addition[index_sign_down] = -result_addition[index_sign_down]
		}
		final_matrix[rowindex,] <- as.numeric(result_addition)

		}
		
		final_matrix[final_matrix == 2] = 3
		final_matrix[final_matrix == 1 ] = 0
		
		#somma = rowSums(final_matrix)
		#index_val = which(somma ==3*ncol(final_matrix))
		#if (length(index_val)>0){
		#	cat( 'Number of row removed', length(index_val),'\n')
		#	final_matrix = final_matrix[-index_val,]
		#}
		#	somma = rowSums(final_matrix)
		#index_val = which(somma ==4*ncol(final_matrix))
		#if (length(index_val)>0){
	#		cat( 'Number of row removed', length(index_val),'\n')
	#		final_matrix = final_matrix[-index_val,]
	#	}
		} else {
		print('Matrix up and down merging failed: different column or row ')
		}
	return(final_matrix)
}
