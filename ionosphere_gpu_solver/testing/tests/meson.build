
####################################################################################################
# Tests
####################################################################################################

# Structure is that each element is: name_of_test : [source_files_for_tests]
tests = {}

tests += {
    'vector_norm_squared_1' : [files('vector_norm_squared_1.cpp')],
    'vector_norm_squared_2' : [files('vector_norm_squared_2.cpp')],
    'vector_addition_1' : [files('vector_addition_1.cpp')],
    'vector_subtraction_1' : [files('vector_subtraction_1.cpp')],
    'vector_elementwise_multiplication_1' : [files('vector_elementwise_multiplication_1.cpp')],
    'vector_elementwise_division_1' : [files('vector_elementwise_division_1.cpp')],
    'vector_dot_product_1' : [files('vector_dot_product_1.cpp')],
    'muliply_vector_with_scalar_and_add_it_to_another_vector_1' : [files('muliply_vector_with_scalar_and_add_it_to_another_vector_1.cpp')],
    'transpose_matrix_1' : [files('transpose_matrix_1.cpp')],
    'create_sparse_matrix_1' : [files('create_sparse_matrix_1.cpp')],
    'create_sparse_matrix_2' : [files('create_sparse_matrix_2.cpp')],
    'create_sparse_matrix_3' : [files('create_sparse_matrix_3.cpp')],
    'matrix_vector_product_cuda_1' : [files('matrix_vector_product_cuda_1.cpp')],
    'sparse_matrix_vector_product_cuda_1' : [files('sparse_matrix_vector_product_cuda_1.cpp')],
    'sparse_matrix_vector_product_cuda_2' : [files('sparse_matrix_vector_product_cuda_2.cpp')],
    'sparse_BiGC_1' : [files('sparse_BiGC_1.cpp')],
    'sparse_BiGC_2' : [files('sparse_BiGC_2.cpp')],
    'sparse_BiGCSTAB_1' : [files('sparse_BiGCSTAB_1.cpp')],
    'sparse_BiGCSTAB_2' : [files('sparse_BiGCSTAB_2.cpp')]
}

test_names_and_exes = {}
foreach test_name, test_sources : tests
    test_names_and_exes += {test_name : executable(
        test_name,
        test_sources,
        include_directories: project_include_directories,
        dependencies: dependencies)
    }
    
endforeach

####################################################################################################
# Running tests
####################################################################################################

foreach test_name, test_exe : test_names_and_exes
    test(test_name, test_exe, is_parallel: false)
endforeach
