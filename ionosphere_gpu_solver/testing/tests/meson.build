
####################################################################################################
# Tests
####################################################################################################

# Structure is that each element is: name_of_test : [source_files_for_tests]
tests = {}

tests += {
    'create_sparse_matrix_1' : [files('create_sparse_matrix_1.cpp')],
    'create_sparse_matrix_2' : [files('create_sparse_matrix_2.cpp')],
    'matrix_vector_product_cuda_1' : [files('matrix_vector_product_cuda_1.cpp')],
    'sparse_matrix_vector_product_cuda_1' : [files('sparse_matrix_vector_product_cuda_1.cpp')],
    'sparse_matrix_vector_product_cuda_2' : [files('sparse_matrix_vector_product_cuda_2.cpp')]
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
    test(test_name, test_exe)
endforeach
