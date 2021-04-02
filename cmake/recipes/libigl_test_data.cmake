## Test data
function(igl_download_test_data)
    igl_download_project_aux(test_data
        "${LIBIGL_EXTERNAL}/../tests/data"
        GIT_REPOSITORY https://github.com/libigl/libigl-tests-data
        GIT_TAG        19cedf96d70702d8b3a83eb27934780c542356fe
    )
endfunction()
