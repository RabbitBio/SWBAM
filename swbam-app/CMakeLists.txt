cmake_minimum_required(VERSION 3.13)
project(SWBAMApp LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT TARGET swbam::swbam)
    find_package(SWBAM 0.1 CONFIG REQUIRED)
endif()
add_library(swbam_cpe_algorithm_objects OBJECT
    slave/sort_mpi.cpp slave/collate_mpi.cpp slave/markdup_mpi.cpp slave/fixmate_mpi.cpp)
target_link_libraries(swbam_cpe_algorithm_objects PRIVATE swbam::cpe)
target_include_directories(swbam_cpe_algorithm_objects PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(swbam_app_commands OBJECT
    CmdInfo.cpp
    commands/swbam_mpi_common.cpp
    commands/swbam_flagstat_mpi.cpp
    commands/swbam_bam2bam_mpi.cpp
    commands/swbam_bam2sam_mpi.cpp
    commands/swbam_sam2bam_mpi.cpp
    commands/swbam_fixmate_mpi.cpp
    algorithms/sort/swbam_sort_mpi.cpp
    algorithms/sort/swbam_sort_optimized_mpi.cpp
    algorithms/collate/swbam_collate_mpi.cpp
    algorithms/markdup/swbam_markdup_mpi.cpp
    algorithms/markdup/swbam_markdup_stream_mpi.cpp
    pipelines/swbam_dedup_workflow_mpi.cpp)
target_link_libraries(swbam_app_commands PRIVATE
    swbam::host swbam::mpi_platform swbam::htslib swbam::deflate)
target_compile_definitions(swbam_app_commands PRIVATE ENABLE_MPI_BAM2BAM)
target_compile_options(swbam_app_commands PRIVATE -Werror=format)
target_include_directories(swbam_app_commands PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/third_party)
add_executable(SWBAM main.cpp
    $<TARGET_OBJECTS:swbam_app_commands>
    $<TARGET_OBJECTS:swbam_cpe_algorithm_objects>)
target_compile_definitions(SWBAM PRIVATE ENABLE_MPI_BAM2BAM)
target_compile_options(SWBAM PRIVATE -Werror=format)
target_include_directories(SWBAM PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/third_party)
target_link_libraries(SWBAM PRIVATE swbam::swbam)
# The Sunway launcher requires this OS ABI marker on hybrid executables.
add_custom_command(TARGET SWBAM POST_BUILD
    COMMAND /bin/sh -c "printf '\\002' | dd of='$<TARGET_FILE:SWBAM>' bs=1 seek=9 count=1 conv=notrunc status=none"
    VERBATIM)
install(TARGETS SWBAM RUNTIME DESTINATION bin)
