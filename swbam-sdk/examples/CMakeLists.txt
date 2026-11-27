cmake_minimum_required(VERSION 3.13)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(SWBAMExamples LANGUAGES C CXX)
    find_package(SWBAM 0.1 CONFIG REQUIRED)
endif()
add_executable(swbam-sdk-count sdk_count.cpp)
add_executable(swbam-sdk-bam2bam sdk_bam2bam.cpp)
foreach(target swbam-sdk-count swbam-sdk-bam2bam)
    target_link_libraries(${target} PRIVATE swbam::swbam)
    add_custom_command(TARGET ${target} POST_BUILD
        COMMAND /bin/sh -c "printf '\\002' | dd of='$<TARGET_FILE:${target}>' bs=1 seek=9 count=1 conv=notrunc status=none"
        VERBATIM)
endforeach()
