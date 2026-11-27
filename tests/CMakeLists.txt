add_executable(swbam-smoke-data smoke_data.cpp)
target_link_libraries(swbam-smoke-data PRIVATE
    swbam::host swbam::mpi_platform swbam::htslib)
target_link_options(swbam-smoke-data PRIVATE -mhybrid -Wl,-zmuldefs)
add_custom_command(TARGET swbam-smoke-data POST_BUILD
    COMMAND /bin/sh -c "printf '\\002' | dd of='$<TARGET_FILE:swbam-smoke-data>' bs=1 seek=9 count=1 conv=notrunc status=none"
    VERBATIM)
