# bench utils library
add_library(bench_utils STATIC tools/bench_utils.c tools/bench_timing.c)

foreach(target IN ITEMS F_2_63 hashMiMC bench)
    add_executable("${target}_test" "${target}.cpp")
    target_link_libraries("${target}_test" HashMIMC bench_utils)
    add_test(NAME ${target} COMMAND "${target}_test")
endforeach()
