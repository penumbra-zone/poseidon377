
set(LIB_SRCS
    HadesMIMC.cpp
    )

add_library(HashMIMC STATIC ${LIB_SRCS})
target_include_directories(HashMIMC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
