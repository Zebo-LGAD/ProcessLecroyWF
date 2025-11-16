
# Search all *.h files inside include, but exclude linkdef.h
file(GLOB_RECURSE HEADER_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
)
list(FILTER HEADER_FILES EXCLUDE REGEX "linkdef.h$")

# Search all *.cpp files inside src
file(GLOB_RECURSE SOURCE_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

# If not found any source files, just skip creating the library
if(NOT SOURCE_FILES)
    message(WARNING "No source files found in ${CMAKE_CURRENT_SOURCE_DIR}/src. Skipping analyze library creation.")
    return()
endif()

add_library(analyze SHARED
    ${HEADER_FILES}
    ${SOURCE_FILES}
)

# Set include directories
target_include_directories(analyze PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Set target properties
set_target_properties(analyze PROPERTIES
    OUTPUT_NAME "analyze"
    POSITION_INDEPENDENT_CODE ON
)

target_link_libraries(analyze
    PUBLIC lcparser
    PUBLIC ROOT::Core ROOT::Matrix ROOT::MathCore
)
