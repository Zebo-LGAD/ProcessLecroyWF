# Target Library

# Search all *.cpp files inside src
file(GLOB_RECURSE SOURCE_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_library(lcparser SHARED
    ${HEADER_FILES}
    ${SOURCE_FILES}
    # "${CMAKE_CURRENT_SOURCE_DIR}/src/lcparserUX.cxx"
)

# Set include directories
target_include_directories(lcparser PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Set target properties
set_target_properties(lcparser PROPERTIES
    OUTPUT_NAME "lcparser"
    POSITION_INDEPENDENT_CODE ON
)

target_link_libraries(lcparser
    PUBLIC ROOT::Core
    ROOT::Gui
)


# GENERATE DICTIONARIES

# Search all *.h files inside include, but exclude linkdef.h
file(GLOB_RECURSE HEADER_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
)
list(FILTER HEADER_FILES EXCLUDE REGEX "linkdef.h$")

ROOT_GENERATE_DICTIONARY(G__lcparser
    ${HEADER_FILES}
    LINKDEF "${CMAKE_CURRENT_SOURCE_DIR}/include/linkdef.h"
    MODULE lcparser
    # OPTIONS -f "lcparser/src/lcparserUX.cxx"
)
