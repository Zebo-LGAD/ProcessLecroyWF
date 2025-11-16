cmake_minimum_required(VERSION 3.1)
project(ProcessLecroyWF
    VERSION 1.0
    LANGUAGES CXX
    DESCRIPTION "Process Lecroy Waveform Files"
)

find_package(ROOT REQUIRED)

# Set default build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Define CMAKE_MODULE_PATH to include custom modules
# list(APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
# include (StandardProjectSettings)
# include (CompilerSettings)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add subdirectories
# add_subdirectory(lcparser)
add_subdirectory(libs)
add_subdirectory(tests)
add_subdirectory(apps)
# add_subdirectory(analyze)

# Installation configuration
include (GNUInstallDirs)

# Packaging configuration
install (
    TARGETS lcparser 
    EXPORT ProcessLecroyWFTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install (
    DIRECTORY libs/lcparser/include/ 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Exporting targets for CMake find_package
install (EXPORT ProcessLecroyWFTargets
    FILE ProcessLecroyWFTargets.cmake
    NAMESPACE ProcessLecroyWF::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ProcessLecroyWF
)


# Link rootlogon.C to build folder
configure_file(${CMAKE_SOURCE_DIR}/rootlogon.C ${CMAKE_BINARY_DIR}/rootlogon.C COPYONLY)
install (FILES ${CMAKE_BINARY_DIR}/rootlogon.C DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})