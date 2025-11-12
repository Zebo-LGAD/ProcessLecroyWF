# Waveform data processor

## Aims:
- Process waveform data from lecroy oscilloscope, usually ends with ".trc" extensions.
- [optional] Converte waveform data to ROOT file.
- [optional] Extract waveform infomation, and fill into ROOT file.
- [TODO] Give some algorithm to analyze temporal & spatial resolution of sensors
- [TODO] Add Qt5 GUI for choosing folders to be processed.

## Usage:
- This program is organized by CMakeLists.txt, you can add your own modules easily.
  - Just create a few directories, i.e. "userModule/", "userModule/src/", "userModule/include".
  - Add a CMakeLists.txt to your module directory "userModule/"

    ``` cmake
    file(GLOB_RECURSE HEADER_FILES "{CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
    file(GLOB_RECURSE SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
    
    # Generate a library named as "libuserModule.so", which can be find using "userModule" in CMake
    add_library(userModule ${HEADER_FILES} ${SOURCE_FILES})
    
    target_include_directories(userModule PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> 
        $<INSTALL_INTERFACE:include>)
    target_link_libraries(userModule ROOT::Core ...)
    ...
    ```
  - If you want to add an executable program (`run` for e.g.), add cmake codes and main.cpp at directory "userModule/"
    ``` cmake
    add_executable (run main.cpp)
    target_link_libraries(run PRIVATE
        userModule ROOT:Core ...
    )
    install(TARGETS run 
        DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
    ```
  - If you want to add a sub directory and sub module, just add cmake codes:
    ``` cmake
        add_subdirectory(sub1)
        add_subdirectory(sub2)
    ```

- If you just want to add your own program quickly, just add some cpp file in "tests/" directory, and add these codes in "tests/CMakeLists.txt".
    ``` cmake
    add_executable(your_program your_program.cpp)
    target_link_libraries(your_program PUBLIC lcparser ROOT::Core ROOT::Tree ROOT::Graf ...)
    add_test(NAME your_program COMMAND your_program)
    ```

## Module introduction:
### lcparser
- Process raw data from lecroy oscilloscope.
- lcparser.h & lcparser.cpp provide a method to read raw ".trc" file
- WFDataConverter.h provide a method to convert to ROOT file, or just extract useful information directly.

## Problems not saved yet:
- May have problem in the first configuration, just time "cmake ." again in build directory. This may due to unset "CMAKE_INSTALL_PREFIX" in "ROOT_GENERATE_DICTIONARY" provieded by ROOT package.