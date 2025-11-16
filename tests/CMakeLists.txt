
add_executable(test_wfdata test_wfdata.cpp)
target_link_libraries(test_wfdata PUBLIC lcparser ROOT::Core ROOT::Tree ROOT::Graf)
add_test(NAME test_wfdata COMMAND test_wfdata)

add_executable(process_sr90 process_sr90.cpp)
target_link_libraries(process_sr90 PUBLIC lcparser ROOT::Core ROOT::Tree ROOT::Graf)
add_test(NAME process_sr90 COMMAND process_sr90)

add_executable(test_analyze test_analyze.cpp)
target_link_libraries(test_analyze PUBLIC analyze lcparser ROOT::Core ROOT::Tree ROOT::Graf)
add_test(NAME test_analyze COMMAND test_analyze)
