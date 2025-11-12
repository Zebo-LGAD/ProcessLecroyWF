
add_executable(test_wfdata test_wfdata.cpp)

target_link_libraries(test_wfdata PUBLIC lcparser ROOT::Core ROOT::Tree ROOT::Graf)
add_test(NAME test_wfdata COMMAND test_wfdata)