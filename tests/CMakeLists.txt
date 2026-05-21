find_package(GTest REQUIRED)
include(GoogleTest)

file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/urqmd_input ${CMAKE_CURRENT_BINARY_DIR}/urqmd_input)
file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/tables.dat ${CMAKE_CURRENT_BINARY_DIR}/tables.dat)

add_executable(UrQMD_Test pdg_quantum_numbers.cc test.cc)

set_target_properties(UrQMD_Test PROPERTIES CXX_CLANG_TIDY "${COLA_CMAKE_CXX_CLANG_TIDY}")

target_link_libraries(UrQMD_Test PRIVATE COLA_UrQMD GTest::Main GTest::GTest)

gtest_discover_tests(UrQMD_Test)
