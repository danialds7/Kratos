# This function automatically configures a given application to build its tests
macro(kratos_add_gtests application_core_name application_test_sources)
    include(GoogleTest)
    add_executable(${application_core_name}Test ${application_test_sources})
    target_link_libraries(${application_core_name}Test ${application_core_name} KratosCoreTestUtilities GTest::gtest_main GTest::gmock_main)
    set_target_properties(${application_core_name}Test PROPERTIES COMPILE_DEFINITIONS "KRATOS_TEST_CORE=IMPORT,API")
    install(TARGETS ${application_core_name}Test DESTINATION test)
    gtest_discover_tests(${application_core_name}Test DISCOVERY_MODE PRE_TEST)
endmacro(kratos_add_gtests)