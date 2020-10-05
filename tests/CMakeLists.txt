set(TEST_FILES
	main_tests.cpp
	const_function.cpp
	linear_function.cpp
	quadratic_function.cpp
	trig_function.cpp
)

enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(verlet-test ${TEST_FILES})
target_link_libraries(verlet-test ${GTEST_LIBRARIES} pthread)
GTEST_ADD_TESTS(verlet-test "" ${TEST_FILES})
message(STATUS "CXX flags : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
