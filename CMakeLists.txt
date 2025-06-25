cmake_minimum_required(VERSION 3.10)
project(CrossFlow VERSION 1.0.0 LANGUAGES C CXX)

# Set global C++ settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release)

# Find dependencies
find_package(MPI COMPONENTS C CXX REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(VTK 9.0 REQUIRED
    CommonCore
    CommonDataModel
    IOLegacy
    IOPLY
    IOXML
    InteractionStyle
    RenderingFreeType
    RenderingOpenGL2
    IOGeometry
    IOExport
    RenderingCore
    RenderingAnnotation
    RenderingVolume
    RenderingVolumeOpenGL2
    jsoncpp
)

# Create interface library for headers
file(GLOB_RECURSE CROSSFLOW_HEADERS "src/*.hpp")
add_library(crossflow INTERFACE)
target_sources(crossflow INTERFACE ${CROSSFLOW_HEADERS})

# Set includes and links
target_include_directories(crossflow INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${EIGEN3_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
)

target_link_libraries(crossflow INTERFACE
    Eigen3::Eigen
    ${VTK_LIBRARIES}
)

target_compile_definitions(crossflow INTERFACE
    EIGEN_NO_DEBUG=1
)

# Add examples if present
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/examples/CMakeLists.txt)
    add_subdirectory(examples)
endif()

enable_testing()
find_package(GTest REQUIRED)
file(GLOB TEST_SOURCES "tests/*.cpp")
foreach(test_source ${TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)
    add_executable(${test_name} ${test_source})
    target_link_libraries(${test_name} PRIVATE crossflow GTest::gtest GTest::gtest_main)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
