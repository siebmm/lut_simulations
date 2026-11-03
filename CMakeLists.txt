cmake_minimum_required(VERSION 3.21)
project(LaserWeldSAFT LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Qt6 REQUIRED COMPONENTS Core)

set(NATIVE saft/src)
set(COMMON ${NATIVE}/saft.cpp ${NATIVE}/screening.cpp
    ${NATIVE}/dataset_geometry.cpp ${NATIVE}/groove_paths.cpp)
add_executable(SAFTFilterStudy saft/cli/filter_cli.cpp ${COMMON})
add_executable(SAFTReview saft/cli/review_cli.cpp ${COMMON})
add_executable(SAFTAdmissionAudit saft/cli/admission_cli.cpp ${COMMON})
add_executable(SAFTBatch saft/cli/batch_cli.cpp ${COMMON})
foreach(target SAFTFilterStudy SAFTReview SAFTAdmissionAudit SAFTBatch)
  target_include_directories(${target} PRIVATE ${NATIVE})
  target_link_libraries(${target} PRIVATE Qt6::Core)
  target_compile_definitions(${target} PRIVATE QT_NO_KEYWORDS)
  set_target_properties(${target} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endforeach()
