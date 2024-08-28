cmake_minimum_required(VERSION 3.16)
project(Neutron_reconstruction_improvements)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED)

include_directories(${ROOT_INCLUDE_DIRS})

add_executable(Neutron_reconstruction_improvements
        main.cpp
        TMVAClassification.C
        TMVAClassificationApplication.C
        TMVARegression.C
        Li_macros/trainBDT.C
        Li_macros/applyBDT.C
        Li_macros/Original/trainBDT.C
        Li_macros/Original/applyBDT.C
        TMVAClassification_original.C
)

target_link_libraries(Neutron_reconstruction_improvements ${ROOT_LIBRARIES})