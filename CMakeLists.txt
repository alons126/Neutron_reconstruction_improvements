cmake_minimum_required(VERSION 3.16)
project(Neutron_reconstruction_improvements)

set(CMAKE_CXX_STANDARD 17)

find_package(ROOT REQUIRED)

include_directories(${ROOT_INCLUDE_DIRS})

set(CLAS12ROOT_INCLUDE_DIR /home/alon/clas12root/Clas12Root)
set(CLAS12BANKS_INCLUDE_DIR /home/alon/clas12root/Clas12Banks)
set(HIPO4_INCLUDE_DIR /home/alon/clas12root/hipo4)

include_directories(${CLAS12ROOT_INCLUDE_DIR})
include_directories(${CLAS12BANKS_INCLUDE_DIR})
include_directories(${HIPO4_INCLUDE_DIR})

set(CLAS12ROOT_LIBRARY /home/alon/clas12root/lib/libClas12Root.so)
set(CLAS12BANKS_LIBRARY /home/alon/clas12root/lib/libClas12Banks.so)
set(HIPO4_LIBRARY /home/alon/clas12root/lib/libHipo4.so)

file(GLOB CLAS12ROOT_SOURCES "${CLAS12ROOT_INCLUDE_DIR}/*.cpp")
file(GLOB CLAS12BANKS_SOURCES "${CLAS12BANKS_INCLUDE_DIR}/*.cpp")
file(GLOB HIPO4_SOURCES "${HIPO4_INCLUDE_DIR}/*.cpp")

#add_executable(Neutron_reconstruction_improvements ${CLAS12ROOT_SOURCES} ${CLAS12BANKS_SOURCES} ${HIPO4_SOURCES})

add_executable(Neutron_reconstruction_improvements
        main.cpp
        TMVAClassification.C
        TMVAClassificationApplication.C
        TMVARegression.C
        Original/TMVAClassification.C
        Original/TMVAClassificationApplication.C
        Original/TMVARegression.C)

target_link_libraries(Neutron_reconstruction_improvements ${ROOT_LIBRARIES} ${CLAS12ROOT_LIBRARY} ${CLAS12BANKS_LIBRARY} ${HIPO4_LIBRARY})