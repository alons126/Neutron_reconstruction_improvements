cmake_minimum_required(VERSION 3.16)
project(Neutron_reconstruction_improvements LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES
    "/Users/alon/Projects/Neutron_reconstruction_improvements/deepn_noBanks.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/D_getfeatures_ppim.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Older_code_phases/D_getfeatures_Phase3 copy 2.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Older_code_phases/ManualVeto_Phase6.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Older_code_phases/D_getfeatures_Phase6 copy 2.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Older_code_phases/D_getfeatures_Phase6 copy.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Older_code_phases/D_getfeatures_Phase3.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Older_code_phases/D_getfeatures_Phase2.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Older_code_phases/D_getfeatures_Phase1.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Older_code_phases/D_getfeatures_Phase5.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Older_code_phases/D_getfeatures_Phase2 copy.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Older_code_phases/D_getfeatures_Phase3 copy.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Older_code_phases/D_getfeatures_Phase4.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Older_code_phases/D_getfeatures_Phase5 copy.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/ManualVeto_Phase7.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Other/cmake-build-debug/CMakeFiles/3.29.6/CompilerIdCXX/CMakeCXXCompilerId.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Other/Erin_features_original/D_getfeatures_ppim.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Other/Erin_features_original/D_getfeatures.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Other/Erin_features_original/N_getfeatures.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Other/Erin_features_original/veto_functions.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Other/Andrew_manual_cuts_original/deepn_noBanks.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/Other/Andrew_manual_cuts_original/deepn_withVeto.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/N_getfeatures.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/build/CMakeFiles/3.31.2/CompilerIdCXX/CMakeCXXCompilerId.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/deepn_withVeto.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/main.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/src/classes/clas12ana/clas12ana.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/src/classes/clas12debug/clas12debug.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/src/functions/Andrews_functions/Andrews_functions.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/src/functions/neutron-veto/veto_functions.cpp"
    "/Users/alon/Projects/Neutron_reconstruction_improvements/src/functions/HipoChain_config.cpp"
)

add_executable(${PROJECT_NAME} ${SOURCES})
