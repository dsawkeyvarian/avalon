cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(linac)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

#As far as I know there is no target defined in Geant4 package
add_library(Geant4 INTERFACE)
target_include_directories(Geant4
  INTERFACE
    ${Geant4_INCLUDE_DIR}
)
target_link_libraries(Geant4
  INTERFACE
    ${Geant4_LIBRARIES}
)

add_library(fmt STATIC src/external/fmt/format.cc)
target_include_directories(fmt PUBLIC src/external)

add_library(linac_core STATIC)
target_sources(linac_core
  PRIVATE
    src/logger.cpp
    src/machine_information.cpp
    src/linac.cpp 
    src/simulation.cpp 
    src/source.cpp
    src/parse.cpp
    src/plan.cpp
    src/truebeam.cpp 
    src/halcyon.cpp
    src/nohead.cpp
    src/math_utils.cpp
    src/stl.cpp
	src/avalon.cpp
  src/avalon_electron.cpp
  )
target_link_libraries(linac_core 
  PUBLIC
    Geant4
    fmt
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>
)

add_executable(linac src/linac_cli.cpp)
target_link_libraries(linac PUBLIC fmt linac_core)

add_executable(LinacTool)
target_sources(LinacTool
  PRIVATE
    src/linac_tool.cpp
)
target_link_libraries(LinacTool
  PRIVATE
    Geant4
    $<$<AND:$<CXX_COMPILER_ID:GNU>,$<VERSION_LESS:$<CXX_COMPILER_VERSION>,9.0>>:stdc++fs>)

add_executable(test)
target_sources(test
  PRIVATE
    src/test.cpp
    src/system_tests.cpp
    src/test_utils.cpp
)
target_link_libraries(test PUBLIC fmt linac_core)
set_target_properties(fmt linac_core linac test LinacTool PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS OFF)

install(TARGETS linac LinacTool test DESTINATION bin)
