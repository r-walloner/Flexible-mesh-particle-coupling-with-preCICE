cmake_minimum_required(VERSION 3.13.4)

set(TARGET "step-68-with-precice")

find_package(deal.II 9.6.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
deal_ii_initialize_cached_variables()

project(${TARGET})
add_executable(${TARGET} ${TARGET}.cc)

deal_ii_setup_target(${TARGET})
deal_ii_initialize_cached_variables()

find_package(precice REQUIRED CONFIG)
message(STATUS "Using the preCICE version found at ${precice_DIR}")
target_link_libraries(${TARGET} precice::precice)
