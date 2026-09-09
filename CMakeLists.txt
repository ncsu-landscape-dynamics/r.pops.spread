cmake_minimum_required(VERSION 3.22)

project(r.pops.spread)

include(build_addon)

find_M()
find_OpenMP()

# Steering runs the TCP client in a separate thread.
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(sources
  graster.hpp
  tcp_client.h
  tcp_client.cpp
  main.cpp)

build_addon(
  ${PROJECT_NAME}
  SOURCES ${sources}
  GRASSLIBS
    datetime
    gis
    raster
    vector
  DEPENDS
    LIBM
  OPTIONAL_DEPENDS
    OpenMP::OpenMP_CXX
  DOCFILES
    pops_logo.png
    r_pops_spread.png
    ${PROJECT_NAME}.html
    ${PROJECT_NAME}.md)

# build_addon() uses the plain target_link_libraries signature, so match it.
target_link_libraries(${PROJECT_NAME} Threads::Threads)

target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include>)

target_sources(${PROJECT_NAME} PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/config.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/date.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/deterministic_kernel.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/host_pool.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/kernel_types.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/kernel.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/model.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/model_type.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/multi_host_pool.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/natural_anthropogenic_kernel.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/neighbor_kernel.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/pest_pool.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/quarantine.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/radial_kernel.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/raster.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/scheduling.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/simulation.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/spread_rate.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/statistics.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/switch_kernel.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/treatments.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/uniform_kernel.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/pops-core/include/pops/utils.hpp)
