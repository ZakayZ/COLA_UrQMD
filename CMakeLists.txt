cmake_minimum_required(VERSION 3.14)
project(COLA_UrQMD VERSION 0.1.0 LANGUAGES CXX Fortran)

include(GNUInstallDirs)

find_package(COLA_Fortran REQUIRED)

set(CMAKE_INSTALL_PREFIX ${COLA_DIR})
set(COLA_MODULE_NAME COLA_UrQMD)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(URQMD_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/urqmd-4.0)

set(URQMD_F_SOURCES
  ${URQMD_SOURCE_DIR}/1fluid.f
  ${URQMD_SOURCE_DIR}/addpart.f
  ${URQMD_SOURCE_DIR}/angdis.f
  ${URQMD_SOURCE_DIR}/anndec.f
  ${URQMD_SOURCE_DIR}/bessel.f
  ${URQMD_SOURCE_DIR}/blockres.f
  ${URQMD_SOURCE_DIR}/boxprg.f
  ${URQMD_SOURCE_DIR}/cascinit.f
  ${URQMD_SOURCE_DIR}/CFmax.f90
  ${URQMD_SOURCE_DIR}/coload.f
  ${URQMD_SOURCE_DIR}/cornelius.f90
  ${URQMD_SOURCE_DIR}/dectim.f
  ${URQMD_SOURCE_DIR}/delpart.f
  ${URQMD_SOURCE_DIR}/detbal.f
  ${URQMD_SOURCE_DIR}/dwidth.f
  ${URQMD_SOURCE_DIR}/error.f
  ${URQMD_SOURCE_DIR}/getmass.f
  ${URQMD_SOURCE_DIR}/getspin.f
  ${URQMD_SOURCE_DIR}/gnuranf.f
  ${URQMD_SOURCE_DIR}/hepchg.f
  ${URQMD_SOURCE_DIR}/hepcmp.f
  ${URQMD_SOURCE_DIR}/hepnam.f
  ${URQMD_SOURCE_DIR}/init.f
  ${URQMD_SOURCE_DIR}/input.f
  ${URQMD_SOURCE_DIR}/iso.f
  ${URQMD_SOURCE_DIR}/ityp2pdg.f
  ${URQMD_SOURCE_DIR}/jdecay2.f
  ${URQMD_SOURCE_DIR}/make22.f
  ${URQMD_SOURCE_DIR}/numrec.f
  ${URQMD_SOURCE_DIR}/output.f
  ${URQMD_SOURCE_DIR}/paulibl.f
  ${URQMD_SOURCE_DIR}/proppot.f
  ${URQMD_SOURCE_DIR}/pythia6409.f
  ${URQMD_SOURCE_DIR}/quadri.f90
  ${URQMD_SOURCE_DIR}/saveinfo.f
  ${URQMD_SOURCE_DIR}/scatter.f
  ${URQMD_SOURCE_DIR}/siglookup.f
  ${URQMD_SOURCE_DIR}/string.f
  ${URQMD_SOURCE_DIR}/tabinit.f
  ${URQMD_SOURCE_DIR}/uhmerge.f
  ${URQMD_SOURCE_DIR}/upmerge.f
  ${URQMD_SOURCE_DIR}/whichres.f
)

add_cola_fortran_library(
  ${COLA_MODULE_NAME}
  SOURCES
    urqmd_generator.f90
    urqmd_ref.f
    ${URQMD_F_SOURCES}
)

set_source_files_properties(
  urqmd_ref.f
  PROPERTIES COMPILE_OPTIONS "-ffixed-form;-ffixed-line-length-120"
)

target_compile_options(${COLA_MODULE_NAME} PRIVATE -fPIC)

target_include_directories(${COLA_MODULE_NAME} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${URQMD_SOURCE_DIR}>
)

target_include_directories(${COLA_MODULE_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_include_directories(${COLA_MODULE_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${COLA_MODULE_NAME}/modules>
)

# ---------------------------------------------------------------------------
# Installation
# ---------------------------------------------------------------------------

set(GENERATED_HEADERS_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated/${COLA_MODULE_NAME})

install(TARGETS ${COLA_MODULE_NAME}
    EXPORT ${COLA_MODULE_NAME}Export
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/${COLA_MODULE_NAME}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${COLA_MODULE_NAME}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${COLA_MODULE_NAME}
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/data/${COLA_MODULE_NAME}Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${COLA_MODULE_NAME}
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}ConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(EXPORT ${COLA_MODULE_NAME}Export DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${COLA_MODULE_NAME})

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${COLA_MODULE_NAME}/modules
    FILES_MATCHING PATTERN "*.mod"
)

install(DIRECTORY ${GENERATED_HEADERS_DIR}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${COLA_MODULE_NAME}
    FILES_MATCHING PATTERN "*.hh"
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${COLA_MODULE_NAME}ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${COLA_MODULE_NAME}
)

enable_testing()
if (BUILD_TESTING)
    add_subdirectory(tests)
endif ()
