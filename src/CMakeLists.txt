# Generate placeholder "Nflexpd_version" module (versions all set to "unknown" - updated at build time if FABM_EMBED_VERSION is set)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.F90.in ${CMAKE_CURRENT_BINARY_DIR}/version.F90)

add_library(fabm_models_nflexpd OBJECT
            nflexpd_model_library.F90
            ${CMAKE_CURRENT_BINARY_DIR}/version.F90
            common.F90
            abio.F90
            phy_classic.F90
            phy_DOQ.F90
            phy_IOQ.F90
            lambert_modified.f90
           )

add_dependencies(fabm_models_nflexpd fabm_base)

if(FABM_EMBED_VERSION)
  # Dynamically generate "nflexpd_version" module with git commit/branch information at every build.
  add_custom_target(fabm_models_nflexpd_version
     ${CMAKE_COMMAND} -DINFILE=${CMAKE_CURRENT_SOURCE_DIR}/version.F90.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/version.F90 -P ${GET_GIT_INFO_SCRIPT}
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
     COMMENT "Retrieving description of last nflexpd commit..." 
     SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/version.F90.in ${GET_GIT_INFO_SCRIPT}
     VERBATIM
  )
  add_dependencies(fabm_models_nflexpd fabm_models_nflexpd_version)
endif()

