include_directories( ${SPF_INCLUDE_Dir})
target_sources(${PROJECT_NAME} PRIVATE
   ${CMAKE_CURRENT_SOURCE_DIR}/check_for_failure.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/readHDF5c.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/read_parameter_file.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/spf_3d.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/spf_communication.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/stochastic_rates.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/voxel_dynamics.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/writeHDF5c.cpp
)
