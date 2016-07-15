# compile everything except the main source file into a statis lib to be linked 
# to the main executable
add_library(${LIBRARY_NAME} STATIC
  scattering_2to2.cpp
  scattering_2to3.cpp
  sample_2to2.cpp
)

set_target_properties(${LIBRARY_NAME} PROPERTIES PREFIX "")


# compile the actual executable
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${LIBRARY_NAME} ${Boost_LIBRARY} ${GSL_LIBRARIES} ${GSLCALAS_LIBRARIES})

install(TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_NAME})
