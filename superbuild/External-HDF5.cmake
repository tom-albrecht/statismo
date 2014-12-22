message( "External project - HDF5" )

message(STATUS "CMAKE_DEBUG_POSTFIX  = " ${CMAKE_DEBUG_POSTFIX}) 


ExternalProject_add(HDF5
  SOURCE_DIR ${CMAKE_BINARY_DIR}/HDF5
  BINARY_DIR ${CMAKE_BINARY_DIR}/HDF5-build
  URL http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.13/src/hdf5-1.8.13.tar.gz
  UPDATE_COMMAND ""
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}
    -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=OFF
    -DHDF5_BUILD_CPP_LIB:BOOL=ON
    -DBUILD_SHARED_LIBS:BOOL=ON
    -DHDF5_BUILD_TOOLS:BOOL=OFF
    #-DCMAKE_DEBUG_POSTFIX:STRING=""
    -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DEPENDENCIES_DIR}
  INSTALL_DIR ${INSTALL_DEPENDENCIES_DIR}
)

set( ENV{HDF5_ROOT_DIR_HINT} )

# We call find_package once so the variables are available
find_package(HDF5)
# translate between the two versions of HDF5
set(HDF5_DIR ${HDF_ROOT_DIR})
get_filename_component(HDF5_LIBRARY_DIRS hdf5 DIRECTORY)
get_target_property(HDF5_C_LIBRARY hdf5 IMPORTED_IMPLIB_RELEASE)
get_target_property(HDF5_CXX_LIBRARY hdf5_cpp IMPORTED_IMPLIB_RELEASE)

message(STATUS "HDF5_C_LIBRARY = " ${HDF5_C_LIBRARY}) 
message(STATUS "BUILD_TYPE = " ${CMAKE_CFG_INTDIR})

#if (WIN32)
#  # On Windows, find_package(HDF5) with cmake 2.8.[8,9] always ends up finding
#  # the dlls instead of the libs. So setting the variables explicitly for
#  # dependent projects.
#  set(cmake_hdf5_c_lib    -DHDF5_C_LIBRARY:FILEPATH=${INSTALL_DEPENDENCIES_DIR}/lib/hdf5.lib)
#  set(cmake_hdf5_cxx_lib  -DHDF5_CXX_LIBRARY:FILEPATH=${INSTALL_DEPENDENCIES_DIR}/lib/hdf5_cpp.lib)
#  set(cmake_hdf5_libs     ${cmake_hdf5_c_lib} ${cmake_hdf5_cxx_lib})
#endif ()

