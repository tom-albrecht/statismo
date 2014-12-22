message( "External project - VTK" )

find_package(Git)
if(NOT GIT_FOUND)
  message(ERROR "Cannot find git. git is required for Superbuild")
endif()

option( USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)

set(git_protocol "git")
if(NOT USE_GIT_PROTOCOL)
  set(git_protocol "http")
endif()

set( VTK_DEPENDENCIES )

if( ${USE_SYSTEM_HDF5} MATCHES "OFF" )
  set( VTK_DEPENDENCIES HDF5 )
endif()

set( _vtkOptions )
if( APPLE )
  set( _vtkOptions -DVTK_REQUIRED_OBJCXX_FLAGS:STRING="" )
endif()

MESSAGE(STATUS "cmake module path in external VTK = " ${CMAKE_MODULE_PATH})
MESSAGE(STATUS "cmake hdf5 root_dir = " ${HDF5_ROOT_DIR})
MESSAGE(STATUS "cmake hdf5 dir = " ${HDF5_DIR})
MESSAGE(STATUS "hdf5 c library = " ${HDF5_C_LIBRARY})
MESSAGE(STATUS "hdf5 cxx library = " ${HDF5_CXX_LIBRARY})
MESSAGE(STATUS "hdf5 library dirs = " ${HDF5_LIBRARY_DIRS})


ExternalProject_Add(VTK
  DEPENDS ${VTK_DEPENDENCIES}
  GIT_REPOSITORY ${git_protocol}://vtk.org/VTK.git
  GIT_TAG v6.1.0
  SOURCE_DIR VTK
  BINARY_DIR VTK-build
  UPDATE_COMMAND ""
  PATCH_COMMAND ""
  CMAKE_GENERATOR ${EP_CMAKE_GENERATOR}
  CMAKE_ARGS
    ${ep_common_args}
    ${_vtkOptions}
    -DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}  
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_SHARED_LIBS:BOOL=ON
    -DBUILD_TESTING:BOOL=OFF
    -DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}
    -DVTK_BUILD_ALL_MODULES:BOOL=OFF
    -DVTK_USE_SYSTEM_HDF5:BOOL=ON
    -DHDF5_DIR:PATH=${HDF5_ROOT_DIR}    
    -DHDF5_C_LIBRARY:PATH=${HDF5_C_LIBRARY}    
    -DHDF5_CXX_LIBRARY:PATH=${HDF5_CXX_LIBRARY}
    #-DHDF5_LIBRARY_DIRS:PATH=${HDF5_LIBRARY_DIRS}        
    -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DEPENDENCIES_DIR}
)

set( VTK_DIR ${INSTALL_DEPENDENCIES_DIR}/lib/cmake/vtk-6.1/ )
