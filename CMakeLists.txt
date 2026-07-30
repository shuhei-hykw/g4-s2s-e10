#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.17)
cmake_host_system_information(RESULT host QUERY HOSTNAME)
cmake_host_system_information(RESULT ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
project(S2S)

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
#find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
option(WITH_GEANT4_UIVIS "Build HypTPC with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# macOS: Geant4Config pulls in EXPAT, and CMake's FindEXPAT resolves it
# to the Xcode SDK.  The resulting -isystem <SDK>/usr/include is placed
# ahead of libc++'s own math.h wrapper, so <cmath> aborts with
# "didn't find libc++'s <math.h> header".  That directory is the default
# sysroot include anyway, so dropping it from the imported target costs
# nothing and keeps the library itself linked.
#
if(APPLE AND TARGET EXPAT::EXPAT)
  set_target_properties(EXPAT::EXPAT PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "")
endif()

# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS EG RIO)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
		    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:${PROJECT_BINARY_DIR}")
# add_compile_options(-D MemoryLeak)

#----------------------------------------------------------------------------
# Copy all libraries to the install directory,
set(G4S2S_LIBS
  ${PROJECT_BINARY_DIR}/libG4S2S_lib.so
  ${PROJECT_BINARY_DIR}/libG4S2S_rdict.pcm
)

foreach(_lib ${G4S2SLIBS})
  configure_file(
    ${PROJECT_BINARY_DIR}/${_lib}
    bin/${_lib}
    COPYONLY
    )
endforeach()

add_custom_target(HYPTPC DEPENDS G4S2S)

#----------------------------------------------------------------------------
# Add the executables, and link it to the Geant4 libraries
add_executable(G4S2S main.cc ${sources} ${headers})
target_link_libraries(G4S2S ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

# Visualization-only binary (MT run manager, forced to 1 thread) so
# the RayTracer vis driver works: this Geant4 build registers it as
# G4TheMTRayTracer at compile time, which segfaults under G4S2S's
# Serial run manager. main.cc/G4S2S is untouched; production physics
# runs are unaffected. See analysis-note.md.
add_executable(G4S2SVis main_vis.cc ${sources} ${headers})
target_link_libraries(G4S2SVis ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

if(UNIX AND NOT APPLE)
  root_generate_dictionary(G4S2S_dict TParticle.h TVector3.h
		LINKDEF ${PROJECT_SOURCE_DIR}/linkdef/LinkDef.h
		MODULE G4S2S)
  add_library(G4S2S_lib SHARED G4S2S_dict)
  target_link_libraries(G4S2S G4S2S_lib)
  target_link_libraries(G4S2SVis G4S2S_lib)
  install(FILES
    ${PROJECT_BINARY_DIR}/libG4S2S_lib.so
    ${PROJECT_BINARY_DIR}/libG4S2S_rdict.pcm
    DESTINATION bin)
endif()

#----------------------------------------------------------------------------
install(TARGETS G4S2S G4S2SVis RUNTIME DESTINATION bin)
