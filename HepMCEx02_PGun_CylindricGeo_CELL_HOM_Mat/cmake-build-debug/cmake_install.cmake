# Install script for directory: /Users/antonc/Documents/GitHub/Master/InOne/HepMCEx02_PGun_CylindricGeo_CELL_HOM_Mat

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/antonc/Documents/GitHub/Master/InOne/HepMCEx02_PGun_CylindricGeo_CELL_HOM_Mat/cmake-build-debug/HepMCEx02")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Applications/root_v6.18.04/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/antonc/Documents/geant4.10.5-install/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/antonc/Programs/HEP_SOFTWARES/HEPMC/HepMC-install/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/antonc/Programs/HEP_SOFTWARES/PYTHIA/pythia8243/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/antonc/Documents/GitHub/Master/InOne/HepMCEx02_PGun_CylindricGeo_CELL_HOM_Mat/cmake-build-debug"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/antonc/Documents/GitHub/Master/InOne/HepMCEx02_PGun_CylindricGeo_CELL_HOM_Mat/cmake-build-debug/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
