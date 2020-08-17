# Install script for directory: /Users/sanmay/GEANT/Tutorial/HepMC/HepMCEx02_PGun_RectangularGeo_CELL

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
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/sanmay/GEANT/Tutorial/HepMC/HepMCEx02_PGun_RectangularGeo_CELL/build/HepMCEx02")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/sanmay/NEW_ROOT/root-install/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/sanmay/GEANT/geant4.10.05-install/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/sanmay/GEANT/Tutorial/HepMC/HepMCEx02_PGun_RectangularGeo_CELL/build"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/anaconda3/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/sanmay/CLHEP/CLHEP_INSTALL/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HepMCEx02")
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
file(WRITE "/Users/sanmay/GEANT/Tutorial/HepMC/HepMCEx02_PGun_RectangularGeo_CELL/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
