cmake_minimum_required(VERSION 3.5)

# list of source files
set(dgpTest2a_files dgpTest2a.cpp dgpPrt.cpp)
set(dgpTest2b_files dgpTest2b.cpp dgpPrt.cpp)
set(dgpTest2c_files dgpTest2c.cpp dgpPrt.cpp)

# define the executable
if(WIN32)
  add_executable(dgpTest2a WIN32 ${dgpTest2a_files})
  add_executable(dgpTest2b WIN32 ${dgpTest2b_files})
  add_executable(dgpTest2c WIN32 ${dgpTest2c_files})
else()
  add_executable(dgpTest2a ${dgpTest2a_files})
  add_executable(dgpTest2b ${dgpTest2b_files})
  add_executable(dgpTest2c ${dgpTest2c_files})
endif()

# in Windows + Visual Studio we need this to make it a console application
if(WIN32)
  if(MSVC)
    set_target_properties(dgpTest2a PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
    set_target_properties(dgpTest2b PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
    set_target_properties(dgpTest2c PROPERTIES LINK_FLAGS "/SUBSYSTEM:CONSOLE")
  endif(MSVC)
endif(WIN32)

# install application
set(BIN_DIR ${CMAKE_INSTALL_PREFIX}/bin)

# target_link_libraries(dgpTest2 core io util wrl)
target_link_libraries(dgpTest2a ${LIB_LIST})
target_link_libraries(dgpTest2b ${LIB_LIST})
target_link_libraries(dgpTest2c ${LIB_LIST})

install(TARGETS dgpTest2a DESTINATION ${BIN_DIR})
install(TARGETS dgpTest2b DESTINATION ${BIN_DIR})
install(TARGETS dgpTest2c DESTINATION ${BIN_DIR})

