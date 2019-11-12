cmake_minimum_required(VERSION 2.8)

set(SOURCE tzap.c testtzap.c backtrace.c)
set(HEADER tzap.h dbg_info.h backtrace.h)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_DEBUG")

add_executable(testtzap ${SOURCE} ${HEADER})