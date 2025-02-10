cmake_minimum_required(VERSION 2.8)

set(SOURCE tzap.c integrate.c calc.c test_v0.c test_v1.c backtrace.c main.c)
set(HEADER tzap.h integrate.h calc.h test_v0.h test_v1.h dbg_info.h backtrace.h)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_DEBUG")

add_executable(testtzap ${SOURCE} ${HEADER})