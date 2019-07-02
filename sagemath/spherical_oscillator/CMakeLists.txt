cmake_minimum_required(VERSION 2.8)

set(SOURCE tzap.c testtzap.c)
set(HEADER tzap.h dbg_info.h)

add_executable(testtzap ${SOURCE} ${HEADER})