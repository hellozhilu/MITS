cmake_minimum_required(VERSION 3.17)

# 设置项目名称
set(PROJECT_NAME BKPWC)
project(${PROJECT_NAME})

# 设置代码路径
set(CODE_PATH src)

# 设置编译器
set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

# 设置编译选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

#设置可执行文件输出路径
set(EXECUTABLE_OUTPUT_PATH "../")

#包含头文件路径
include_directories(${CODE_PATH})

# 递归查找源代码文件
# file(GLOB_RECURSE CODE_FILES ${CODE_PATH}/*.*)
file(GLOB_RECURSE CODE_FILES
    ${CODE_PATH}/*.cpp
    ${CODE_PATH}/*.h
#    ${CODE_PATH}/*.hpp
#    ${CODE_PATH}/*.c
)

# 添加可执行文件
add_executable(${PROJECT_NAME} ${CODE_FILES})