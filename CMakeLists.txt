cmake_minimum_required(VERSION 3.14)
project(Coalescence VERSION 1.0 DESCRIPTION "Chunzheng's Coalescence Model for Heavy-ion Collisions")

# 使用 C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 设置优化选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")

# 寻找ROOT
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

# 添加源文件和头文件的路径
set(SOURCES
    src/main.cxx
    src/EventsReader.cxx
    src/EventsWriter.cxx
    src/Coalescence.cxx
    src/CalculateObvs.cxx
    src/Par.cxx
)

# 设置可执行文件的输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

# 生成可执行文件并链接ROOT库
add_executable(Coalescence ${SOURCES})

# 添加 include 目录到包含路径
target_include_directories(Coalescence PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

# 链接 ROOT 库
target_link_libraries(Coalescence PRIVATE ${ROOT_LIBRARIES})