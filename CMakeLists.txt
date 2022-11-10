cmake_minimum_required(VERSION 3.1.0)
project(lpmatrixgrounder VERSION 0.1.0)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Og")
set(CMAKE_CXX_FLAGS_RELEASE "-march=native -O3 -pipe")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CTest)
enable_testing()

add_library(${PROJECT_NAME} SHARED
  src/lpmatrixgrounder.cpp
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)


add_subdirectory(externals/clingo)

target_include_directories(${PROJECT_NAME}
  PUBLIC ${PROJECT_BINARY_DIR}
  PUBLIC externals/clingo/libclingo/
)

target_link_directories(${PROJECT_NAME}
# looking for built libclingo.a binary
  PRIVATE bin/
)

target_link_libraries(${PROJECT_NAME}
PUBLIC
  libclingo
)


# a testing prgram for the lib
add_executable(testlpmatrixgrounder
    src/main.cpp 
)
target_link_libraries(testlpmatrixgrounder
PUBLIC
  libclingo
  ${PROJECT_NAME}
)