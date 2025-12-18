cmake_minimum_required(VERSION 3.0.2)
project(quadrotor_df)

## C++14标准
add_compile_options(-std=c++14)

## 查找catkin和所需组件
find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  geometry_msgs
)

## 查找Eigen3
find_package(Eigen3 REQUIRED)

###################################
## catkin specific configuration ##
###################################
catkin_package(
  CATKIN_DEPENDS roscpp std_msgs geometry_msgs
  DEPENDS EIGEN3
)

###########
## Build ##
###########

## 指定包含目录
include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
)

## 声明C++可执行文件
add_executable(differential_flatness_node 
  src/differential_flatness.cpp
)

## 链接库
target_link_libraries(differential_flatness_node
  ${catkin_LIBRARIES}
)

## 确保消息头文件在编译前生成
add_dependencies(differential_flatness_node 
  ${${PROJECT_NAME}_EXPORTED_TARGETS} 
  ${catkin_EXPORTED_TARGETS}
)

#############
## Install ##
#############

## 安装可执行文件
install(TARGETS differential_flatness_node
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
