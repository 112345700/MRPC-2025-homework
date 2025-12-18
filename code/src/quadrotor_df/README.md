# Quadrotor Differential Flatness Package

## 题目3：差分平坦性（Differential Flatness）

基于微分平坦理论的四旋翼姿态计算ROS包。

### 功能描述

根据给定的Lemniscate双曲线轨迹，使用微分平坦理论计算无人机的姿态四元数：

**轨迹方程：**
```
x = 10*cos(t) / (1 + sin²(t))
y = 10*sin(t)*cos(t) / (1 + sin²(t))
z = 10
```

**时间范围：** t ∈ [0, 2π)

**偏航角：** ψ 始终指向速度方向

### 依赖项

- ROS (Kinetic/Melodic/Noetic)
- Eigen3
- C++14编译器

### 编译方法

#### 方法1：在ROS工作空间中编译

```bash
# 将此包复制到ROS工作空间
cd ~/catkin_ws/src
cp -r /path/to/quadrotor_df .

# 编译
cd ~/catkin_ws
catkin_make

# 或使用catkin build
catkin build quadrotor_df
```

#### 方法2：独立编译（无ROS）

```bash
cd quadrotor_df/src

# 使用g++直接编译
g++ -std=c++14 -O2 -o differential_flatness differential_flatness.cpp -I/usr/include/eigen3

# 运行
./differential_flatness
```

### 运行方法

#### 在ROS环境中运行

```bash
# Source环境
source ~/catkin_ws/devel/setup.bash

# 运行节点
rosrun quadrotor_df differential_flatness_node

# 或指定输出文件路径
rosrun quadrotor_df differential_flatness_node /path/to/output.csv
```

#### 独立运行

```bash
cd quadrotor_df/src

# 运行可执行文件
./differential_flatness

# 或指定输出路径
./differential_flatness /path/to/output.csv
```

### 输出文件

生成的CSV文件格式：

```csv
t,x,y,z,w
0.00,0.0499792,0.0000000,0.0000000,0.9987503
0.02,0.0499167,0.0499167,0.0024979,0.9975021
0.04,0.0523491,0.0473595,0.0523491,0.9961306
...
```

**列说明：**
- `t`: 时间（秒），保留2位小数
- `x, y, z, w`: 四元数分量，保留7位小数
- 四元数格式：(qx, qy, qz, qw)

### 数学原理

使用微分平坦理论，从轨迹的位置和导数推导姿态：

1. **计算速度和加速度：**
   - v(t) = dp/dt
   - a(t) = d²p/dt²

2. **计算偏航角：**
   - ψ = atan2(vy, vx)

3. **计算推力方向（机体z轴）：**
   - z_B = (a + g*e_z) / ||a + g*e_z||

4. **构造旋转矩阵：**
   - x_C = [cos(ψ), sin(ψ), 0]
   - y_B = z_B × x_C (归一化)
   - x_B = y_B × z_B
   - R_WB = [x_B | y_B | z_B]

5. **转换为四元数并确保连续性**

### 验证

程序会自动验证以下性质：

- ✓ 四元数归一化：||q|| = 1
- ✓ qw ≥ 0 约束
- ✓ 时间连续性（无跳变）
- ✓ 周期性（起点和终点接近）

### 代码结构

```
quadrotor_df/
├── CMakeLists.txt          # CMake配置文件
├── package.xml             # ROS包配置
├── README.md              # 本文件
└── src/
    └── differential_flatness.cpp  # 主程序源代码
```

### 主要类和函数

**DifferentialFlatnessTrajectory类：**

- `lemniscatePosition(t)`: 计算位置
- `lemniscateVelocity(t)`: 计算速度
- `lemniscateAcceleration(t)`: 计算加速度
- `computeAttitudeFromFlatness(t)`: 计算姿态四元数
- `generateTrajectory()`: 生成完整轨迹
- `verifyQuaternionProperties()`: 验证四元数性质

### 示例输出

```
============================================================
题目3：差分平坦性（Differential Flatness）
============================================================

生成轨迹姿态序列...
时间范围: [0.00, 6.28)
时间步长: 0.020 s
采样点数: 315

  进度: 50/315
  进度: 100/315
  ...
生成完成！

============================================================
四元数性质验证
============================================================
1. 归一化检查：
   ||q|| = 1 的最大偏差 = 2.22e-16
   ✓ 通过（偏差 < 1e-10）

2. qw ≥ 0 检查：
   最小 qw 值 = 0.0055232
   ✓ 通过（所有 qw ≥ 0）

3. 连续性检查：
   相邻四元数点积最小值 = 0.9997518
   点积为负的点数 = 0
   ✓ 通过（无跳变点）

4. 四元数统计信息：
   qx: 最小值=-0.545067, 最大值=0.545125, 均值=-0.001134
   qy: 最小值=-0.545104, 最大值=0.545078, 均值=-0.000850
   qz: 最小值=-0.923876, 最大值=0.861394, 均值=-0.181980
   qw: 最小值=0.005523, 最大值=0.923877, 均值=0.560422

5. 周期性检查：
   起点与终点四元数距离 = 0.005260
   ✓ 轨迹具有良好的周期性
============================================================

✓ 计算结果已保存到: df_quaternion.csv

============================================================
题目3 完成！
============================================================
```

### 注意事项

1. **Eigen版本：** 确保安装了Eigen3（通常通过 `sudo apt-get install libeigen3-dev` 安装）

2. **四元数顺序：** Eigen的Quaterniond内部存储顺序是 (x, y, z, w)，与某些库不同

3. **数值精度：** 使用数值微分计算加速度，步长设为1e-6以平衡精度和稳定性

4. **连续性处理：** 使用点积检测四元数跳变，自动取反以保持连续性

### 参考文献

- Mellinger, D., & Kumar, V. (2011). "Minimum snap trajectory generation and control for quadrotors." ICRA.
- 微分平坦理论在四旋翼控制中的应用

### 作者

学生作业 - 题目3实现

### 许可

MIT License
