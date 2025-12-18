#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
题目3：差分平坦性（Differential Flatness）
使用微分平坦理论从给定轨迹推导四旋翼姿态

给定轨迹（Lemniscate双曲线）：
x = 10*cos(t) / (1 + sin²(t))
y = 10*sin(t)*cos(t) / (1 + sin²(t))
z = 10

时间：t ∈ [0, 2π)
偏航角：ψ 始终指向速度方向
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R

# 设置matplotlib支持中文
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

class DifferentialFlatnessTrajectory:
    """
    基于微分平坦理论的四旋翼轨迹姿态计算
    """
    
    def __init__(self, g=9.81):
        """
        初始化
        
        参数：
        g: 重力加速度 (m/s²)
        """
        self.g = g
    
    def lemniscate_position(self, t):
        """
        计算Lemniscate曲线的位置
        
        x = 10*cos(t) / (1 + sin²(t))
        y = 10*sin(t)*cos(t) / (1 + sin²(t))
        z = 10
        
        参数：
        t: 时间 (标量或数组)
        
        返回：
        position: (3,) 或 (N, 3) 数组
        """
        cos_t = np.cos(t)
        sin_t = np.sin(t)
        denominator = 1 + sin_t**2
        
        x = 10 * cos_t / denominator
        y = 10 * sin_t * cos_t / denominator
        z = 10 * np.ones_like(t)
        
        if np.isscalar(t):
            return np.array([x, y, z])
        else:
            return np.column_stack([x, y, z])
    
    def lemniscate_velocity(self, t):
        """
        计算速度（一阶导数）
        
        使用链式法则求导：
        dx/dt = d/dt [10*cos(t) / (1 + sin²(t))]
        dy/dt = d/dt [10*sin(t)*cos(t) / (1 + sin²(t))]
        dz/dt = 0
        """
        cos_t = np.cos(t)
        sin_t = np.sin(t)
        denominator = 1 + sin_t**2
        
        # x的导数
        # 分子: d/dt[10*cos(t)] = -10*sin(t)
        # 分母: d/dt[1 + sin²(t)] = 2*sin(t)*cos(t)
        # 商法则: (u/v)' = (u'v - uv') / v²
        numerator_x = -10 * sin_t
        denominator_x = denominator
        derivative_denominator = 2 * sin_t * cos_t
        
        dx_dt = (numerator_x * denominator_x - 10 * cos_t * derivative_denominator) / (denominator**2)
        
        # y的导数
        # y = 10*sin(t)*cos(t) / (1 + sin²(t))
        # 分子导数: d/dt[10*sin(t)*cos(t)] = 10*(cos²(t) - sin²(t))
        numerator_y = 10 * (cos_t**2 - sin_t**2)
        
        dy_dt = (numerator_y * denominator - 10 * sin_t * cos_t * derivative_denominator) / (denominator**2)
        
        # z的导数
        dz_dt = np.zeros_like(t)
        
        if np.isscalar(t):
            return np.array([dx_dt, dy_dt, dz_dt])
        else:
            return np.column_stack([dx_dt, dy_dt, dz_dt])
    
    def lemniscate_acceleration(self, t):
        """
        计算加速度（二阶导数）
        
        对速度再求导得到加速度
        """
        # 使用数值微分
        dt = 1e-6
        v_plus = self.lemniscate_velocity(t + dt)
        v_minus = self.lemniscate_velocity(t - dt)
        
        acceleration = (v_plus - v_minus) / (2 * dt)
        
        return acceleration
    
    def compute_attitude_from_flatness(self, t):
        """
        使用微分平坦理论计算姿态
        
        步骤：
        1. 计算位置 p(t)
        2. 计算速度 v(t) = dp/dt
        3. 计算加速度 a(t) = d²p/dt²
        4. 计算偏航角 ψ = atan2(vy, vx)
        5. 计算机体z轴方向：z_B = (a + g*e_z) / ||a + g*e_z||
        6. 构造旋转矩阵 R_WB
        7. 转换为四元数
        
        参数：
        t: 时间
        
        返回：
        quaternion: [qx, qy, qz, qw]
        """
        # 1. 计算位置
        position = self.lemniscate_position(t)
        
        # 2. 计算速度
        velocity = self.lemniscate_velocity(t)
        
        # 3. 计算加速度
        acceleration = self.lemniscate_acceleration(t)
        
        # 4. 计算偏航角（速度方向）
        yaw = np.arctan2(velocity[1], velocity[0])
        
        # 5. 计算推力方向（机体z轴）
        # 加上重力补偿
        thrust_vector = acceleration + np.array([0, 0, self.g])
        
        # 归一化得到机体z轴方向
        z_B = thrust_vector / np.linalg.norm(thrust_vector)
        
        # 6. 构造旋转矩阵
        # 期望的x轴方向（偏航方向）
        x_C = np.array([np.cos(yaw), np.sin(yaw), 0])
        
        # 机体y轴：z_B × x_C（右手定则）
        y_B = np.cross(z_B, x_C)
        y_B = y_B / np.linalg.norm(y_B)
        
        # 机体x轴：y_B × z_B
        x_B = np.cross(y_B, z_B)
        
        # 旋转矩阵 R_WB = [x_B | y_B | z_B]
        R_WB = np.column_stack([x_B, y_B, z_B])
        
        # 7. 转换为四元数
        rot = R.from_matrix(R_WB)
        quat = rot.as_quat()  # [qx, qy, qz, qw]
        
        return quat
    
    def generate_trajectory(self, t_start=0.0, t_end=2*np.pi, dt=0.02):
        """
        生成完整轨迹的姿态序列
        
        参数：
        t_start: 起始时间
        t_end: 结束时间
        dt: 时间步长
        
        返回：
        DataFrame with columns: t, qx, qy, qz, qw
        """
        # 时间采样
        time_samples = np.arange(t_start, t_end, dt)
        
        print(f"生成轨迹姿态序列...")
        print(f"时间范围: [{t_start:.2f}, {t_end:.2f})")
        print(f"时间步长: {dt:.3f} s")
        print(f"采样点数: {len(time_samples)}")
        
        results = []
        q_prev = None
        
        for i, t in enumerate(time_samples):
            # 计算姿态四元数
            q = self.compute_attitude_from_flatness(t)
            
            # 归一化
            q = q / np.linalg.norm(q)
            
            # 确保连续性
            if q_prev is not None:
                if np.dot(q, q_prev) < 0:
                    q = -q
            
            # 确保 qw >= 0
            if q[3] < 0:
                q = -q
            
            q_prev = q.copy()
            
            # 保存结果
            results.append({
                't': t,
                'qx': q[0],
                'qy': q[1],
                'qz': q[2],
                'qw': q[3]
            })
            
            # 显示进度
            if (i + 1) % 50 == 0:
                print(f"  进度: {i + 1}/{len(time_samples)}")
        
        print("生成完成！\n")
        
        return pd.DataFrame(results)


def plot_quaternion_evolution(df, save_path='problem3_quaternion_plot.png'):
    """
    绘制四元数随时间的变化曲线
    """
    fig, axes = plt.subplots(4, 1, figsize=(14, 12))
    fig.suptitle('基于微分平坦理论的四旋翼姿态四元数变化', 
                 fontsize=18, fontweight='bold', y=0.995)
    
    components = ['qx', 'qy', 'qz', 'qw']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    labels = ['$q_x$', '$q_y$', '$q_z$', '$q_w$']
    
    for i, (comp, color, label) in enumerate(zip(components, colors, labels)):
        axes[i].plot(df['t'], df[comp], color=color, linewidth=2.5, label=label)
        axes[i].set_ylabel(label, fontsize=14, fontweight='bold')
        axes[i].grid(True, alpha=0.4, linestyle='--', linewidth=0.8)
        axes[i].legend(loc='upper right', fontsize=12, framealpha=0.9)
        axes[i].tick_params(labelsize=11)
        
        # 设置y轴范围
        y_min, y_max = df[comp].min(), df[comp].max()
        y_margin = max((y_max - y_min) * 0.1, 0.05)
        axes[i].set_ylim(y_min - y_margin, y_max + y_margin)
        
        # 只在最后一个子图显示x轴标签
        if i == 3:
            axes[i].set_xlabel('时间 (s)', fontsize=14, fontweight='bold')
        else:
            axes[i].set_xticklabels([])
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"✓ 四元数变化曲线图已保存到: {save_path}")


def plot_trajectory_3d(df_traj, save_path='problem3_trajectory_3d.png'):
    """
    绘制三维轨迹
    """
    from mpl_toolkits.mplot3d import Axes3D
    
    # 创建微分平坦对象计算轨迹
    df_obj = DifferentialFlatnessTrajectory()
    positions = df_obj.lemniscate_position(df_traj['t'].values)
    
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # 绘制轨迹
    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], 
            'b-', linewidth=2.5, label='Lemniscate轨迹')
    
    # 标记起点和终点
    ax.scatter(positions[0, 0], positions[0, 1], positions[0, 2], 
              c='g', s=100, marker='o', label='起点')
    ax.scatter(positions[-1, 0], positions[-1, 1], positions[-1, 2], 
              c='r', s=100, marker='s', label='终点')
    
    ax.set_xlabel('X (m)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Y (m)', fontsize=12, fontweight='bold')
    ax.set_zlabel('Z (m)', fontsize=12, fontweight='bold')
    ax.set_title('Lemniscate双曲线轨迹（俯视平面内）', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # 设置相等的坐标轴比例
    max_range = np.array([positions[:, 0].max() - positions[:, 0].min(),
                          positions[:, 1].max() - positions[:, 1].min(),
                          1.0]).max() / 2.0
    mid_x = (positions[:, 0].max() + positions[:, 0].min()) * 0.5
    mid_y = (positions[:, 1].max() + positions[:, 1].min()) * 0.5
    mid_z = 10.0
    
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - 1, mid_z + 1)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"✓ 三维轨迹图已保存到: {save_path}")


def verify_quaternion_properties(df):
    """
    验证四元数的数学性质
    """
    print("=" * 60)
    print("四元数性质验证")
    print("=" * 60)
    
    # 1. 归一化检查
    norms = np.sqrt(df['qx']**2 + df['qy']**2 + df['qz']**2 + df['qw']**2)
    max_norm_error = np.max(np.abs(norms - 1.0))
    print(f"1. 归一化检查：")
    print(f"   ||q|| = 1 的最大偏差 = {max_norm_error:.2e}")
    if max_norm_error < 1e-10:
        print("   ✓ 通过（偏差 < 1e-10）")
    else:
        print("   ⚠ 注意：存在较大偏差")
    
    # 2. qw >= 0 检查
    min_qw = df['qw'].min()
    print(f"\n2. qw ≥ 0 检查：")
    print(f"   最小 qw 值 = {min_qw:.7f}")
    if min_qw >= -1e-10:
        print("   ✓ 通过（所有 qw ≥ 0）")
    else:
        print(f"   ✗ 失败：存在 {np.sum(df['qw'] < 0)} 个负值")
    
    # 3. 连续性检查
    q_array = df[['qx', 'qy', 'qz', 'qw']].values
    dot_products = np.sum(q_array[:-1] * q_array[1:], axis=1)
    min_dot = np.min(dot_products)
    discontinuities = np.sum(dot_products < 0)
    
    print(f"\n3. 连续性检查：")
    print(f"   相邻四元数点积最小值 = {min_dot:.7f}")
    print(f"   点积为负的点数 = {discontinuities}")
    if discontinuities == 0:
        print("   ✓ 通过（无跳变点）")
    else:
        print(f"   ⚠ 存在 {discontinuities} 个可能的跳变点")
    
    # 4. 统计信息
    print(f"\n4. 四元数统计信息：")
    for comp in ['qx', 'qy', 'qz', 'qw']:
        print(f"   {comp}: 最小值={df[comp].min():.6f}, "
              f"最大值={df[comp].max():.6f}, "
              f"均值={df[comp].mean():.6f}")
    
    # 5. 周期性检查
    print(f"\n5. 周期性检查：")
    q_start = df[['qx', 'qy', 'qz', 'qw']].iloc[0].values
    q_end = df[['qx', 'qy', 'qz', 'qw']].iloc[-1].values
    
    # 四元数距离（考虑q和-q等价）
    dist1 = np.linalg.norm(q_end - q_start)
    dist2 = np.linalg.norm(q_end + q_start)
    dist = min(dist1, dist2)
    
    print(f"   起点与终点四元数距离 = {dist:.6f}")
    if dist < 0.1:
        print("   ✓ 轨迹具有良好的周期性")
    else:
        print(f"   ⚠ 起点终点差异较大，可能不是完美周期")
    
    print("=" * 60)


def main():
    """
    主函数
    """
    print("=" * 60)
    print("题目3：差分平坦性（Differential Flatness）")
    print("=" * 60)
    print()
    
    # 创建差分平坦对象
    df_obj = DifferentialFlatnessTrajectory(g=9.81)
    
    # 生成姿态序列
    df_result = df_obj.generate_trajectory(
        t_start=0.0,
        t_end=2*np.pi,
        dt=0.02
    )
    
    print(f"生成结果（前10行）：")
    print(df_result.head(10).to_string(index=False))
    print(f"\n总共 {len(df_result)} 个时间点")
    
    # 验证四元数性质
    print()
    verify_quaternion_properties(df_result)
    
    # 保存结果到CSV
    output_csv = 'C:/Users/luhown/Documents/GitHub/MRPC-2025-homework/documents/df_quaternion.csv'
    
    # 格式化输出：时间2位小数，四元数7位小数
    df_output = df_result.copy()
    df_output['t'] = df_output['t'].round(2)
    df_output[['qx', 'qy', 'qz', 'qw']] = df_output[['qx', 'qy', 'qz', 'qw']].round(7)
    
    df_output.to_csv(output_csv, index=False, float_format='%.7f')
    print(f"\n✓ 计算结果已保存到: {output_csv}")
    
    # 绘制四元数变化曲线
    print("\n正在绘制四元数变化曲线...")
    plot_path_quat = 'C:/Users/luhown/Documents/GitHub/MRPC-2025-homework/documents/problem3_quaternion_plot.png'
    plot_quaternion_evolution(df_result, save_path=plot_path_quat)
    
    # 绘制三维轨迹
    print("正在绘制三维轨迹...")
    plot_path_traj = 'C:/Users/luhown/Documents/GitHub/MRPC-2025-homework/documents/problem3_trajectory_3d.png'
    plot_trajectory_3d(df_result, save_path=plot_path_traj)
    
    print("\n" + "=" * 60)
    print("题目3 完成！")
    print("=" * 60)
    
    return df_result


if __name__ == "__main__":
    results = main()
