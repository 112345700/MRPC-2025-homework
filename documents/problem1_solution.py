#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
题目1：坐标系转换
计算执行器在世界坐标系下的姿态（四元数表示）

已知参数：
- ω = 0.5 rad/s  (角频率)
- α = π/12       (圆锥半角)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R

# 设置matplotlib支持中文
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def quaternion_multiply(q1, q2):
    """
    四元数乘法: q1 * q2
    四元数格式: [qx, qy, qz, qw]
    
    使用Hamilton乘积规则
    """
    x1, y1, z1, w1 = q1
    x2, y2, z2, w2 = q2
    
    w = w1*w2 - x1*x2 - y1*y2 - z1*z2
    x = w1*x2 + x1*w2 + y1*z2 - z1*y2
    y = w1*y2 - x1*z2 + y1*w2 + z1*x2
    z = w1*z2 + x1*y2 - y1*x2 + z1*w2
    
    return np.array([x, y, z, w])

def normalize_quaternion(q):
    """
    归一化四元数，确保 ||q|| = 1
    """
    norm = np.linalg.norm(q)
    if norm < 1e-10:
        print("警告：四元数模长接近零")
        return np.array([0, 0, 0, 1])
    return q / norm

def ensure_continuous_quaternion(q, q_prev):
    """
    确保四元数的连续性
    
    四元数 q 和 -q 代表相同的旋转，但为了保持连续性，
    如果 q·q_prev < 0，需要将 q 取反
    """
    if np.dot(q, q_prev) < 0:
        return -q
    return q

def rotation_matrix_to_quaternion(R_matrix):
    """
    旋转矩阵转四元数
    使用scipy的Rotation类确保稳定性
    
    返回格式: [qx, qy, qz, qw]
    """
    rot = R.from_matrix(R_matrix)
    quat = rot.as_quat()  # scipy返回 [x, y, z, w]
    return quat

def compute_end_effector_rotation(t, omega=0.5, alpha=np.pi/12):
    """
    计算末端执行器相对于机体系的旋转矩阵
    
    根据题目给出的公式(1)：
    
    B_R_D = ⎡cos(ωt)   -sin(ωt)cos(α)   sin(ωt)sin(α) ⎤
            ⎢sin(ωt)    cos(ωt)cos(α)  -cos(ωt)sin(α) ⎥
            ⎣   0           sin(α)           cos(α)     ⎦
    
    参数：
    t: 时间 (s)
    omega: 角频率 (rad/s)，题目要求 ω = 0.5
    alpha: 圆锥半角 (rad)，题目要求 α = π/12
    
    返回：
    B_R_D: 机体系到执行器系的旋转矩阵 (3x3)
    """
    wt = omega * t
    cos_wt = np.cos(wt)
    sin_wt = np.sin(wt)
    cos_alpha = np.cos(alpha)
    sin_alpha = np.sin(alpha)
    
    # 严格按照题目公式构造旋转矩阵
    B_R_D = np.array([
        [cos_wt,  -sin_wt * cos_alpha,   sin_wt * sin_alpha],
        [sin_wt,   cos_wt * cos_alpha,  -cos_wt * sin_alpha],
        [0,        sin_alpha,             cos_alpha]
    ])
    
    return B_R_D

def compute_world_frame_orientation(tracking_file, omega=0.5, alpha=np.pi/12):
    """
    计算执行器在世界坐标系下的姿态
    
    参数：
    tracking_file: 包含无人机姿态的CSV文件路径
    omega: 角频率，默认 0.5 rad/s
    alpha: 圆锥半角，默认 π/12
    
    返回：
    results_df: 包含时间和世界系下四元数的DataFrame
    """
    print("="*60)
    print("题目1：坐标系转换")
    print("="*60)
    print(f"参数设置：ω = {omega:.1f} rad/s, α = π/12 = {np.rad2deg(alpha):.1f}°")
    print("="*60)
    
    # 读取无人机姿态数据
    df = pd.read_csv(tracking_file)
    print(f"\n读取到 {len(df)} 个时间点的姿态数据")
    print(f"时间范围：{df['t'].min():.2f}s ~ {df['t'].max():.2f}s")
    
    # 初始化结果列表
    results = []
    q_prev = None
    
    print("\n正在计算执行器在世界坐标系下的姿态...")
    
    for idx, row in df.iterrows():
        t = row['t']
        
        # 步骤1：读取无人机姿态四元数 (世界系到机体系)
        # q_WB = [qx, qy, qz, qw]
        q_WB = np.array([row['qx'], row['qy'], row['qz'], row['qw']])
        q_WB = normalize_quaternion(q_WB)
        
        # 步骤2：计算执行器相对于机体系的旋转矩阵
        B_R_D = compute_end_effector_rotation(t, omega, alpha)
        
        # 步骤3：旋转矩阵转换为四元数 (机体系到执行器系)
        q_BD = rotation_matrix_to_quaternion(B_R_D)
        
        # 步骤4：四元数复合计算世界系到执行器系的四元数
        # 旋转复合：W_R_D = W_R_B * B_R_D
        # 四元数复合：q_WD = q_WB ⊗ q_BD
        q_WD = quaternion_multiply(q_WB, q_BD)
        
        # 步骤5：归一化
        q_WD = normalize_quaternion(q_WD)
        
        # 步骤6：确保连续性
        if q_prev is not None:
            q_WD = ensure_continuous_quaternion(q_WD, q_prev)
        
        # 步骤7：确保 qw ≥ 0
        if q_WD[3] < 0:
            q_WD = -q_WD
        
        # 更新前一个四元数
        q_prev = q_WD.copy()
        
        # 保存结果
        results.append({
            't': t,
            'qx': q_WD[0],
            'qy': q_WD[1],
            'qz': q_WD[2],
            'qw': q_WD[3]
        })
        
        # 显示进度（每100个点）
        if (idx + 1) % 100 == 0:
            print(f"  进度: {idx + 1}/{len(df)}")
    
    print("计算完成！\n")
    
    return pd.DataFrame(results)

def plot_quaternion_evolution(df, save_path='problem1_quaternion_plot.png'):
    """
    绘制四元数随时间的变化曲线
    
    参数：
    df: 包含时间和四元数的DataFrame
    save_path: 图片保存路径
    """
    fig, axes = plt.subplots(4, 1, figsize=(14, 12))
    fig.suptitle('执行器在世界坐标系下的姿态四元数变化曲线', 
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
        
        # 设置合适的y轴范围
        y_min, y_max = df[comp].min(), df[comp].max()
        y_margin = (y_max - y_min) * 0.1
        axes[i].set_ylim(y_min - y_margin, y_max + y_margin)
        
        # 只在最后一个子图显示x轴标签
        if i == 3:
            axes[i].set_xlabel('时间 (s)', fontsize=14, fontweight='bold')
        else:
            axes[i].set_xticklabels([])
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"✓ 四元数变化曲线图已保存到: {save_path}")
    
    return fig

def verify_quaternion_properties(df):
    """
    验证四元数的各项性质
    """
    print("="*60)
    print("四元数性质验证")
    print("="*60)
    
    # 1. 归一化检查
    norms = np.sqrt(df['qx']**2 + df['qy']**2 + df['qz']**2 + df['qw']**2)
    max_norm_error = np.max(np.abs(norms - 1.0))
    print(f"1. 归一化检查：")
    print(f"   ||q|| = 1 的最大偏差 = {max_norm_error:.2e}")
    if max_norm_error < 1e-10:
        print("   ✓ 通过（偏差 < 1e-10）")
    else:
        print("   注意：存在较大偏差")
    
    # 2. qw ≥ 0 检查
    min_qw = df['qw'].min()
    print(f"\n2. qw ≥ 0 检查：")
    print(f"   最小 qw 值 = {min_qw:.7f}")
    if min_qw >= -1e-10:  # 允许微小数值误差
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
        print(f"   存在 {discontinuities} 个可能的跳变点")
    
    # 4. 统计信息
    print(f"\n4. 四元数统计信息：")
    for comp in ['qx', 'qy', 'qz', 'qw']:
        print(f"   {comp}: 最小值={df[comp].min():.6f}, "
              f"最大值={df[comp].max():.6f}, "
              f"均值={df[comp].mean():.6f}")
    
    print("="*60)

def main():
    """
    主函数
    """
    # 输入文件路径
    tracking_file = 'C:/Users/luhown/Documents/GitHub/MRPC-2025-homework/documents/tracking.csv'
    
    # 题目参数
    omega = 0.5  # rad/s
    alpha = np.pi / 12  # rad
    
    # 计算世界系下的执行器姿态
    results_df = compute_world_frame_orientation(tracking_file, omega, alpha)
    
    # 输出前10行查看结果
    print("\n计算结果（前10行）：")
    print(results_df.head(10).to_string(index=False))
    print(f"\n总共 {len(results_df)} 个时间点")
    
    # 验证四元数性质
    print()
    verify_quaternion_properties(results_df)
    
    # 绘制四元数变化曲线
    print("\n正在绘制四元数变化曲线...")
    plot_path = 'C:/Users/luhown/Documents/GitHub/MRPC-2025-homework/documents/problem1_quaternion_plot.png'
    plot_quaternion_evolution(results_df, save_path=plot_path)
    
    # 保存结果到CSV文件
    output_csv = 'C:/Users/luhown/Documents/GitHub/MRPC-2025-homework/documents/problem1_world_frame_quaternion.csv'
    results_df.to_csv(output_csv, index=False, float_format='%.10f')
    print(f"✓ 计算结果已保存到: {output_csv}")
    
    print("\n" + "="*60)
    print("题目1 完成！")
    print("="*60)
    
    return results_df

if __name__ == "__main__":
    results = main()
