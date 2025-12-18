/**
 * @file differential_flatness.cpp
 * @brief 基于微分平坦理论的四旋翼姿态计算
 * 
 * 题目3：差分平坦性（Differential Flatness）
 * 
 * 给定轨迹（Lemniscate双曲线）：
 * x = 10*cos(t) / (1 + sin²(t))
 * y = 10*sin(t)*cos(t) / (1 + sin²(t))
 * z = 10
 * 
 * 时间：t ∈ [0, 2π)
 * 偏航角：ψ 始终指向速度方向
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace Eigen;
using namespace std;

/**
 * @brief 微分平坦轨迹类
 */
class DifferentialFlatnessTrajectory {
private:
    double g_;  // 重力加速度 (m/s²)
    
public:
    /**
     * @brief 构造函数
     * @param g 重力加速度，默认9.81 m/s²
     */
    DifferentialFlatnessTrajectory(double g = 9.81) : g_(g) {}
    
    /**
     * @brief 计算Lemniscate曲线的位置
     * @param t 时间
     * @return 位置向量 [x, y, z]
     */
    Vector3d lemniscatePosition(double t) {
        double cos_t = cos(t);
        double sin_t = sin(t);
        double denominator = 1.0 + sin_t * sin_t;
        
        Vector3d position;
        position(0) = 10.0 * cos_t / denominator;              // x
        position(1) = 10.0 * sin_t * cos_t / denominator;      // y
        position(2) = 10.0;                                     // z
        
        return position;
    }
    
    /**
     * @brief 计算速度（一阶导数）
     * @param t 时间
     * @return 速度向量 [vx, vy, vz]
     */
    Vector3d lemniscateVelocity(double t) {
        double cos_t = cos(t);
        double sin_t = sin(t);
        double denominator = 1.0 + sin_t * sin_t;
        double derivative_denominator = 2.0 * sin_t * cos_t;
        
        Vector3d velocity;
        
        // dx/dt 使用商法则
        double numerator_x = -10.0 * sin_t;
        velocity(0) = (numerator_x * denominator - 10.0 * cos_t * derivative_denominator) 
                      / (denominator * denominator);
        
        // dy/dt
        double numerator_y = 10.0 * (cos_t * cos_t - sin_t * sin_t);
        velocity(1) = (numerator_y * denominator - 10.0 * sin_t * cos_t * derivative_denominator)
                      / (denominator * denominator);
        
        // dz/dt
        velocity(2) = 0.0;
        
        return velocity;
    }
    
    /**
     * @brief 计算加速度（二阶导数）
     * 使用数值微分方法
     * @param t 时间
     * @return 加速度向量 [ax, ay, az]
     */
    Vector3d lemniscateAcceleration(double t) {
        double dt = 1e-6;
        Vector3d v_plus = lemniscateVelocity(t + dt);
        Vector3d v_minus = lemniscateVelocity(t - dt);
        
        Vector3d acceleration = (v_plus - v_minus) / (2.0 * dt);
        
        return acceleration;
    }
    
    /**
     * @brief 使用微分平坦理论计算姿态四元数
     * 
     * 步骤：
     * 1. 计算速度 v(t) = dp/dt
     * 2. 计算加速度 a(t) = d²p/dt²
     * 3. 计算偏航角 ψ = atan2(vy, vx)
     * 4. 计算机体z轴方向：z_B = (a + g*e_z) / ||a + g*e_z||
     * 5. 构造旋转矩阵 R_WB = [x_B | y_B | z_B]
     * 6. 转换为四元数
     * 
     * @param t 时间
     * @return 四元数 Quaterniond (qw, qx, qy, qz) - 注意Eigen的顺序
     */
    Quaterniond computeAttitudeFromFlatness(double t) {
        // 1. 计算速度
        Vector3d velocity = lemniscateVelocity(t);
        
        // 2. 计算加速度
        Vector3d acceleration = lemniscateAcceleration(t);
        
        // 3. 计算偏航角（速度方向）
        double yaw = atan2(velocity(1), velocity(0));
        
        // 4. 计算推力方向（机体z轴）
        // 加上重力补偿
        Vector3d thrust_vector = acceleration + Vector3d(0, 0, g_);
        
        // 归一化得到机体z轴方向
        Vector3d z_B = thrust_vector.normalized();
        
        // 5. 构造旋转矩阵
        // 期望的x轴方向（偏航方向）
        Vector3d x_C(cos(yaw), sin(yaw), 0);
        
        // 机体y轴：z_B × x_C（右手定则）
        Vector3d y_B = z_B.cross(x_C);
        y_B.normalize();
        
        // 机体x轴：y_B × z_B
        Vector3d x_B = y_B.cross(z_B);
        
        // 旋转矩阵 R_WB = [x_B | y_B | z_B]
        Matrix3d R_WB;
        R_WB.col(0) = x_B;
        R_WB.col(1) = y_B;
        R_WB.col(2) = z_B;
        
        // 6. 转换为四元数
        Quaterniond quat(R_WB);
        
        // 归一化
        quat.normalize();
        
        return quat;
    }
    
    /**
     * @brief 生成完整轨迹的姿态序列
     * @param t_start 起始时间
     * @param t_end 结束时间
     * @param dt 时间步长
     * @param output_file 输出CSV文件名
     */
    void generateTrajectory(double t_start, double t_end, double dt, 
                           const string& output_file) {
        
        cout << "生成轨迹姿态序列..." << endl;
        cout << "时间范围: [" << fixed << setprecision(2) 
             << t_start << ", " << t_end << ")" << endl;
        cout << "时间步长: " << setprecision(3) << dt << " s" << endl;
        
        // 计算采样点数
        int num_samples = static_cast<int>((t_end - t_start) / dt);
        cout << "采样点数: " << num_samples << endl;
        cout << endl;
        
        // 打开输出文件
        ofstream file(output_file);
        if (!file.is_open()) {
            cerr << "错误：无法打开文件 " << output_file << endl;
            return;
        }
        
        // 写入CSV头
        file << "t,x,y,z,w" << endl;
        
        // 存储结果用于验证
        vector<Quaterniond> quaternions;
        vector<double> times;
        
        Quaterniond q_prev(1, 0, 0, 0);  // 初始化为单位四元数
        
        // 生成轨迹
        for (int i = 0; i < num_samples; ++i) {
            double t = t_start + i * dt;
            
            // 计算姿态四元数
            Quaterniond q = computeAttitudeFromFlatness(t);
            
            // 确保连续性：如果与前一个四元数点积为负，取反
            if (i > 0 && q.coeffs().dot(q_prev.coeffs()) < 0) {
                q.coeffs() = -q.coeffs();
            }
            
            // 确保 qw >= 0
            if (q.w() < 0) {
                q.coeffs() = -q.coeffs();
            }
            
            // 保存当前四元数
            q_prev = q;
            quaternions.push_back(q);
            times.push_back(t);
            
            // 写入文件
            // 格式：时间2位小数，四元数7位小数
            // 注意：Eigen的Quaterniond存储顺序是 (x, y, z, w)
            // 但输出要求是 t, x, y, z, w
            file << fixed << setprecision(2) << t << ","
                 << setprecision(7) << q.x() << ","
                 << setprecision(7) << q.y() << ","
                 << setprecision(7) << q.z() << ","
                 << setprecision(7) << q.w() << endl;
            
            // 显示进度
            if ((i + 1) % 50 == 0) {
                cout << "进度: " << (i + 1) << "/" << num_samples << endl;
            }
        }
        
        file.close();
        cout << "生成完成" << endl;
        cout << endl;
        
        // 验证四元数性质
        verifyQuaternionProperties(quaternions);
        
        cout << "计算结果已保存到: " << output_file << endl;
        cout << endl;
    }
    
private:
    /**
     * @brief 验证四元数性质
     */
    void verifyQuaternionProperties(const vector<Quaterniond>& quaternions) {
        cout << "============================================================" << endl;
        cout << "四元数性质验证" << endl;
        cout << "============================================================" << endl;
        
        // 1. 归一化检查
        double max_norm_error = 0.0;
        for (const auto& q : quaternions) {
            double norm = q.norm();
            double error = fabs(norm - 1.0);
            if (error > max_norm_error) {
                max_norm_error = error;
            }
        }
        
        cout << "归一化检查：" << endl;
        cout << "   ||q|| = 1 的最大偏差 = " << scientific << setprecision(2) 
             << max_norm_error << endl;
        if (max_norm_error < 1e-10) {
            cout << "通过（偏差 < 1e-10）" << endl;
        } else {
            cout << "注意：存在较大偏差" << endl;
        }
        
        // 2. qw >= 0 检查
        double min_qw = quaternions[0].w();
        for (const auto& q : quaternions) {
            if (q.w() < min_qw) {
                min_qw = q.w();
            }
        }
        
        cout << endl << "qw ≥ 0 检查：" << endl;
        cout << "最小 qw 值 = " << fixed << setprecision(7) << min_qw << endl;
        if (min_qw >= -1e-10) {
            cout << "通过（所有 qw ≥ 0）" << endl;
        } else {
            cout << "失败：存在负值" << endl;
        }
        
        // 3. 连续性检查
        double min_dot = 1.0;
        int discontinuities = 0;
        for (size_t i = 1; i < quaternions.size(); ++i) {
            double dot = quaternions[i].coeffs().dot(quaternions[i-1].coeffs());
            if (dot < min_dot) {
                min_dot = dot;
            }
            if (dot < 0) {
                discontinuities++;
            }
        }
        
        cout << endl << "连续性检查：" << endl;
        cout << "相邻四元数点积最小值 = " << setprecision(7) << min_dot << endl;
        cout << "点积为负的点数 = " << discontinuities << endl;
        if (discontinuities == 0) {
            cout << "通过（无跳变点）" << endl;
        } else {
            cout << "存在" << discontinuities << "个可能的跳变点" << endl;
        }
        
        // 4. 统计信息
        cout << endl << "四元数统计信息：" << endl;
        
        double min_qx = quaternions[0].x(), max_qx = quaternions[0].x();
        double min_qy = quaternions[0].y(), max_qy = quaternions[0].y();
        double min_qz = quaternions[0].z(), max_qz = quaternions[0].z();
        double max_qw = quaternions[0].w();
        
        double sum_qx = 0, sum_qy = 0, sum_qz = 0, sum_qw = 0;
        
        for (const auto& q : quaternions) {
            min_qx = min(min_qx, q.x()); max_qx = max(max_qx, q.x());
            min_qy = min(min_qy, q.y()); max_qy = max(max_qy, q.y());
            min_qz = min(min_qz, q.z()); max_qz = max(max_qz, q.z());
            max_qw = max(max_qw, q.w());
            
            sum_qx += q.x(); sum_qy += q.y(); sum_qz += q.z(); sum_qw += q.w();
        }
        
        int n = quaternions.size();
        cout << "   qx: 最小值=" << setprecision(6) << min_qx 
             << ", 最大值=" << max_qx << ", 均值=" << sum_qx/n << endl;
        cout << "   qy: 最小值=" << min_qy 
             << ", 最大值=" << max_qy << ", 均值=" << sum_qy/n << endl;
        cout << "   qz: 最小值=" << min_qz 
             << ", 最大值=" << max_qz << ", 均值=" << sum_qz/n << endl;
        cout << "   qw: 最小值=" << min_qw 
             << ", 最大值=" << max_qw << ", 均值=" << sum_qw/n << endl;
        
        // 5. 周期性检查
        cout << endl << "周期性检查：" << endl;
        Vector4d q_start = quaternions.front().coeffs();
        Vector4d q_end = quaternions.back().coeffs();
        
        double dist1 = (q_end - q_start).norm();
        double dist2 = (q_end + q_start).norm();
        double dist = min(dist1, dist2);
        
        cout << "   起点与终点四元数距离 = " << setprecision(6) << dist << endl;
        if (dist < 0.1) {
            cout << " 轨迹具有良好的周期性" << endl;
        } else {
            cout << " 起点终点差异较大" << endl;
        }
        
        cout << "============================================================" << endl;
        cout << endl;
    }
};

/**
 * @brief 主函数
 */
int main(int argc, char** argv) {
    // 创建差分平坦轨迹对象
    DifferentialFlatnessTrajectory df_traj(9.81);
    
    // 生成轨迹姿态序列
    // 时间范围: [0, 2π)
    // 时间步长: 0.02 s
    double t_start = 0.0;
    double t_end = 2.0 * M_PI;
    double dt = 0.02;
    
    string output_file = "df_quaternion.csv";
    
    // 如果提供了命令行参数，使用指定的输出路径
    if (argc > 1) {
        output_file = argv[1];
    }
    
    df_traj.generateTrajectory(t_start, t_end, dt, output_file);
    
    return 0;
}
