#include "Astar_searcher.h"

using namespace std;
using namespace Eigen;


static const int    kInflateZLayers = 0;
static const int    kObsCostRadius = 1;
static const double kObsCostWeight = 0.2;
static const bool   kEnableCornerCuttingCheck = true;
static const double kTieBreaker = 1.0 + 1e-4;

// =====================================================

Astarpath::~Astarpath() {
  if (data_raw)  delete[] data_raw;
  if (data_plan) delete[] data_plan;

  data_raw = nullptr;
  data_plan = nullptr;

  if (Map_Node) {
    for (int i = 0; i < GRID_X_SIZE; i++) {
      for (int j = 0; j < GRID_Y_SIZE; j++) {
        for (int k = 0; k < GRID_Z_SIZE; k++) {
          delete Map_Node[i][j][k];
        }
        delete[] Map_Node[i][j];
      }
      delete[] Map_Node[i];
    }
    delete[] Map_Node;
    Map_Node = nullptr;
  }
}

void Astarpath::begin_grid_map(double _resolution,
                               Vector3d global_xyz_l,
                               Vector3d global_xyz_u,
                               int max_x_id, int max_y_id, int max_z_id) {
  gl_xl = global_xyz_l(0);
  gl_yl = global_xyz_l(1);
  gl_zl = global_xyz_l(2);

  gl_xu = global_xyz_u(0);
  gl_yu = global_xyz_u(1);
  gl_zu = global_xyz_u(2);

  GRID_X_SIZE = max_x_id;
  GRID_Y_SIZE = max_y_id;
  GRID_Z_SIZE = max_z_id;
  GLYZ_SIZE = GRID_Y_SIZE * GRID_Z_SIZE;
  GLXYZ_SIZE = GRID_X_SIZE * GLYZ_SIZE;

  resolution = _resolution;
  inv_resolution = 1.0 / _resolution;

  data_raw  = new uint8_t[GLXYZ_SIZE];
  data_plan = new uint8_t[GLXYZ_SIZE];
  memset(data_raw,  0, GLXYZ_SIZE);
  memset(data_plan, 0, GLXYZ_SIZE);

  Map_Node = new MappingNodePtr**[GRID_X_SIZE];
  for (int i = 0; i < GRID_X_SIZE; i++) {
    Map_Node[i] = new MappingNodePtr*[GRID_Y_SIZE];
    for (int j = 0; j < GRID_Y_SIZE; j++) {
      Map_Node[i][j] = new MappingNodePtr[GRID_Z_SIZE];
      for (int k = 0; k < GRID_Z_SIZE; k++) {
        Vector3i idx(i, j, k);
        Map_Node[i][j][k] = new MappingNode(idx, gridIndex2coord(idx));
      }
    }
  }
}

void Astarpath::resetGrid(MappingNodePtr ptr) {
  ptr->id = 0;
  ptr->Father = nullptr;
  ptr->g_score = std::numeric_limits<double>::infinity();
  ptr->f_score = std::numeric_limits<double>::infinity();
}

void Astarpath::resetUsedGrids() {
  for (auto n : used_nodes) resetGrid(n);
  used_nodes.clear();
  while (!open_pq.empty()) open_pq.pop();
  terminatePtr = nullptr;
}
void Astarpath::resetOccupy() {
  memset(data_raw,  0, GLXYZ_SIZE);
  memset(data_plan, 0, GLXYZ_SIZE);
}

// ------------------- 坐标转换 -------------------

Vector3d Astarpath::gridIndex2coord(const Vector3i &index) {
  Vector3d pt;
  pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
  pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
  pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;
  return pt;
}

Vector3i Astarpath::coord2gridIndex(const Vector3d &pt) {
  Vector3i idx;
  idx << min(max(int((pt(0) - gl_xl) * inv_resolution), 0), GRID_X_SIZE - 1),
         min(max(int((pt(1) - gl_yl) * inv_resolution), 0), GRID_Y_SIZE - 1),
         min(max(int((pt(2) - gl_zl) * inv_resolution), 0), GRID_Z_SIZE - 1);
  return idx;
}

Vector3i Astarpath::c2i(const Vector3d &pt) { return coord2gridIndex(pt); }

Vector3d Astarpath::coordRounding(const Vector3d &coord) {
  return gridIndex2coord(coord2gridIndex(coord));
}


// ================= Occupancy =================

bool Astarpath::isOccupiedRaw(const int &idx_x, const int &idx_y,
                                        const int &idx_z) const {
  return (idx_x >= 0 && idx_x < GRID_X_SIZE && idx_y >= 0 && idx_y < GRID_Y_SIZE &&
          idx_z >= 0 && idx_z < GRID_Z_SIZE &&
          (data_raw[idx_x * GLYZ_SIZE + idx_y * GRID_Z_SIZE + idx_z] == 1));
}
inline bool Astarpath::isOccupied(const Eigen::Vector3i &idx) const {
  return isOccupied(idx(0), idx(1), idx(2));
}

bool Astarpath::is_occupy(const Vector3i &idx) {
  return isOccupiedRaw(idx(0), idx(1), idx(2));
}

inline bool Astarpath::isOccupiedRaw(const Eigen::Vector3i &idx) const {
  return isOccupiedRaw(idx(0), idx(1), idx(2));
}

bool Astarpath::isOccupied(int x, int y, int z) const {
  return (data_plan[x * GLYZ_SIZE + y * GRID_Z_SIZE + z] == 1);
}
// ------------------- 起飞竖直通道检查 -------------------
// 检查 (x,y) 固定时，从 z_from 到 z_to 的竖直柱状通道是否无障碍
bool Astarpath::verticalCorridorFree(const Eigen::Vector3i& xy_idx,
                                     int z_from, int z_to) const {
  int z0 = std::min(z_from, z_to);
  int z1 = std::max(z_from, z_to);

  int x = xy_idx(0);
  int y = xy_idx(1);

  if (x < 0 || x >= GRID_X_SIZE || y < 0 || y >= GRID_Y_SIZE) return false;

  z0 = std::max(0, z0);
  z1 = std::min(GRID_Z_SIZE - 1, z1);

  for (int z = z0; z <= z1; ++z) {
    if (isOccupied(x, y, z)) return false;
  }
  return true;
}

// ------------------- 障碍设置-------------------
// 关键：固定高度规划下，障碍膨胀主要做 XY，Z 方向不要跟 r 一起膨胀
void Astarpath::set_barrier(double x, double y, double z, int inflation_r) {
  if (x < gl_xl || y < gl_yl || z < gl_zl ||
      x >= gl_xu || y >= gl_yu || z >= gl_zu)
    return;

  int ix = (x - gl_xl) * inv_resolution;
  int iy = (y - gl_yl) * inv_resolution;
  int iz = (z - gl_zl) * inv_resolution;

  if (ix < 0 || iy < 0 || iz < 0 ||
      ix >= GRID_X_SIZE || iy >= GRID_Y_SIZE || iz >= GRID_Z_SIZE)
    return;

  // 真实障碍（不膨胀）
  data_raw[ix * GLYZ_SIZE + iy * GRID_Z_SIZE + iz] = 1;

  // 规划障碍（膨胀）
  int r = max(0, inflation_r);
  for (int dx = -r; dx <= r; dx++)
    for (int dy = -r; dy <= r; dy++)
      for (int dz = -kInflateZLayers; dz <= kInflateZLayers; dz++) {
        int nx = ix + dx, ny = iy + dy, nz = iz + dz;
        if (nx >= 0 && ny >= 0 && nz >= 0 &&
            nx < GRID_X_SIZE && ny < GRID_Y_SIZE && nz < GRID_Z_SIZE) {
          data_plan[nx * GLYZ_SIZE + ny * GRID_Z_SIZE + nz] = 1;
        }
      }
}

// ------------------- 防穿角-------------------
bool Astarpath::validMove(const Vector3i& from, const Vector3i& to) const {
  if (isOccupied(to)) return false;
  if (!kEnableCornerCuttingCheck) return true;

  Vector3i d = to - from;
  int adx = std::abs(d(0)), ady = std::abs(d(1)), adz = std::abs(d(2));
  int s = adx + ady + adz;
  if (s <= 1) return true;

  if (d(0) != 0) {
    Vector3i mid = from + Vector3i(d(0), 0, 0);
    if (isOccupied(mid)) return false;
  }
  if (d(1) != 0) {
    Vector3i mid = from + Vector3i(0, d(1), 0);
    if (isOccupied(mid)) return false;
  }
  if (d(2) != 0) {
    Vector3i mid = from + Vector3i(0, 0, d(2));
    if (isOccupied(mid)) return false;
  }
  return true;
}

// ------------------- 取邻居：固定高度 2D A*-------------------
void Astarpath::AstarGetSucc(MappingNodePtr currentPtr,
                             vector<MappingNodePtr> &neighborPtrSets,
                             vector<double> &edgeCostSets) {
  neighborPtrSets.clear();
  edgeCostSets.clear();

  const Vector3i cur = currentPtr->index;
  const int z = cur(2);  // 固定高度层

  for (int dx = -1; dx <= 1; dx++) {
    for (int dy = -1; dy <= 1; dy++) {
      if (dx == 0 && dy == 0) continue;

      Vector3i nb(cur(0) + dx, cur(1) + dy, z);

      if (nb(0) < 0 || nb(0) >= GRID_X_SIZE ||
          nb(1) < 0 || nb(1) >= GRID_Y_SIZE)
        continue;

      // dz=0，所以 validMove 不会触发 z 方向的 mid-check（仍可保留）
      if (!validMove(cur, nb)) continue;

      neighborPtrSets.push_back(Map_Node[nb(0)][nb(1)][nb(2)]);

      // 2D 步长
      double step = std::sqrt(double(dx*dx + dy*dy));
      edgeCostSets.push_back(step * resolution);
    }
  }
}

// ------------------- 启发式-------------------
double Astarpath::getHeu(MappingNodePtr node1, MappingNodePtr node2) {
  Vector3d delta = node2->coord - node1->coord;
  return delta.norm() * kTieBreaker;
}

// ------------------- 软避障代价-------------------
// 思路：只看“最近障碍距离”，代价 ∈ [0, kObsCostWeight]
double Astarpath::obstacleProximityCost(const Eigen::Vector3i& idx) const {
  const int z = idx(2);  // 固定高度层
  int best_d2 = 1e9;

  for (int dx = -kObsCostRadius; dx <= kObsCostRadius; ++dx) {
    for (int dy = -kObsCostRadius; dy <= kObsCostRadius; ++dy) {
      if (dx == 0 && dy == 0) continue;

      int x = idx(0) + dx;
      int y = idx(1) + dy;
      if (x < 0 || x >= GRID_X_SIZE || y < 0 || y >= GRID_Y_SIZE) continue;

      if (isOccupied(Eigen::Vector3i(x, y, z))) {
        int d2 = dx*dx + dy*dy;
        if (d2 < best_d2) best_d2 = d2;
      }
    }
  }

  if (best_d2 == 1e9) return 0.0;

  double d = std::sqrt((double)best_d2); // grid distance
  double normalized = (double(kObsCostRadius) - d) / double(kObsCostRadius);
  if (normalized < 0) normalized = 0;
  if (normalized > 1) normalized = 1;

  return kObsCostWeight * normalized;
}

// ------------------- A* 主函数：先起飞到规划高度，再平面搜索 -------------------
bool Astarpath::AstarSearch(Vector3d start_pt, Vector3d end_pt) {
  ros::Time time_1 = ros::Time::now();

  resetUsedGrids();

  // 1) 原始栅格索引（包含当前高度）
  Vector3i start_idx_raw = coord2gridIndex(start_pt);
  Vector3i end_idx_raw   = coord2gridIndex(end_pt);

  // 2) 规划高度层：默认取 end_pt 的高度（巡航高度）
  int z_plan_idx = end_idx_raw(2);

  // 3) 起飞可行性检查：在 (start_x, start_y) 竖直上升到 z_plan_idx 不能撞
  if (!verticalCorridorFree(start_idx_raw, start_idx_raw(2), z_plan_idx)) {
    ROS_ERROR("[A*] Takeoff corridor blocked: cannot climb to planning altitude!");
    return false;
  }

  // 4) 固定高度规划：start/goal 都强制在 z_plan_idx
  Vector3i start_idx = start_idx_raw;
  Vector3i end_idx   = end_idx_raw;
  start_idx(2) = z_plan_idx;
  end_idx(2)   = z_plan_idx;
  goalIdx = end_idx;

  // 5) 对齐到格心（用于节点坐标）
  start_pt = gridIndex2coord(start_idx);
  end_pt   = gridIndex2coord(end_idx);

  // 6) 起点终点占据检查（在规划高度层）
  if (isOccupied(start_idx)) {
    ROS_ERROR("[A*] Start position at planning altitude is occupied!");
    return false;
  }
  if (isOccupied(end_idx)) {
    ROS_ERROR("[A*] Goal position at planning altitude is occupied!");
    return false;
  }

  MappingNodePtr startPtr = Map_Node[start_idx(0)][start_idx(1)][start_idx(2)];
  MappingNodePtr endPtr   = Map_Node[end_idx(0)][end_idx(1)][end_idx(2)];

  // init start
  touchNode(startPtr);
  startPtr->g_score = 0.0;
  startPtr->f_score = getHeu(startPtr, endPtr);
  startPtr->Father  = nullptr;
  startPtr->id      = 1;

  open_pq.push({startPtr->f_score, startPtr});

  MappingNodePtr currentPtr = nullptr;
  vector<MappingNodePtr> neighborPtrSets;
  vector<double> edgeCostSets;

  int expanded_nodes = 0;

  while (!open_pq.empty()) {
    auto top = open_pq.top();
    open_pq.pop();

    currentPtr = top.node;
    double f_in_queue = top.f;

    // lazy deletion：过滤旧条目
    if (currentPtr->id == -1) continue;
    if (std::fabs(f_in_queue - currentPtr->f_score) > 1e-9) continue;

    // close
    currentPtr->id = -1;
    expanded_nodes++;

    // goal test
    if (currentPtr->index == end_idx) {
      ros::Time time_2 = ros::Time::now();
      terminatePtr = currentPtr;

      ROS_INFO("[A*] Path found successfully!");
      ROS_INFO("  - Time: %.2f ms", (time_2 - time_1).toSec() * 1000.0);
      ROS_INFO("  - Path cost: %.2f m", currentPtr->g_score);
      ROS_INFO("  - Expanded nodes: %d", expanded_nodes);
      return true;
    }

    // expand neighbors
    AstarGetSucc(currentPtr, neighborPtrSets, edgeCostSets);

    for (size_t i = 0; i < neighborPtrSets.size(); i++) {
      MappingNodePtr nb = neighborPtrSets[i];
      double edge_cost = edgeCostSets[i];

      if (nb->id == -1) continue;

      // 软避障：有上界，不会在窄通道里“把路炸没”
      double obs_cost = obstacleProximityCost(nb->index);

      double tentative_g = currentPtr->g_score + edge_cost + obs_cost;

      if (nb->id == 0) {
        touchNode(nb);
        nb->g_score = tentative_g;
        nb->f_score = tentative_g + getHeu(nb, endPtr);
        nb->Father  = currentPtr;
        nb->id      = 1;
        open_pq.push({nb->f_score, nb});
      } else if (nb->id == 1) {
        if (tentative_g < nb->g_score) {
          nb->g_score = tentative_g;
          nb->f_score = tentative_g + getHeu(nb, endPtr);
          nb->Father  = currentPtr;
          open_pq.push({nb->f_score, nb});
        }
      }
    }
  }

  ros::Time time_2 = ros::Time::now();
  ROS_ERROR("[A*] Path NOT found!");
  ROS_ERROR("  - Time: %.2f ms", (time_2 - time_1).toSec() * 1000.0);
  ROS_ERROR("  - Expanded nodes: %d", expanded_nodes);
  return false;
}

// ------------------- 输出路径/访问节点 -------------------

vector<Vector3d> Astarpath::getPath() {
  vector<Vector3d> path;
  if (terminatePtr == nullptr) {
    ROS_WARN("[A*] terminatePtr is NULL, no path available!");
    return path;
  }

  vector<MappingNodePtr> temp;
  MappingNodePtr cur = terminatePtr;
  while (cur != nullptr) {
    temp.push_back(cur);
    cur = cur->Father;
  }

  for (int i = (int)temp.size() - 1; i >= 0; i--) {
    path.push_back(temp[i]->coord);
  }

  ROS_INFO("[A*] Path extracted: %d waypoints", (int)path.size());
  return path;
}

vector<Vector3d> Astarpath::getVisitedNodes() {
  vector<Vector3d> visited_nodes;
  visited_nodes.reserve(used_nodes.size());
  for (auto n : used_nodes) {
    if (n->id == -1) visited_nodes.push_back(n->coord);
  }
  ROS_WARN("visited_nodes size : %d", (int)visited_nodes.size());
  return visited_nodes;
}

// ------------------- Douglas-Peucker 简化 -------------------

std::vector<Vector3d> Astarpath::pathSimplify(const vector<Vector3d> &path,
                                              double path_resolution) {
  if (path.size() < 3) return path;

  double dmax = 0;
  int index = 0;
  int end = (int)path.size() - 1;

  for (int i = 1; i < end; i++) {
    double d = perpendicularDistance(path[i], path[0], path[end]);
    if (d > dmax) {
      index = i;
      dmax = d;
    }
  }

  vector<Vector3d> resultPath;
  if (dmax > path_resolution) {
    vector<Vector3d> sub1(path.begin(), path.begin() + index + 1);
    vector<Vector3d> sub2(path.begin() + index, path.end());

    vector<Vector3d> rec1 = pathSimplify(sub1, path_resolution);
    vector<Vector3d> rec2 = pathSimplify(sub2, path_resolution);

    resultPath.insert(resultPath.end(), rec1.begin(), rec1.end() - 1);
    resultPath.insert(resultPath.end(), rec2.begin(), rec2.end());
  } else {
    resultPath.push_back(path[0]);
    resultPath.push_back(path[end]);
  }
  return resultPath;
}

double Astarpath::perpendicularDistance(const Eigen::Vector3d point_insert,
                                        const Eigen::Vector3d point_st,
                                        const Eigen::Vector3d point_end) {
  Vector3d line1 = point_end - point_st;
  Vector3d line2 = point_insert - point_st;
  return double(line2.cross(line1).norm() / line1.norm());
}

// ------------------- 多项式相关（保留你原逻辑）-------------------

Vector3d Astarpath::getPosPoly(MatrixXd polyCoeff, int k, double t) {
  Vector3d ret;
  int poly_num1D = (int)polyCoeff.cols() / 3;

  for (int dim = 0; dim < 3; dim++) {
    VectorXd coeff = (polyCoeff.row(k)).segment(dim * poly_num1D, poly_num1D);
    VectorXd time = VectorXd::Zero(poly_num1D);

    for (int j = 0; j < poly_num1D; j++) {
      time(j) = (j == 0) ? 1.0 : std::pow(t, j);
    }
    ret(dim) = coeff.dot(time);
  }
  return ret;
}

int Astarpath::safeCheck(MatrixXd polyCoeff, VectorXd time) {
  int unsafe_segment = -1;
  double delta_t = resolution;
  double t = delta_t;
  Vector3d pos;

  for (int i = 0; i < polyCoeff.rows(); i++) {
    while (t < time(i)) {
      pos = getPosPoly(polyCoeff, i, t);
      Vector3i idx = coord2gridIndex(pos);
      if (isOccupiedRaw(idx(0), idx(1), idx(2))) {
        unsafe_segment = i;
        return unsafe_segment;
      }
      t += delta_t;
    }
    t = delta_t;
  }
  return -1;
}

