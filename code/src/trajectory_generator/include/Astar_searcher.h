#pragma once

#include <ros/ros.h>
#include <Eigen/Eigen>
#include <vector>
#include <queue>
#include <cstdint>
#include <cstring>
#include <limits>
#include <cmath>

class MappingNode;
typedef MappingNode* MappingNodePtr;

class MappingNode {
public:
  Eigen::Vector3i index;
  Eigen::Vector3d coord;

  // A* states: 0 = unvisited, 1 = in open, -1 = in closed
  int id = 0;
  double g_score = std::numeric_limits<double>::infinity(); // meters
  double f_score = std::numeric_limits<double>::infinity(); // meters
  MappingNodePtr Father = nullptr;

  MappingNode(const Eigen::Vector3i& idx, const Eigen::Vector3d& pos)
      : index(idx), coord(pos) {}
};

class Astarpath {
public:
  Astarpath() = default;
  ~Astarpath();

  void begin_grid_map(double _resolution,
                      Eigen::Vector3d global_xyz_l,
                      Eigen::Vector3d global_xyz_u,
                      int max_x_id, int max_y_id, int max_z_id);

  void resetOccupy();
  void set_barrier(double coord_x, double coord_y, double coord_z, int inflation_r = 1);

  bool AstarSearch(Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);

  std::vector<Eigen::Vector3d> getPath();
  std::vector<Eigen::Vector3d> getVisitedNodes();
  void resetUsedGrids();

  // optional: simplify
  std::vector<Eigen::Vector3d> pathSimplify(const std::vector<Eigen::Vector3d> &path,
                                            double path_resolution);
  double perpendicularDistance(const Eigen::Vector3d point_insert,
                               const Eigen::Vector3d point_st,
                               const Eigen::Vector3d point_end);

  // poly related (kept from your code)
  Eigen::Vector3d getPosPoly(Eigen::MatrixXd polyCoeff, int k, double t);
  int safeCheck(Eigen::MatrixXd polyCoeff, Eigen::VectorXd time);

  // coordinate helpers
  Eigen::Vector3d gridIndex2coord(const Eigen::Vector3i &index);
  Eigen::Vector3i coord2gridIndex(const Eigen::Vector3d &pt);
  Eigen::Vector3i c2i(const Eigen::Vector3d &pt);
  Eigen::Vector3d coordRounding(const Eigen::Vector3d &coord);

  // occupancy queries
  bool is_occupy(const Eigen::Vector3i &index);
  inline bool isOccupied(const Eigen::Vector3i &idx) const;
  inline bool isOccupiedRaw(const Eigen::Vector3i &idx) const;

  //inline bool isOccupied(const Eigen::Vector3i &index) const;
  //inline bool isFree(const Eigen::Vector3i &index) const;

private:
  // map bounds
  double gl_xl = 0, gl_yl = 0, gl_zl = 0;
  double gl_xu = 0, gl_yu = 0, gl_zu = 0;

  // grid
  int GRID_X_SIZE = 0, GRID_Y_SIZE = 0, GRID_Z_SIZE = 0;
  int GLYZ_SIZE = 0, GLXYZ_SIZE = 0;

  // resolution
  double resolution = 1.0;
  double inv_resolution = 1.0;
  

  // occupancy data
  uint8_t* data_raw = nullptr;   // 真实障碍：给 issafe / safeCheck
  uint8_t* data_plan = nullptr;  // 规划障碍：给 A*


  // nodes
  MappingNodePtr*** Map_Node = nullptr;

  // goal
  Eigen::Vector3i goalIdx;
  MappingNodePtr terminatePtr = nullptr;

  // for fast reset
  std::vector<MappingNodePtr> used_nodes;

  // --- A* open set (priority queue) ---
  struct OpenEntry {
    double f;
    MappingNodePtr node;
  };
  struct OpenCmp {
    bool operator()(const OpenEntry& a, const OpenEntry& b) const {
      // min-heap behavior with priority_queue by reversing
      return a.f > b.f;
    }
  };
  std::priority_queue<OpenEntry, std::vector<OpenEntry>, OpenCmp> open_pq;

private:
  void resetGrid(MappingNodePtr ptr);
  

  bool isOccupied(int idx_x, int idx_y, int idx_z) const;
  //inline bool isFree(int idx_x, int idx_y, int idx_z) const;
  double obstacleProximityCost(const Eigen::Vector3i& idx) const;
  bool verticalCorridorFree(const Eigen::Vector3i& xy_idx, int z_from, int z_to) const;
  bool isOccupiedRaw(const int &idx_x, const int &idx_y, const int &idx_z) const;




  // successors
  void AstarGetSucc(MappingNodePtr currentPtr,
                    std::vector<MappingNodePtr> &neighborPtrSets,
                    std::vector<double> &edgeCostSets);

  // real heuristic in meters
  double getHeu(MappingNodePtr node1, MappingNodePtr node2);

  // optional: prevent corner-cutting in 26-neighborhood
  bool validMove(const Eigen::Vector3i& from, const Eigen::Vector3i& to) const;

  inline void touchNode(MappingNodePtr n) {
    if (n->id == 0) used_nodes.push_back(n);
  }
};

