#pragma once

#include <chrono>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace kal_a_star {

using Duration = std::chrono::duration<double, std::ratio<1>>;
using Time = std::chrono::time_point<std::chrono::steady_clock, Duration>;

using Position = Eigen::Vector2d;
using Path = std::vector<Position>;
struct StampedPosition {
    Position position;
    Time stamp;
};
using Trajectory = std::vector<StampedPosition>;
using Pose = Eigen::Isometry2d;

struct Point {
    double x;
    double y;
};


class Node {
public:
    double x;            // x坐标
    double y;            // y坐标
    double gCost;        // 从起始节点到当前节点的路径代价
    double hCost;        // 从当前节点到目标节点的估算代价
    Node* parent;        // 父节点指针
    Node(double x, double y) : x(x), y(y), gCost(0), hCost(0), parent(nullptr) {}
    // 计算总代价
    double getFCost() const {
        return gCost + hCost;
    }
};

} // namespace a_star