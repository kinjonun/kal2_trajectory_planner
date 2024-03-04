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
    double x;
    double y;
    double gCost;        // Path cost from the start node to the current node
    double hCost;        // Estimated cost from current node to target node
    Node* parent;        // Parent Node Pointer
    Node(double x, double y) : x(x), y(y), gCost(0), hCost(0), parent(nullptr) {}
    // Calculation of the total cost
    double getFCost() const {
        return gCost + hCost;
    }

    double  getX() const {
        return x;
    }
};

} // namespace a_star