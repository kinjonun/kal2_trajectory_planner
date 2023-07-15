#pragma once

#include "types.hpp"

namespace kal_a_star {
    const std::vector<Point> conesright = {{0.8, 0.93}, {1.0, 0.95},{1.1, 1.0}, {1.2, 1.1}, {1.4, 1.2}, {1.5, 1.3},{1.6, 1.5}, {1.8, 1.7},
                                {2.0,1.9},{2.2, 1.95}, {2.4, 1.8}, {2.5,2.25}, {2.6, 1.7}, {2.8, 1.6},{3, 1.5} };
    const std::vector<Point> coneslinks = {{1.0, 1.5}, {1.2, 1.7}, {1.4, 2.0},{1.6, 2.1}, {1.8, 2.2}, {2.0, 2.4},
                               {2.2, 2.5}, {2.4, 2.5}, {2.6, 2.2}, {2.8, 2.1}, {3.0, 2.0},{3.2, 1.85} };

    extern std::vector<Point> linkcones;
    extern std::vector<Point> rightcones;

    extern Point targetpoint;
    extern Point startpoint;

    const Point currentposition = {0, 0.8};

    const Node startNode(1.0, 1.25);
    const Node targetNode(3.0, 1.7); // startNode, targetNode, right, link

    const double MAP_WIDTH = 15;
    const double MAP_HEIGHT = 15;
    const double radius = 0.15;
    const double distanztocenter = 0.05;

    const int COST_FREE = 1;       // 自由路径代价
    const double COST_COLLISION = 2; // 碰撞代价
    const double COST_DEVIATION = 200;  // 偏离中心线代价



 class A_star
    {
    private:
        std::vector<Point> conesposition_;

    public:

        void cones(const std::vector<Point>& cones);

        std::vector<Node*> aStarSearch(const Node& startNode, const Node& targetNode, const std::vector<Point>& link, const std::vector<Point>& right );
        std::vector<Node*> aStar(const Pose currentPosition);
        const std::vector<Point>& getConesPosition() const { return conesposition_; }

        // double distanceToPoint(const Point& point1, const Point& point2);
        // double sumDistancesToTarget(const std::vector<Point>& coords, const Point& target);
        // void printPath(const std::vector<Node*>& path);
        // void print_vector(const std::vector<Node*>& path);
        // 计算从当前节点到目标节点的欧几里得距离作为估算代价
        // double calculateHCost(double currentX, double currentY, double targetX, double targetY);
        // bool isObstacle( double x,  double y, const std::vector<Point> &coords);
        // Eigen::VectorXd fitPolynomial(const std::vector<Point>& points, int degree);
        // bool isWithinMap(double x, double y);
        // bool isNodeInList(const std::vector<Node*> &nodeList, double x, double y);
        // Node* getMinFCostNode(const std::vector <Node*>&nodeList) ;
        // std::vector<Node*> buildPath(Node* endNode);
        // double polynomial_center(const Eigen::VectorXd polynomialParameters, double x);
        // std::vector<Position> convertToEigenVector(const std::vector<Node*>& path);
        // void linkrightcones(std::vector<Point>& ConePoints, const Point& currentPosition, std::vector<Point>& linkcones, std::vector<Point>& rightcones);
        // void startandtraget(const std::vector<Point>& linkcones, const std::vector<Point>& rightcones,  Point& startpoint, Point& targetpoint);

    };



}
