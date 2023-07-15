#pragma once

#include <cmath>
#include <chrono>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <glog/logging.h>


#include "../types.hpp"

namespace kal_a_star::utils {


inline double distanceToPoint(const Point& point1, const Point& point2) {
    double x1 = point1.x;
    double y1 = point1.y;
    double x2 = point2.x;
    double y2 = point2.y;
    return std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}

inline double sumDistancesToTarget(const std::vector<Point>& coords, const Point& target) {
    double sum = 0.0;
    for (const auto& point : coords) {
        double dist = distanceToPoint(target, point);
        sum += (1/(dist));
    }

    return sum;
}

// 检查节点是否为障碍物
inline bool isObstacle( double x,  double y, const std::vector<Point> &coords) {
    Point target1 = {x,y};
    std::vector<Point> obstacle = coords;
    for(const auto& point : obstacle) {
        double dist = distanceToPoint(target1, point);
        if (dist <= radius)
        {
            return true;
        }
    }
    return false;
}


inline double polynomial_center(const Eigen::VectorXd polynomialParameters, double x) {
    int64_t polynomialDegree = polynomialParameters.rows() - 1;
    double y = 0;
    for (int i = polynomialDegree; i >= 0; i--) {
        y += polynomialParameters(i) * std::pow(x, i);
    }
    return y;
}


inline double calculateDerivative(const Eigen::VectorXd polynomialParameters, double x) {
    int64_t polynomialDegree = polynomialParameters.rows() - 1;
    double y = 0;
    for (int i = polynomialDegree; i >= 1; i--) {
        y += polynomialParameters(i) * std::pow(x, i-1) * polynomialDegree;
    }
    return y;
}


// 计算从当前节点到目标节点的欧几里得距离作为估算代价
inline double calculateHCost(double currentX, double currentY, double targetX, double targetY) {
    double dx = std::abs(currentX - targetX);
    double dy = std::abs(currentY - targetY);
    return static_cast<double>(std::sqrt(dx * dx + dy * dy));
}

// 拟合多项式的函数
inline Eigen::VectorXd fitPolynomial(const std::vector<Point>& points, int degree) {

    int numPoints = points.size();
    Eigen::MatrixXd A(numPoints, degree + 1);
    Eigen::VectorXd y(numPoints);
    for (int i = 0; i < numPoints; i++) {
        double xi = points[i].x;
        double yi = points[i].y;
        for (int j = 0; j <= degree; j++) {
            A(i, j) = std::pow(xi, j);
        }
        y(i) = yi;
    }

    // 使用最小二乘法拟合多项式曲线
    Eigen::VectorXd coeffs = A.colPivHouseholderQr().solve(y);
    return coeffs;
}
// 检查节点是否在地图范围内
inline bool isWithinMap(double x, double y) {
    return (x >= -10 && x < MAP_WIDTH && y >= -10 && y < MAP_HEIGHT);
}


// 检查节点是否已经在开放列表或封闭列表中
inline bool isNodeInList(const std::vector<Node*> &nodeList, double x, double y) {
    for (const Node* node : nodeList) {
        if (node->x == x && node->y == y) {
            return true;
        }
    }
    return false;
}

// 获取开放列表中代价最小的节点
inline Node* getMinFCostNode(const std::vector <Node*>&nodeList)  {
    Node* minNode = nodeList[0];
    for (const Node* node : nodeList) {
        if (node->getFCost() < minNode->getFCost()) {
            minNode = const_cast<Node*>(node);
        }
    }
    return minNode;
}

// 构建路径
inline std::vector<Node*> buildPath(Node* endNode) {
    std::vector<Node*> path;
    Node* currentNode = endNode;
    while (currentNode != nullptr) {
        path.push_back(currentNode);
        currentNode = currentNode->parent;
    }
    std::reverse(path.begin(), path.end());
    return path;
}


inline bool compareByDistance(const Point& point1, const Point& point2, const Point& currentPosition) {
    double distance1 = distanceToPoint(point1, currentPosition);
    double distance2 = distanceToPoint(point2, currentPosition);
    return distance1 < distance2;
}


inline double angleBetweenPoints(const Point& currentPosition, const Point& point1, const Point& point2) {
    double dx1 = point1.x - currentPosition.x;
    double dy1 = point1.y - currentPosition.y;
    double dx2 = point2.x - currentPosition.x;
    double dy2 = point2.y - currentPosition.y;


    // 计算两个向量的内积
    double dotProduct = dx1 * dx2 + dy1 * dy2;

    double magnitude1 = std::sqrt(dx1 * dx1 + dy1 * dy1);
    double magnitude2 = std::sqrt(dx2 * dx2 + dy2 * dy2);

    // 计算夹角的余弦值
    double cosAngle = dotProduct / (magnitude1 * magnitude2);

    // 使用反余弦函数得到夹角的弧度
    double angleRad = std::acos(cosAngle);

    // 将弧度转换为角度并返回
    double angleDeg = angleRad * 180.0 / M_PI;
    return angleDeg;
}


inline void linkrightcones( std::vector<Point>& ConePoints, const Point& currentPosition,std::vector<Point>& linkcones,std::vector<Point>& rightcones){


    std::sort(ConePoints.begin(), ConePoints.end(), [&currentPosition](const Point& p1, const Point& p2) {
        return compareByDistance(p1, p2, currentPosition);
    });

    std::vector<Point>::iterator it;

    if (ConePoints[0].x >0) {
        rightcones.push_back(ConePoints[0]);
        for (it = ConePoints.begin() + 1; it != ConePoints.end(); ++it){
            if (it->x < 0){
                linkcones.push_back(*it);
                break;
            }
        }

    } else{
        linkcones.push_back(ConePoints[0]);
        for (it = ConePoints.begin() + 1; it != ConePoints.end(); ++it){
            if (it->x > 0){
                rightcones.push_back(*it);
                break;
            }
        }
    }

    std::vector<Point> destinationPoints;

    for (const auto& point : ConePoints) {
        if (linkcones.empty()){
            if((point.x == rightcones[0].x && point.y == rightcones[0].y)){
                    continue;
            }
        }else if (rightcones.empty()){
            if((point.x == linkcones[0].x && point.y == linkcones[0].y)){
                    continue;
            }
        }else{
            if ((point.x == linkcones[0].x && point.y == linkcones[0].y)||(point.x == rightcones[0].x && point.y == rightcones[0].y)||(point.x == 0 && point.y == 0))
                    continue;
        }

        destinationPoints.push_back(point);
    }


    // for (const Point& Position : destinationPoints) {
    //         std::cout << "destinationPoints (" << Position.x << ", " << Position.y << ")" << std::endl;
    //     }



    if(destinationPoints.size() < 2){
        throw std::runtime_error("destinationPoints ist klein________ _________ ");}

    if (linkcones.empty()&&rightcones.empty())
        { throw std::runtime_error("rightcones and linkcones must contain points!"); }

     for (const Point& point : destinationPoints) {
        if ((!linkcones.empty()) && (!rightcones.empty()) ){
            unsigned int  n = linkcones.size();
            unsigned int  m = rightcones.size();
            double b = angleBetweenPoints(currentPosition,linkcones[n-1],point);
            double c = angleBetweenPoints(currentPosition,rightcones[m-1],point);
            if(b < c){
                linkcones.push_back(point);
            } else{
                rightcones.push_back(point);
            }
        }

        if(linkcones.empty()){
            rightcones.push_back(point);
        }

        if(rightcones.empty()){
            linkcones.push_back(point);
        }

    }

    // std::cout << linkcones.size()<< " tmd"<<std::endl;
    // std::cout << rightcones.size()<< " tmd"<<std::endl;

    // for (const Point& Position : rightcones) {
    //         std::cout << "rightcones (" << Position.x << ", " << Position.y << ")" << std::endl;
    //     }
    // for (const Point& Position : linkcones) {
    //         std::cout << "linkcones (" << Position.x << ", " << Position.y << ")" << std::endl;
    //     }

}




inline void startandtraget(const std::vector<Point>& linkcones, const std::vector<Point>& rightcones,  Point& startpoint, Point& targetpoint) {

    unsigned int i = linkcones.size();
    unsigned int j = rightcones.size();

    if (i <= 2 && j <= 2) {
        throw std::runtime_error("cones数量太少 ");
        }

     else {

        if (i >= 2 && j >= 2) {
            startpoint = {(linkcones[0].x + rightcones[0].x) / 2, (linkcones[0].y + rightcones[0].y) / 2};
            targetpoint = {(linkcones[linkcones.size() - 1].x + rightcones[rightcones.size() - 1].x) / 2,
                           (linkcones[linkcones.size() - 1].y + rightcones[rightcones.size() - 1].y) / 2};

        } else if (j >= 2) {

            double dx = rightcones[rightcones.size() - 1].x - rightcones[rightcones.size() - 2].x;
            double dy = rightcones[rightcones.size() - 1].y - rightcones[rightcones.size() - 2].y;
            double nx = 0.4 * dy / sqrt(dx * dx + dy * dy);
            double ny = -0.4 * dx / sqrt(dx * dx + dy * dy);
            targetpoint = {rightcones[rightcones.size() - 1].x - nx, rightcones[rightcones.size() - 1].y - ny};

            if(linkcones.empty()){
                double dxr = rightcones[1].x - rightcones[0].x;
                double dyr = rightcones[1].y - rightcones[0].y;
                double nxr = 0.4 * dyr / sqrt(dxr * dxr + dyr * dyr);
                double nyr = -0.4 * dxr / sqrt(dxr * dxr + dyr * dyr);
                startpoint = {rightcones[0].x - nxr, rightcones[0].y - nyr};
            } else {
                startpoint = {(linkcones[0].x + rightcones[0].x) / 2, (linkcones[0].y + rightcones[0].y) / 2};
            }


        } else {

            double dx = linkcones[linkcones.size() - 1].x - linkcones[linkcones.size() - 2].x;
            double dy = linkcones[linkcones.size() - 1].y - linkcones[linkcones.size() - 2].y;
            double nx = 0.4 * dy / sqrt(dx * dx + dy * dy);
            double ny = -0.4 * dx / sqrt(dx * dx + dy * dy);
            targetpoint = {linkcones[linkcones.size() - 1].x + nx, linkcones[linkcones.size() - 1].y + ny};

            if(rightcones.empty()){
                double dxl = linkcones[1].x - linkcones[0].x;
                double dyl = linkcones[1].y - linkcones[0].y;
                double nxl = 0.4 * dyl / sqrt(dxl * dxl + dyl * dyl);
                double nyl = -0.4 * dxl / sqrt(dxl * dxl + dyl * dyl);
                startpoint = {linkcones[0].x + nxl, linkcones[0].y + nyl};
            } else {
                startpoint = {(linkcones[0].x + rightcones[0].x) / 2, (linkcones[0].y + rightcones[0].y) / 2};
            }
        }
    }

}


} //namespace kal_a_star::utils