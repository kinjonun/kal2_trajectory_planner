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

inline double conetoconesClosestDistance(const Point& cone, const std::vector<Point>& points) {
    double closestDistance = std::numeric_limits<double>::max();

    for (const Point& point : points) {
        double distance = distanceToPoint(cone,point);
        if (distance < closestDistance) {
            closestDistance = distance;
        }
    }

    return closestDistance;
}

inline double maxDistanceToPoints(const Point& cone, const std::vector<Point>& points) {
    double Maxdist = 0.0;

    for (const Point& point : points ) {
        double distance = distanceToPoint(cone, point);
        if (distance > Maxdist) {
            Maxdist = distance;
        }
    }

    return Maxdist;
}


inline double sumDistancesToTarget(const std::vector<Point>& coords, const Point& target) {
    double sum = 0.0;
    for (const auto& point : coords) {
        double dist = distanceToPoint(target, point);
        sum += (1/(dist));
    }

    return sum;
}

// Check if the node is an obstacle
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


// Calculate the Euclidean distance from the current node to the target node as the estimated cost
inline double calculateHCost(double currentX, double currentY, double targetX, double targetY) {
    double dx = std::abs(currentX - targetX);
    double dy = std::abs(currentY - targetY);
    return static_cast<double>(std::sqrt(dx * dx + dy * dy));
}

// Functions that fit polynomials
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

    // Fitting polynomial curves using the least squares method
    Eigen::VectorXd coeffs = A.colPivHouseholderQr().solve(y);
    return coeffs;
}
// Check if the node is within the map
inline bool isWithinMap(double x, double y, const Point& start, const Point& traget ) {
    return (x >= std::min(start.x, traget.x) && x < std::max(start.x, traget.x) && y >= std::min(start.y, traget.y) && y < std::max(start.y, traget.y));
}


// Check if the node is already in the open or closed list
inline bool isNodeInList(const std::vector<Node*> &nodeList, double x, double y) {
    for (const Node* node : nodeList) {
        if (node->x == x && node->y == y) {
            return true;
        }
    }
    return false;
}

// Get the least costly node in the open list
inline Node* getMinFCostNode(const std::vector <Node*>&nodeList)  {
    Node* minNode = nodeList[0];
    for (const Node* node : nodeList) {
        if (node->getFCost() < minNode->getFCost()) {
            minNode = const_cast<Node*>(node);
        }
    }
    return minNode;
}

// Building paths
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

inline void removeAllZeroPoints(std::vector<Point>& ConePoints) {
    ConePoints.erase(
        std::remove_if(ConePoints.begin(), ConePoints.end(), [](const Point& p) {
            return p.x == 0 && p.y == 0;
        }),
        ConePoints.end()
    );
}


inline void print_position(const std::vector<Position>& position){
    for (auto point : position){
        std::cout << "( " << point[0] << ", " << point[1] << ")  ";
    }
    std::cout  << std::endl;
}

inline double angleBetweenPoints(const Point& currentPosition, const Point& point1, const Point& point2) {
    double dx1 = point1.x - currentPosition.x;
    double dy1 = point1.y - currentPosition.y;
    double dx2 = point2.x - currentPosition.x;
    double dy2 = point2.y - currentPosition.y;


    // Compute the inner product of two vectors
    double dotProduct = dx1 * dx2 + dy1 * dy2;

    double magnitude1 = std::sqrt(dx1 * dx1 + dy1 * dy1);
    double magnitude2 = std::sqrt(dx2 * dx2 + dy2 * dy2);

    // Calculate the cosine of the angle
    double cosAngle = dotProduct / (magnitude1 * magnitude2);

    // Use the inverse cosine function to get the radians of the angle of intersection
    double angleRad = std::acos(cosAngle);

    // Converts radians to angles and return
    double angleDeg = angleRad * 180.0 / M_PI;
    return angleDeg;
}

inline void conessortieren(std::vector<Point>& ConePoints, std::vector<Point>& erstreihe, std::vector<Point>& zweireihe){
    if(ConePoints.size() < 1){
        std::cout << "conesanzahl ist klein als ________ 1,kann nicht in Reihe hinzufügen " <<std::endl;
        return;
    }

// in welt koord. müssen

    removeAllZeroPoints(ConePoints);

    if (erstreihe.size()==0){
        erstreihe.push_back(ConePoints[0]);
    }

    Point bezugpoint = {erstreihe[0].x , erstreihe[0].y};
    std::sort(erstreihe.begin(), erstreihe.end(), [&bezugpoint](const Point& p1, const Point& p2) {
        return utils::compareByDistance(p1, p2, bezugpoint);
    });

    if (!zweireihe.empty()){
        std::sort(zweireihe.begin(), zweireihe.end(), [&bezugpoint](const Point& p1, const Point& p2) {
            return utils::compareByDistance(p1, p2, bezugpoint);
        });
    }



    for( auto point : ConePoints ){
        double m = conetoconesClosestDistance(point, erstreihe);

        if(zweireihe.size() == 0){
            if(m < 0.7 && m > 0.2){
                erstreihe.push_back(point);
            }else if( m > 0.7){
                zweireihe.push_back(point);
            }
        }
        else{

            double n = conetoconesClosestDistance(point, zweireihe);
            // 第一个不小于的
            auto it = std::lower_bound(erstreihe.begin(), erstreihe.end(), point, [](const Point& p1, const Point& p2) {
                return p1.x < p2.x;
            });

            auto it2 = std::lower_bound(zweireihe.begin(), zweireihe.end(), point, [](const Point& p1, const Point& p2) {
                return p1.x < p2.x;
            });

            double erstwinkel;
            double zweiwinkel;

            if (n > m && m >0.2 && m < 0.8 ) {
                // if (erstreihe.size() > 2) {
                //     if (it != erstreihe.begin() && it != erstreihe.end()) {
                //         Point p1 = *(it - 1);
                //         Point p2 = *it;
                //         erstwinkel = angleBetweenPoints(point, p1, p2);
                //         if (erstwinkel > 90) {
                //             erstreihe.insert(it, point);
                //         }

                //     } else if (it == erstreihe.begin()) {
                //         Point p1 = *(it + 1);
                //         Point p2 = *it;
                //         erstwinkel = angleBetweenPoints(point, p1, p2);
                //         if (erstwinkel > 90) {
                //             erstreihe.insert(it, point);
                //         }

                //     } else {
                //         Point p1 = *(it - 1);
                //         Point p2 = *(it - 2);
                //         erstwinkel = angleBetweenPoints(point, p1, p2);
                //         if (erstwinkel > 90) {
                //             erstreihe.push_back(point);
                //            }
                //     }


                // } else {
                    erstreihe.push_back(point);
               // }



            }else if( n < m && n > 0.2 && n < 0.8 && m > 0.6 && m < 1.5 ){
                // if (zweireihe.size() >2){
                //     if (it != zweireihe.begin() && it != zweireihe.end()) {
                //         Point p1 = *(it - 1);
                //         Point p2 = *it;
                //         zweiwinkel = angleBetweenPoints(point, p1, p2);
                //         if (erstwinkel > 100) {
                //             erstreihe.insert(it, point);
                //         }
                //     } else if (it == zweireihe.begin()) {
                //         Point p1 = *(it + 1);
                //         Point p2 = *it;
                //         zweiwinkel = angleBetweenPoints(point, p1, p2);
                //         if (erstwinkel > 100) {
                //             erstreihe.insert(it, point);
                //         }
                //     } else {
                //         Point p1 = *(it - 1);
                //         Point p2 = *(it - 2);
                //         zweiwinkel = angleBetweenPoints(point, p1, p2);
                //         erstreihe.push_back(point);
                //     }

                // } else{
                    zweireihe.push_back(point);
              //  }

            }
        }
    }

    for (const Point& Position : erstreihe) {
            std::cout << "erstreihe (" << Position.x << ", " << Position.y << ")" << std::endl;
        }
    for (const Point& Position : zweireihe) {
            std::cout << "zweireihe (" << Position.x << ", " << Position.y << ")" << std::endl;
        }
}


inline void calculateMittelpoint(std::vector<Point>& mittelpoints, std::vector<Point>& erstreihe, std::vector<Point>& zweireihe){
    if (zweireihe.empty() || erstreihe.empty()){
        return;
    }
    std::cout << "erst.size " << erstreihe.size() << "  zweireihe   " << zweireihe.size() << std::endl;
    mittelpoints.clear();
    for (int i =0; i < zweireihe.size(); i++){
        for (int j =0; j < erstreihe.size(); j++ ){
            mittelpoints.push_back( Point{ (zweireihe[i].x +erstreihe[j].x)/2, (zweireihe[i].y + erstreihe[j].y)/2 });
        }
    }
    Point currentPosition = {erstreihe[0].x,erstreihe[0].y};
    std::sort(mittelpoints.begin(), mittelpoints.end(), [&currentPosition](const Point& p1, const Point& p2) {
        return compareByDistance(p1, p2, currentPosition);
    });


}

inline int getSign(double value){
    return (value < 0)? -1:1;
}

// inline void startandzielpoint(const std::vector<Point>& erstreihe, const std::vector<Point>& zweireihe,  Point& startpoint, Point& targetpoint){
//     unsigned int i = erstreihe.size();
//     unsigned int j = zweireihe.size();

//     if (i == 0 || j == 0 ) {
//         std::cout << "erstreihe oder zweireihe ist leer---------------------------------- "<<std::endl<<std::endl;
//         return;
//     }
//     if (i > 1 && j > 1 ) {

//         startpoint = {(erstreihe[0].x + zweireihe[0].x) / 2, (erstreihe[0].y + zweireihe[0].y) / 2};
//         auto resulterstreihe = maxDistanceToPoints(startpoint, erstreihe);
//         double m =resulterstreihe.first;
//         std::vector<Point>::const_iterator MaxdistIterst = resulterstreihe.second;

//         auto resultzweireihe = maxDistanceToPoints(startpoint, zweireihe);
//         double n = resultzweireihe.first;
//         std::vector<Point>::const_iterator MaxdistIzwei = resultzweireihe.second;

//         int conessign = getSign(erstreihe[0].y - zweireihe[0].y);


//         if(m > n){
//             double dx = erstreihe[i - 1].x - erstreihe[i- 2].x;
//             double dy = erstreihe[i - 1].y - erstreihe[i - 2].y;
//             double nx = 0.4 * dy / sqrt(dx * dx + dy * dy);
//             double ny = -0.4 * dx / sqrt(dx * dx + dy * dy);
//             targetpoint = {erstreihe[i - 1].x - conessign * nx, erstreihe[i - 1].y - conessign * ny};
//         }else{
//             double dx = zweireihe[i - 1].x - zweireihe[i- 2].x;
//             double dy = zweireihe[i - 1].y - zweireihe[i - 2].y;
//             double nx = 0.4 * dy / sqrt(dx * dx + dy * dy);
//             double ny = -0.4 * dx / sqrt(dx * dx + dy * dy);
//             targetpoint = {zweireihe[i - 1].x - conessign * nx, zweireihe[i - 1].y - conessign * ny};

//         }
//     }

//      else if(i > 1){
//         int conessign = getSign(erstreihe[0].y - zweireihe[0].y);
//         startpoint = {(erstreihe[0].x + zweireihe[0].x) / 2, (erstreihe[0].y + zweireihe[0].y) / 2};
//         double dx = erstreihe[i - 1].x - erstreihe[i- 2].x;
//         double dy = erstreihe[i - 1].y - erstreihe[i - 2].y;
//         double nx = 0.4 * dy / sqrt(dx * dx + dy * dy);
//         double ny = -0.4 * dx / sqrt(dx * dx + dy * dy);
//         targetpoint = {erstreihe[i - 1].x - conessign * nx, erstreihe[i - 1].y - conessign * ny};

//      }else{
//         int conessign = getSign(erstreihe[0].y - zweireihe[0].y);
//         startpoint = {(erstreihe[0].x + zweireihe[0].x) / 2, (erstreihe[0].y + zweireihe[0].y) / 2};
//         double dx = zweireihe[i - 1].x - zweireihe[i- 2].x;
//         double dy = zweireihe[i - 1].y - zweireihe[i - 2].y;
//         double nx = 0.4 * dy / sqrt(dx * dx + dy * dy);
//         double ny = -0.4 * dx / sqrt(dx * dx + dy * dy);
//         targetpoint = {zweireihe[i - 1].x - conessign * nx, zweireihe[i - 1].y - conessign * ny};
//      }

// }


inline void vorstartpoint(const std::vector<Point>& erstreihe, const std::vector<Point>& zweireihe,  Point& zusatzpoint){
    unsigned int i = erstreihe.size();
    unsigned int j = zweireihe.size();
    if (i == 0 || j == 0 ) {
        std::cout << "erstreihe oder zweireihe ist leer---------------------------------- " << std::endl << std::endl;
        return;
    }

    int conessign = getSign(erstreihe[0].y - zweireihe[0].y);
    if (i >= j &&  i > 1){
        double dx = erstreihe[1].x - erstreihe[0].x;
        double dy = erstreihe[1].y - erstreihe[0].y;
        double nx = 0.4 * dy / sqrt(dx * dx + dy * dy);
        double ny = -0.4 * dx / sqrt(dx * dx + dy * dy);
        Point vorpoint = {erstreihe[0].x-dx, erstreihe[0].y-dy};
        zusatzpoint = {vorpoint.x - conessign * nx, vorpoint.y - conessign * ny};
    } else if (i <j && j >1){
        double dx = zweireihe[1].x - zweireihe[0].x;
        double dy = zweireihe[1].y - zweireihe[0].y;
        double nx = 0.4 * dy / sqrt(dx * dx + dy * dy);
        double ny = -0.4 * dx / sqrt(dx * dx + dy * dy);
        Point vorpoint = {zweireihe[0].x-dx, zweireihe[0].y-dy};
        zusatzpoint = {vorpoint.x - conessign * nx, vorpoint.y - conessign * ny};
    }
}

inline void endpoint(const std::vector<Point>& erstreihe, const std::vector<Point>& zweireihe,  Point& zusatzpoint){
    unsigned int i = erstreihe.size();
    unsigned int j = zweireihe.size();
    if (i == 0 || j == 0 ) {
        std::cout << "erstreihe oder zweireihe ist leer--------- " << std::endl << std::endl;
        return;
    }

     int conessign = getSign(erstreihe[0].y - zweireihe[0].y);
    if (i >= j &&  i > 10){
        double dx = erstreihe[i-1].x - erstreihe[i-2].x;
        double dy = erstreihe[i-1].y - erstreihe[i-2].y;
        double nx = 0.4 * dy / sqrt(dx * dx + dy * dy);
        double ny = -0.4 * dx / sqrt(dx * dx + dy * dy);
        Point endpoint = {erstreihe[i-1].x + dx, erstreihe[i-1].y + dy};
        zusatzpoint = {endpoint.x - conessign * nx, endpoint.y - conessign * ny};
    } else if (i <j && j >10){
        double dx = zweireihe[i-1].x - zweireihe[i-2].x;
        double dy = zweireihe[i-1].y - zweireihe[i-2].y;
        double nx = 0.4 * dy / sqrt(dx * dx + dy * dy);
        double ny = -0.4 * dx / sqrt(dx * dx + dy * dy);
        Point endpoint = {zweireihe[i-1].x + dx, zweireihe[i-1].y + dy};
        zusatzpoint = {endpoint.x - conessign * nx, endpoint.y - conessign * ny};
    }
}

inline void linkrightcones( std::vector<Point>& ConePoints, const Point& currentPosition,std::vector<Point>& linkcones,std::vector<Point>& rightcones){

    if(ConePoints.size() < 3){
        std::cout << "conesanzahl ist klein als ________ 3 " <<std::endl;
        return;
    }


    std::sort(ConePoints.begin(), ConePoints.end(), [&currentPosition](const Point& p1, const Point& p2) {
        return compareByDistance(p1, p2, currentPosition);
    });


    removeAllZeroPoints(ConePoints);

    std::vector<Point>::iterator it;

    if (ConePoints[0].x > 0.01) {
        rightcones.push_back(ConePoints[0]);
        for (it = ConePoints.begin() + 1; it != ConePoints.end(); ++it){
            if (it->x < 0){
                linkcones.push_back(*it);
                break;
            }
        }

    } else if (ConePoints[0].x < -0.01 ){
        linkcones.push_back(ConePoints[0]);
        for (it = ConePoints.begin() + 1; it != ConePoints.end(); ++it){
            if (it->x > 0){
                rightcones.push_back(*it);
                break;
            }
        }
    }


    if(linkcones.empty() && rightcones.empty() ){
        std::cout << "linkcones und rightcones sind null, können nicht link und right underscheiden " <<std::endl;
        return;
    }

    std::vector<Point> destinationPoints;

    for (const auto& point : ConePoints) {
        if (linkcones.empty()){
            if((point.x == rightcones[0].x && point.y == rightcones[0].y)||(point.x == 0 && point.y == 0)){
                    continue;
            }
        }else if (rightcones.empty()){
            if((point.x == linkcones[0].x && point.y == linkcones[0].y)||(point.x == 0 && point.y == 0)){
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
        std::cout << "destinationPoints ist klein________ _________ " <<std::endl<<std::endl;
        return;
        // throw std::runtime_error("destinationPoints ist klein________ _________ ");
        }

    if (linkcones.empty()&&rightcones.empty()){
        std::cout << "rightcones and linkcones must contain points! "<<std::endl <<std::endl;
        return;
        //    throw std::runtime_error("rightcones and linkcones must contain points!");
        }

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

    std::cout << linkcones.size()<< " tmd"<<std::endl;
    std::cout << rightcones.size()<< " tmd"<<std::endl;

    for (const Point& Position : rightcones) {
            std::cout << "rightcones (" << Position.x << ", " << Position.y << ")" << std::endl;
        }
    for (const Point& Position : linkcones) {
            std::cout << "linkcones (" << Position.x << ", " << Position.y << ")" << std::endl;
        }

}



inline void startandtraget(const std::vector<Point>& linkcones, const std::vector<Point>& rightcones,  Point& startpoint, Point& targetpoint) {

    unsigned int i = linkcones.size();
    unsigned int j = rightcones.size();

    if (i < 2 && j < 2) {
        std::cout << "cones zu wenig,  können nicht startpoint rechnen "<<std::endl<<std::endl;
        return;
        // throw std::runtime_error("cones zu wenig ");
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