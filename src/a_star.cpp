#include "a_star.hpp"
#include <chrono>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <Eigen/Geometry>
#include <glog/logging.h>
#include "types.hpp"
#include "internal/utils.hpp"


namespace kal_a_star {

// std::vector<Point> linkcones;
// std::vector<Point> rightcones(10);
// Point targetpoint;
// Point startpoint;



std::vector<Node*> A_star::aStarSearch(const Node& startNode, const Node& targetNode, const std::vector<Point>& linkcone, const std::vector<Point>& rightcone ){
    std::vector<Node*> openList;
    std::vector<Node*> closedList;
    openList.push_back(new Node(startNode.x, startNode.y));

    int n = linkcone.size();
    int m = rightcone.size();
    Eigen::VectorXd link_keoff;
    Eigen::VectorXd right_keoff;

    if (n > 2 && m > 2){

        link_keoff = utils::fitPolynomial(linkcone,std::min(n,m));
        right_keoff = utils::fitPolynomial(rightcone,std::min(n,m));

    }else if(m > 2){

        right_keoff = utils::fitPolynomial(rightcone,m);

        link_keoff = right_keoff;
        double mittelpointderivative = utils::calculateDerivative(right_keoff,(startNode.x + targetNode.x)/2);
        link_keoff.coeffRef(0) += 0.4 * sqrt(mittelpointderivative * mittelpointderivative + 1) / mittelpointderivative;

    }else{
        link_keoff = utils::fitPolynomial(linkcone,n);

        right_keoff = link_keoff;
        double mittelpointderivative = utils::calculateDerivative(link_keoff,(startNode.x + targetNode.x)/2);
        right_keoff.coeffRef(0) -= 0.4 * sqrt(mittelpointderivative * mittelpointderivative + 1) / mittelpointderivative;
    }


    while (!openList.empty()) {
        Node* currentNode = utils::getMinFCostNode(openList);

        // 如果当前节点为目标节点，则构建路径并返回
        if ( (std::pow( (currentNode->x - targetNode.x),2) + std::pow( (currentNode->y - targetNode.y),2) ) < std::pow(0.03,2) ) {
                return utils::buildPath(currentNode);
        }

        // 将当前节点从开放列表中移除，并添加到封闭列表中
        openList.erase(std::remove(openList.begin(), openList.end(), currentNode), openList.end());
        closedList.push_back(currentNode);


        // 遍历当前节点周围的邻居节点
        for (int dx = -1; dx <= 1; dx += 1) {
            for (int dy = -1; dy <= 1; dy += 1) {
                // if (dx==0 && dy == 0){
                 //    continue;
               //  }
                double neighborX = (currentNode->x) + (0.03 * dx);
               // cout << dx << "r"<< currNode.x <<"ss" <<endl;
                double neighborY = (currentNode->y) + (0.03 * dy);

                Point target = {neighborX, neighborY};

                // std::cout << neighborX << "q" <<endl;
                // 跳过无效的邻居节点 || utils::isWithinMap(neighborX, neighborY)

                if ( utils::isObstacle(neighborX, neighborY, rightcone) || utils::isObstacle(neighborX, neighborY, linkcone) ||
                utils::isNodeInList(closedList, neighborX, neighborY) || !utils::isWithinMap(neighborX, neighborY) ) {
                    // std::cout<< "x ";
                    continue;
                }
                // cout << "aaa";
                double neighborGCost = currentNode->gCost + std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));

                // 计算邻居节点到目标节点的估算代价
                double neighborHCost = COST_FREE * utils::calculateHCost(neighborX, neighborY, targetNode.x, targetNode.y);

                // 计算碰撞代价
                double sumDistances = utils::sumDistancesToTarget(rightcone, target) + utils::sumDistancesToTarget(linkcone, target);
                neighborHCost +=  COST_COLLISION * sumDistances;


                // 计算偏离中心线代价
                double diff = utils::polynomial_center((link_keoff+right_keoff)/2,neighborX);
                // double diff = polynomial(neighborX) + 3 ;
                double deviation = std::abs(neighborY - diff);
                if (deviation > distanztocenter) {
                    neighborHCost += ( deviation * deviation ) * COST_DEVIATION;
                }


                // 检查邻居节点是否已经在开放列表中
                bool isInOpenList = utils::isNodeInList(openList, neighborX, neighborY);
                if (!isInOpenList || neighborGCost < currentNode->gCost) {
                    Node* neighborNode;

                    if (isInOpenList) {
                        neighborNode = *std::find_if(openList.begin(), openList.end(),
                                                     [neighborX, neighborY](const Node* node) {
                                                         return (node->x == neighborX && node->y == neighborY);
                                                     });
                    } else {
                        neighborNode = new Node(neighborX, neighborY);
                    }

                    neighborNode->gCost = neighborGCost;
                    neighborNode->hCost = neighborHCost;
                    neighborNode->parent = currentNode;

                    if (!isInOpenList) {
                        openList.push_back(neighborNode);
                    }
                }
            }
        }
    }

    // 如果找不到路径，则返回空路径
     return std::vector<Node*>();
}




std::vector<Node*> A_star::aStar(const Pose currentPosition){
    std::vector<Point> linkcones;
    std::vector<Point> rightcones;
    Point currentPoint= {currentPosition.translation().x(),currentPosition.translation().y()};



    if(conesposition_.size() < 3){
        throw std::runtime_error("conesanzahl ist klein________ _________ _________");
        // std::cout << "conesanzahl ist klein________ _________ _________" << std::endl;

    }


    utils::linkrightcones(conesposition_, currentPoint,linkcones,rightcones);

    Point startpoint, targetpoint;
    utils::startandtraget(linkcones,rightcones,startpoint,targetpoint);
    Node start(startpoint.x, startpoint.y);
    Node traget(targetpoint.x, targetpoint.y);


    // linkcones = {{-0.4, 1.0}, {-0.4, 1.3}, {-0.4, 1.6}, {-0.4, 1.9}};
    // rightcones ={ };
    std::vector<Node*> path = aStarSearch(start,traget,linkcones,rightcones);
    return path;



}





void A_star::cones(const std::vector<Point>& cones){
    conesposition_ = cones;

}

// void A_star::cones(const std::vector<Point>& cones){
//     conesposition_ = cones;
//     utils::linkrightcones(conesposition_, currentposition, linkcones, rightcones);

//     for (const Point& Position : linkcones) {

//         std::cout << "Point (" << Position.x << ", " << Position.y << ")" << std::endl;
//     }
//     std::cout << "------------------" << std::endl;


//     utils::startandtraget( linkcones, rightcones,   startpoint, targetpoint);
//     std::cout << "startpoint (" << startpoint.x << ", " << startpoint.y << ")" << std::endl;
//     std::cout << "targetpoint (" << targetpoint.x << ", " << targetpoint.y << ")" << std::endl;







// }


// utils::linkrightcones(conesposition_, currentPosition,linkcones,rightcones);

}