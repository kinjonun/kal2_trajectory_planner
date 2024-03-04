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
    Point startpunkt = {startNode.x, startNode.y};
    Point zielpunkt = {targetNode.x, targetNode.y};

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

        // If the current node is the target node, construct the path and return the
        if ( (std::pow( (currentNode->x - targetNode.x),2) + std::pow( (currentNode->y - targetNode.y),2) ) < std::pow(0.15,2) ) {
                return utils::buildPath(currentNode);
        }

        // Remove the current node from the open list and add it to the closed list
        openList.erase(std::remove(openList.begin(), openList.end(), currentNode), openList.end());
        closedList.push_back(currentNode);


        // Iterate through the neighboring nodes around the current node
        for (int dx = -1; dx <= 1; dx += 1) {
            for (int dy = -1; dy <= 1; dy += 1) {
                // if (dx==0 && dy == 0){
                 //    continue;
               //  }
                double neighborX = (currentNode->x) + (0.15 * dx);
               // cout << dx << "r"<< currNode.x <<"ss" <<endl;
                double neighborY = (currentNode->y) + (0.15 * dy);

                Point target = {neighborX, neighborY};

                // std::cout << neighborX << "q" <<endl;
                // Skip invalid neighbor nodes || utils::isWithinMap(neighborX, neighborY)

                if ( utils::isObstacle(neighborX, neighborY, rightcone) || utils::isObstacle(neighborX, neighborY, linkcone) ||
                utils::isNodeInList(closedList, neighborX, neighborY) || !utils::isWithinMap(neighborX, neighborY,startpunkt, zielpunkt) ) {
                    // std::cout<< "x ";
                    continue;
                }
                // cout << "aaa";
                double neighborGCost = currentNode->gCost + std::sqrt(std::pow(0.15*dx, 2) + std::pow(0.15*dy, 2));

                // Calculate the estimated cost of neighboring nodes to the target node
                double neighborHCost = COST_FREE * utils::calculateHCost(neighborX, neighborY, targetNode.x, targetNode.y);

                // Calculating the cost of a collision
                double sumDistances = utils::sumDistancesToTarget(rightcone, target) + utils::sumDistancesToTarget(linkcone, target);
                neighborHCost +=  COST_COLLISION * sumDistances;


                // Calculation of the cost of deviation from the centerline
                double diff = utils::polynomial_center((link_keoff+right_keoff)/2,neighborX);
                // double diff = polynomial(neighborX) + 3 ;
                double deviation = std::abs(neighborY - diff);
                if (deviation > distanztocenter) {
                    neighborHCost += ( deviation * deviation ) * COST_DEVIATION;
                }


                // Check if neighbor nodes are already in the open list
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

    // If the path is not found, the empty path is returned
     return std::vector<Node*>();
}

std::vector<Position> A_star::Polynom(const Point& startpunkt, const Point& zielpunkt, const std::vector<Point>& linkcone, const std::vector<Point>& rightcone ){

    int n = linkcone.size();
    int m = rightcone.size();
    Eigen::VectorXd link_keoff;
    Eigen::VectorXd right_keoff;

    if (n > 2 && m > 2){

        link_keoff = utils::fitPolynomial(linkcone,std::min(n,m)-1);
        right_keoff = utils::fitPolynomial(rightcone,std::min(n,m)-1);

    }else if(m > 2){

        right_keoff = utils::fitPolynomial(rightcone,m-1);

        link_keoff = right_keoff;
        double mittelpointderivative = utils::calculateDerivative(right_keoff,(startNode.x + targetNode.x)/2);
        link_keoff.coeffRef(0) += std::abs(0.9 * sqrt(mittelpointderivative * mittelpointderivative + 1) / mittelpointderivative);

    }else{
        link_keoff = utils::fitPolynomial(linkcone,n-1);

        right_keoff = link_keoff;
        double mittelpointderivative = utils::calculateDerivative(link_keoff,(startNode.x + targetNode.x)/2);
        right_keoff.coeffRef(0) -= std::abs(0.9 * sqrt(mittelpointderivative * mittelpointderivative + 1) / mittelpointderivative);

    }

    for (int i=0; i < link_keoff.size(); i++){
        std::cout << "link_keoff: " << link_keoff[i] << std::endl;
    }
     for (int i=0; i < right_keoff.size(); i++){
        std::cout << "right_keoff: " << right_keoff[i] << std::endl;
    }


   std::vector<Position> interpolar;
    if (startpunkt.x <= zielpunkt.x){
        for (double i= startpunkt.x; i < zielpunkt.x; i= i + 0.03 ){
            double poly_y = utils::polynomial_center((link_keoff+right_keoff)/2,i);
            interpolar.push_back(Position(i, poly_y));
        }
    } else {
        for (double i= startpunkt.x; i > zielpunkt.x; i= i - 0.03 ){
            double poly_y = utils::polynomial_center((link_keoff+right_keoff)/2,i);
            interpolar.push_back(Position(i, poly_y));
        }
    }
    return interpolar;



}


std::vector<Position> A_star::fitPolynom(const Pose currentPosition){

    Point currentPoint= {currentPosition.translation().x(),currentPosition.translation().y()};
    // if(conesposition_.size() < 4){
    //     std::cout << "conesanzahl ist klein________ _________ _________  "<<std::endl;

    //     // throw std::runtime_error("conesanzahl ist klein________ _________ _________");
    //     // std::cout << "conesanzahl ist klein________ _________ _________" << std::endl;

    // }
    Point nullpunkt = {0,0};
    std::vector<Point> linkcones1;
    std::vector<Point> rightcones1;
    utils::linkrightcones(conesposition_, nullpunkt,linkcones1,rightcones1);
    Point startpoint1, targetpoint1;
    utils::startandtraget(linkcones1,rightcones1,startpoint1,targetpoint1);
    if((std::pow(targetpoint1.x, 2) + std::pow(targetpoint1.y, 2)) < 0.01){
        std::cout <<std::endl << "traget ist zu nah          "<<std::endl;
        // return std::vector<Node*>();
        return std::vector<Position>();

    }


    std::cout << "linkcones_.size: " << linkcones_.size()<< " tmd"<<std::endl;
    std::cout  << "rightcones_.size: " << rightcones_.size()<< " tmd"<<std::endl;

    Point startpoint, targetpoint;
    // startpoint = {0,1.0};
    // targetpoint = {0,2.0};

    utils::startandtraget(linkcones_,rightcones_,startpoint,targetpoint);

    std::cout << "startpoint: " << startpoint.x<< " "<< startpoint.y <<std::endl;
    std::cout  << "targetpoint: " << targetpoint.x << " "<< targetpoint.y <<std::endl;

    Node start(currentPoint.x, currentPoint.y);
    Node traget(targetpoint.x, targetpoint.y);


    // linkcones = {{-0.4, 1.0}, {-0.4, 1.3}, {-0.4, 1.6}, {-0.4, 1.9}};
    // rightcones ={ };
    // std::vector<Node*> path = aStarSearch(start,traget,linkcones_,rightcones_);
    std::vector<Position> eigenVectors = Polynom(currentPoint, targetpoint, linkcones_, rightcones_ );


    return eigenVectors;

}

std::vector<Node*> A_star::aStar(const Pose currentPosition){

    Point currentPoint= {currentPosition.translation().x(),currentPosition.translation().y()};



    // if(conesposition_.size() < 4){
    //     std::cout << "conesanzahl ist klein________ _________ _________  "<<std::endl;

    //     // throw std::runtime_error("conesanzahl ist klein________ _________ _________");
    //     // std::cout << "conesanzahl ist klein________ _________ _________" << std::endl;

    // }
    Point nullpunkt = {0,0};
    std::vector<Point> linkcones1;
    std::vector<Point> rightcones1;
    utils::linkrightcones(conesposition_, nullpunkt,linkcones1,rightcones1);
    Point startpoint1, targetpoint1;
    utils::startandtraget(linkcones1,rightcones1,startpoint1,targetpoint1);
    if((std::pow(targetpoint1.x, 2) + std::pow(targetpoint1.y, 2)) < 0.01){
        std::cout <<std::endl << "traget ist zu nah          "<<std::endl;
        return std::vector<Node*>();
    }


    std::cout << "linkcones_.size: " << linkcones_.size()<< " tmd"<<std::endl;
    std::cout  << "rightcones_.size: " << rightcones_.size()<< " tmd"<<std::endl;

    Point startpoint, targetpoint;
    // startpoint = {0,1.0};
    // targetpoint = {0,2.0};

    utils::startandtraget(linkcones_,rightcones_,startpoint,targetpoint);

    std::cout << "startpoint: " << startpoint.x<< " "<< startpoint.y <<std::endl;
    std::cout  << "targetpoint: " << targetpoint.x << " "<< targetpoint.y <<std::endl;

    Node start(currentPoint.x, currentPoint.y);
    Node traget(targetpoint.x, targetpoint.y);




    // linkcones = {{-0.4, 1.0}, {-0.4, 1.3}, {-0.4, 1.6}, {-0.4, 1.9}};
    // rightcones ={ };
    std::vector<Node*> path = aStarSearch(start,traget,linkcones_,rightcones_);

    return path;



}


void A_star::cones(const std::vector<Point>& cones){
    conesposition_ = cones;

}

void A_star::linkandrightcones(const std::vector<Point>& linkcones, const std::vector<Point>& rightcones){
    linkcones_ = linkcones;
    rightcones_ = rightcones;
}

void A_star::erstundzweireihe(const std::vector<Point>& erstreihe, const std::vector<Point>& zweireihe){
    erstreihe_ = erstreihe;
    zweireihe_ = zweireihe;
}

void A_star::fitcones(std::vector<Position>& interpolar){
    Point startpoint, targetpoint;
 //   utils::startandzielpoint( erstreihe_, zweireihe_,  startpoint, targetpoint);

    int n = erstreihe_.size();
    int m = zweireihe_.size();
    Eigen::VectorXd keoff;



    if(n < 1 || m < 1 ){
        return;
        std::cout << "hzuuhkjhdk--------------------------------------------" << std::endl;
    }

    // int conessign = utils::getSign(erstreihe_[0].y - zweireihe_[0].y);
    double distance = utils::distanceToPoint(erstreihe_[0], zweireihe_[0]);


    Point bezugpoint = {(erstreihe_[0].x + zweireihe_[0].x) / 2, (erstreihe_[0].y + zweireihe_[0].y) / 2};

    double distancetoerstreihe = utils::maxDistanceToPoints(bezugpoint, erstreihe_);
    double distancetozweireihe = utils::maxDistanceToPoints(bezugpoint, zweireihe_);

    Eigen::VectorXd koeff;

    if (distancetoerstreihe >= distancetozweireihe){
        koeff = utils::fitPolynomial(erstreihe_, std::min(3, n-1));

        double koeff_null = bezugpoint.y - (utils::polynomial_center(koeff,bezugpoint.x)) + koeff.coeffRef(0);
        koeff.coeffRef(0) = koeff_null;

        for (int i =0; i< koeff.size(); ++i){
            std::cout << "koeff:  " << i << "J= " << koeff[i] << std::endl;
        }


        if (erstreihe_[0].x < erstreihe_[n-1].x ){
            for (double i = erstreihe_[0].x; i < erstreihe_[n-1].x;  i= i + 0.02){
                double poly_y = utils::polynomial_center(koeff,i);
                interpolar.push_back(Position(i, poly_y));
            }
        }else{
            for (double i = erstreihe_[0].x; i > erstreihe_[n-1].x;  i= i - 0.02){
                double poly_y = utils::polynomial_center(koeff,i);
                interpolar.push_back(Position(i, poly_y));
            }

        }
    }
     else{

        koeff = utils::fitPolynomial(zweireihe_, std::min(3, m-1));
        for (int i =0; i< koeff.size(); ++i){
            std::cout << "koeff:  " << i << "J= " << koeff[i] << std::endl;
        }

        double ss = utils::polynomial_center(koeff,bezugpoint.x);
        double coeff_erst = koeff.coeffRef(0);
        double koeff_null = bezugpoint.y - ss + coeff_erst;


        koeff.coeffRef(0) = koeff_null;
        if (zweireihe_[0].x < zweireihe_[m-1].x ){
            for (double i = zweireihe_[0].x; i < zweireihe_[m-1].x;  i= i + 0.02){
                double poly_y = utils::polynomial_center(koeff,i);
                interpolar.push_back(Position(i, poly_y));
            }
        }else{
                for (double i = zweireihe_[0].x; i > zweireihe_[m-1].x;  i= i - 0.02){
                double poly_y = utils::polynomial_center(koeff,i);
                interpolar.push_back(Position(i, poly_y));
                }
            }
        }


    // std::cout << "interpolar.size: -------------------------------------------" << interpolar.size();


}





void A_star::fitmittel(const Pose currentPosition, std::vector<Position>& interpolar){
    Point currentPoint= {currentPosition.translation().x(),currentPosition.translation().y()};

    static std::vector<Point> mittelpoints;

    Point startpoint, targetpoint;
    // utils::startandzielpoint( erstreihe_, zweireihe_,  startpoint, targetpoint);


    utils::calculateMittelpoint(mittelpoints, erstreihe_, zweireihe_);

    mittelpoints.push_back(targetpoint);
    int minttelpunktanzahl = mittelpoints.size() ;

    Eigen::VectorXd koeff;
    int i = getConeszahl();
    if(i < 4){
        koeff = utils::fitPolynomial(mittelpoints, 1);
    }else{
        koeff = utils::fitPolynomial(mittelpoints, std::min( 4 , minttelpunktanzahl ));
    }


    if (minttelpunktanzahl < 1){
        return;
    }

    if (mittelpoints[0].x < mittelpoints[mittelpoints.size()-1].x){
        for (double i= mittelpoints[0].x; i < mittelpoints[mittelpoints.size()-1].x; i= i + 0.02 ){
            double poly_y = utils::polynomial_center(koeff,i);
            interpolar.push_back(Position(i, poly_y));
        }
    }else{
        for (double i= mittelpoints[0].x; i > mittelpoints[mittelpoints.size()-1].x; i= i - 0.02 ){
            double poly_y = utils::polynomial_center(koeff,i);
            interpolar.push_back(Position(i, poly_y));
        }
    }


    std::cout << "interpolar.size: " << interpolar.size();
    //utils::print_position(interpolar);







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