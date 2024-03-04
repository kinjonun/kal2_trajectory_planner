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

    const int COST_FREE = 1;            // free path cost
    const double COST_COLLISION = 2;    // collision cost
    const double COST_DEVIATION = 200;  // Cost of deviation from centerline



 class A_star
    {
    private:
        std::vector<Point> conesposition_;
        std::vector<Point> linkcones_;
        std::vector<Point> rightcones_;
        std::vector<Point> erstreihe_;
        std::vector<Point> zweireihe_;

    public:

        void cones(const std::vector<Point>& cones);
        void linkandrightcones(const std::vector<Point>& linkcones, const std::vector<Point>& rightcones);
        void erstundzweireihe(const std::vector<Point>& erstreihe, const std::vector<Point>& zweireihe);


        std::vector<Node*> aStarSearch(const Node& startNode, const Node& targetNode, const std::vector<Point>& link, const std::vector<Point>& right );
        std::vector<Node*> aStar(const Pose currentPosition);
        const std::vector<Point>& getConesPosition() const { return conesposition_; }
        const int getConeszahl() const { return erstreihe_.size()+zweireihe_.size(); }

        std::vector<Position> Polynom(const Point& startpunkt, const Point& zielpunkt, const std::vector<Point>& linkcone, const std::vector<Point>& rightcone );
        std::vector<Position> fitPolynom(const Pose currentPosition);

        void fitmittel(const Pose currentPosition,  std::vector<Position>& interpolar);

        void fitcones(std::vector<Position>& interpolar);

    };



}
