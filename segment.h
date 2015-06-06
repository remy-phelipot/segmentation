#pragma once

#include <vector>
#include <memory>
#include "boundingbox.h"
#include <list>
#include <queue>

class Segment;

typedef std::shared_ptr<double> CostSharedPtr;
typedef std::shared_ptr<Segment> SegmentSharedPtr;
typedef std::shared_ptr<Segment> SegmentWeakPtr;

struct Neighbor {
    Neighbor(const SegmentSharedPtr& a, const CostSharedPtr& c)
        : m_cost(c), m_NeighborRegion(a) {}

    CostSharedPtr m_cost;
    SegmentWeakPtr m_NeighborRegion;
};


inline bool operator<(const Neighbor& a, const Neighbor& b)
{
    return *a.m_cost < *b.m_cost;
}

inline bool operator>(const Neighbor& a, const Neighbor& b)
{
    return *a.m_cost > *b.m_cost;
}

class Segment
{
public:
    typedef std::vector<double> DoubleList;
    typedef std::vector<unsigned int> IndexList;
    typedef std::list<Neighbor> NeighborList;

    Segment(size_t bandCount, long unsigned int id)
    : m_area(0.0), m_perimeter(0.0),
      m_averageColor(bandCount), m_averageColorSquare(bandCount),
      m_colorSum(bandCount), m_stdColor(bandCount),
      m_id(id), m_connections(1), toDelete(false), merged(false)
    {}

    void merge(const Segment& other);

    double m_area;
    double m_perimeter;
    DoubleList m_averageColor;
    DoubleList m_averageColorSquare;
    DoubleList m_colorSum;
    DoubleList m_stdColor;
    NeighborList m_neighbors;

    unsigned int m_connections;
    unsigned int m_id;
    BoundingBox2D m_boundingBox;
    IndexList m_pixels;

    bool toDelete;
    bool merged;
private:
    Segment(){}
};
