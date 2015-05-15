#pragma once

#include <list>
#include <map>
#include <array>
#include <memory>
#include <boost/heap/fibonacci_heap.hpp>

#include "segment.h"

typedef std::shared_ptr<Segment> SegmentSharedPtr;

struct Edge {
    struct compare_node
    {
        bool operator()(const std::shared_ptr<Edge>& n1, const std::shared_ptr<Edge>& n2) const
        {
            return n1->cost > n2->cost;
        }
    };

    Edge(SegmentSharedPtr a,SegmentSharedPtr b)
        : cost(0.0), segment1(a), segment2(b) {}

    bool operator()(const std::shared_ptr<Edge>& n1, const std::shared_ptr<Edge>& n2) const
    {
        return n1->cost > n2->cost;
    }

    double cost;
    SegmentSharedPtr segment1;
    SegmentSharedPtr segment2;
    boost::heap::fibonacci_heap<std::shared_ptr<Edge>, boost::heap::compare<compare_node>>::handle_type handle;
};
typedef std::shared_ptr<Edge> EdgeSharedPtr;
typedef boost::heap::fibonacci_heap<EdgeSharedPtr, boost::heap::compare<Edge::compare_node>> Heap;
typedef Heap::handle_type HeapHandle;

class BaatzAlgorithm
{
public:
    typedef std::vector<SegmentSharedPtr> SegmentList;

    BaatzAlgorithm(double cw, double cmw, unsigned int s)
        : m_colorWeight(cw), m_compactnessWeight(cmw), m_scale(s)
    {}

    void loadPixelFromArray(const std::vector<unsigned char>& pixels,
                            size_t bandCount, size_t width, size_t height);
    void segmentation();

    void merge(EdgeSharedPtr);

    SegmentList getSegments() const;
private:
    typedef std::list<HeapHandle> EdgeList;
    typedef std::map<SegmentSharedPtr,EdgeList> SegmentMap;

    double colorComponentCost(EdgeSharedPtr);
    double compactnessComponentCost(EdgeSharedPtr);
    double computeCostFunction(EdgeSharedPtr);

    double m_colorWeight;
    double m_compactnessWeight;
    unsigned int m_scale;
    size_t m_bandCount;
    SegmentMap m_segments;
    Heap m_heap;
    size_t m_width;
    size_t m_height;
};
