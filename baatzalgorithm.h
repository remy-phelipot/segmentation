#pragma once

#include <list>
#include <map>
#include <array>
#include <memory>

#include "segment.h"

class BaatzAlgorithm
{
public:
    typedef std::list<SegmentSharedPtr> SegmentList;

    BaatzAlgorithm(double cw, double cmw, unsigned int s)
        : m_colorWeight(cw), m_compactnessWeight(cmw), m_scale(s)
    {}

    void loadPixelFromArray(const std::vector<unsigned char>& pixels,
                            size_t bandCount, size_t width, size_t height);
    void segmentation();

    void merge(SegmentSharedPtr,SegmentSharedPtr);

    SegmentList getSegments() const;
private:
    double colorComponentCost(const SegmentSharedPtr&,const SegmentSharedPtr&);
    double compactnessComponentCost(const SegmentSharedPtr&,const SegmentSharedPtr&);
    double computeCostFunction(const SegmentSharedPtr&,const SegmentSharedPtr&);

    double m_colorWeight;
    double m_compactnessWeight;
    unsigned int m_scale;
    size_t m_bandCount;
    size_t m_width;
    size_t m_height;
    SegmentList m_segments;
};
