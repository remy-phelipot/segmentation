#pragma once

#include <vector>

#include "boundingbox.h"

class Segment
{
public:
    typedef std::vector<double> DoubleList;
    typedef std::vector<long unsigned int> IndexList;

    Segment(size_t bandCount, long unsigned int id)
    : m_area(0.0), m_perimeter(0.0),
      m_averageColor(bandCount), m_averageColorSquare(bandCount),
      m_colorSum(bandCount), m_stdColor(bandCount),
      m_id(id), m_connections(1)
    {}

    void merge(const Segment& other);

    double m_area;
    double m_perimeter;
    DoubleList m_averageColor;
    DoubleList m_averageColorSquare;
    DoubleList m_colorSum;
    DoubleList m_stdColor;

    unsigned int m_connections;
    long unsigned int m_id;
    BoundingBox2D m_boundingBox;
    IndexList m_pixels;
private:
    Segment(){}
};
