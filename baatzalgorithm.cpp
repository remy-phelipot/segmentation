#include <cmath>
#include <set>
#include "baatzalgorithm.h"
#include <iostream>
#include <algorithm>

double BaatzAlgorithm::computeCostFunction(const SegmentSharedPtr& s1,const SegmentSharedPtr& s2)
{
    double spectral_h;
    double spatial_h;
    double cost = 0;

    spectral_h = colorComponentCost(s1, s2);
    cost += m_colorWeight * spectral_h;

    if (cost < m_scale)
    {
        spatial_h = compactnessComponentCost(s1, s2);
        cost += (1 - m_compactnessWeight)*spatial_h;
    }
    return cost;
}


double BaatzAlgorithm::colorComponentCost(const SegmentSharedPtr& s1,const SegmentSharedPtr& s2)
{
    unsigned int bands = s1->m_averageColor.size();
    std::vector<double> mean(bands), colorSum(bands), squarePixels(bands),
                        stddevNew(bands);
    double a_current = s1->m_area;
    double a_neighbor = s2->m_area;
    double a_sum = a_current + a_neighbor;

    for (size_t b = 0; b < bands; b++)
    {
        mean[b] = ((s1->m_averageColor[b] * a_current) + (s2->m_averageColor[b] * a_neighbor)) / a_sum;
        squarePixels[b] = (s1->m_averageColorSquare[b]) + (s2->m_averageColorSquare[b]);
        colorSum[b] = s1->m_colorSum[b] + s2->m_colorSum[b];
        stddevNew[b] = 0.0;
    }

    for (size_t b = 0; b < bands; b++)
    {
        stddevNew[b] = squarePixels[b] - 2 * mean[b] * colorSum[b] + a_sum*mean[b] * mean[b];
    }

    /* calculates color factor per band and total */
    double color_h = 0.0;
    for (size_t b = 0; b < bands; b++)
    {
        const double stddev = std::sqrt(stddevNew[b] / a_sum);
        double color_f = (a_current*s1->m_stdColor[b]) + (a_neighbor*s2->m_stdColor[b]);
        color_f = (a_sum*stddev) - color_f;
        color_h += color_f;
    }
    return color_h;
}

double BaatzAlgorithm::compactnessComponentCost(const SegmentSharedPtr& s1,const SegmentSharedPtr& s2)
{
    double spatial_h, smooth_f, compact_f;
    std::vector<double> area(3), perimeter(3), b_box_len(3); /* 0-current segment; 1-neighbor segment; 2-merged (new) segment */

    /* area */
    area[0] = s1->m_area;
    area[1] = s2->m_area;
    area[2] = area[0] + area[1];

    /* perimeter */
    perimeter[0] = s1->m_perimeter;
    perimeter[1] = s2->m_perimeter;
    perimeter[2] = s1->m_perimeter + s2->m_perimeter - 2 * s2->m_connections;

    /* bounding box lenght */
    auto mbbox = s1->m_boundingBox.merge(s2->m_boundingBox);
    b_box_len[0] = (s1->m_boundingBox.getSize(0)) * 2 + (s1->m_boundingBox.getSize(1)) * 2;
    b_box_len[1] = (s2->m_boundingBox.getSize(0)) * 2 + (s2->m_boundingBox.getSize(1)) * 2;
    b_box_len[2] = (mbbox.getSize(0)) * 2 + (mbbox.getSize(1)) * 2;

    /* smoothness factor */
    smooth_f = (area[2] * perimeter[2] / b_box_len[2] -
        (area[1] * perimeter[1] / b_box_len[1] + area[0] * perimeter[0] / b_box_len[0]));

    /* compactness factor */
    compact_f = (area[2] * perimeter[2] / std::sqrt(area[2]) -
        (area[1] * perimeter[1] / std::sqrt(area[1]) + area[0] * perimeter[0] / std::sqrt(area[0])));

    /* spatial heterogeneity */
    spatial_h = m_compactnessWeight*compact_f + (1 - m_compactnessWeight)*smooth_f;

    return spatial_h;
}

void BaatzAlgorithm::loadPixelFromArray(const std::vector<unsigned char>& pixels,
                                        size_t bandCount, size_t width, size_t height)
{
    std::vector<SegmentSharedPtr> segments;
    for (size_t y = 0; y < height; y++)
    {
        for (size_t x = 0; x < width; x++)
        {
            const size_t array_base = (y * width + x) * bandCount;
            BoundingBox2D bb(BoundingBox2D::IndexArrayType{ { x, y } },
                             BoundingBox2D::SizeArrayType{ { 1, 1 } });
            auto segment = std::make_shared<Segment>(bandCount, y * width + x);

            segment->m_area = 1;
            segment->m_perimeter = 4;
            segment->m_boundingBox = bb;
            segment->m_pixels.push_back(y * width + x);

            for (size_t b = 0; b < bandCount; b++)
            {
                const unsigned int band_id = array_base + b;
                segment->m_averageColor[b] = pixels[band_id];
                segment->m_averageColorSquare[b] = std::pow(pixels[band_id], 2);
                segment->m_colorSum[b] = pixels[band_id];
                segment->m_stdColor[b] = 0.0;
            }

            segments.push_back(segment);
            m_segments.push_back(segment);
        }
    }

    // Init the edges
    for (size_t y = 0; y < height - 1; y++)
        for (size_t x = 0; x < width - 1; x++)
        {
            const unsigned int base = y * width + x;

            auto baseRegion = segments[base];
            auto rightRegion = segments[base+1];
            auto bottomRegion = segments[base+width];

            auto rightCost = std::make_shared<double>(computeCostFunction(baseRegion,rightRegion));
            auto bottomCost = std::make_shared<double>(computeCostFunction(baseRegion,bottomRegion));

            baseRegion->m_neighbors.push_back(Neighbor(rightRegion,rightCost));
            baseRegion->m_neighbors.push_back(Neighbor(bottomRegion,bottomCost));
            rightRegion->m_neighbors.push_back(Neighbor(baseRegion,rightCost));
            bottomRegion->m_neighbors.push_back(Neighbor(baseRegion,bottomCost));
        }

    for (size_t x = 0; x < width - 1; x++)
    {
        const unsigned int base = (height-1) * width + x;

        auto baseRegion = segments[base];
        auto rightRegion = segments[base+1];

        auto rightCost = std::make_shared<double>(computeCostFunction(baseRegion,rightRegion));

        baseRegion->m_neighbors.push_back(Neighbor(rightRegion,rightCost));
        rightRegion->m_neighbors.push_back(Neighbor(baseRegion,rightCost));
    }

    for (size_t y = 0; y < height - 1; y++)
    {
        const unsigned int base = y * width + (width-1);

        auto baseRegion = segments[base];
        auto bottomRegion = segments[base+width];

        auto bottomCost = std::make_shared<double>(computeCostFunction(baseRegion,bottomRegion));

        baseRegion->m_neighbors.push_back(Neighbor(bottomRegion,bottomCost));
        bottomRegion->m_neighbors.push_back(Neighbor(baseRegion,bottomCost));
    }

    for(SegmentSharedPtr ptr : m_segments)
        ptr->m_neighbors.sort();

    m_width = width;
    m_height = height;
}

void BaatzAlgorithm::segmentation()
{
    size_t nb_it = 0;
    size_t previousSegmentCount = m_segments.size()+1;
    double threshold = m_scale;

    std::cout << "segmentation start " << m_segments.size() << " edges" << std::endl;
    while(previousSegmentCount > m_segments.size() && nb_it < 500 && m_segments.size() > 1) {
        std::cout << nb_it << " iterations with " << m_segments.size() << std::endl;
        previousSegmentCount = m_segments.size();
        for(SegmentSharedPtr r : m_segments) {
            if(!r->toDelete && !r->merged){
                Neighbor minE = r->m_neighbors.front();

                if(!minE.m_NeighborRegion->toDelete && !minE.m_NeighborRegion->merged && *(minE.m_cost) < threshold) {
                    Neighbor minNE = minE.m_NeighborRegion->m_neighbors.front();
                    if(minNE.m_NeighborRegion == r)
                        merge(r,minE.m_NeighborRegion);
                }
            }
        }

        m_segments.remove_if([](const SegmentSharedPtr ptr){return ptr->toDelete;});

        for(SegmentSharedPtr r : m_segments)
            r->merged = false;

        nb_it++;
    }
    std::cout << "segmentation end " << m_segments.size() << " edges with "<< nb_it << " iterations" << std::endl;
}

void BaatzAlgorithm::merge(SegmentSharedPtr s1, SegmentSharedPtr s2)
{
    // Merge segment structure
    s1->merge(*s2);

    // Update neighbor
    s1->m_neighbors.pop_front();
    s2->m_neighbors.pop_front();

    for(Neighbor n : s2->m_neighbors) {
        n.m_NeighborRegion->m_neighbors.remove_if([&](const Neighbor& o){return o.m_NeighborRegion == s1 || o.m_NeighborRegion == s2;});
        s1->m_neighbors.remove_if([&](const Neighbor& o){return o.m_NeighborRegion == n.m_NeighborRegion;});
        n.m_NeighborRegion->m_neighbors.push_back(Neighbor(s1,n.m_cost));
    }

    s1->m_neighbors.splice(s1->m_neighbors.end(),s2->m_neighbors);

    // Update costs
    for(Neighbor& n : s1->m_neighbors) {
        *n.m_cost = computeCostFunction(s1,n.m_NeighborRegion);
        n.m_NeighborRegion->m_neighbors.sort();
    }

    s1->m_neighbors.sort();

    s1->merged = true;
    s2->toDelete = true;
}

BaatzAlgorithm::SegmentList BaatzAlgorithm::getSegments() const
{
    return m_segments;
}
