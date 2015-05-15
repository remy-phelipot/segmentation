#include <cmath>
#include <set>
#include "baatzalgorithm.h"
#include <iostream>

double BaatzAlgorithm::computeCostFunction(EdgeSharedPtr e)
{
    double spectral_h;
    double spatial_h;
    double cost = 0;

    spectral_h = colorComponentCost(e);
    cost += m_colorWeight * spectral_h;

    if (cost < m_scale)
    {
        spatial_h = compactnessComponentCost(e);
        cost += (1 - m_compactnessWeight)*spatial_h;
    }
    return cost;
}


double BaatzAlgorithm::colorComponentCost(EdgeSharedPtr edge)
{
    auto r1 = edge->segment1;
    auto r2 = edge->segment2;
    unsigned int bands = r1->m_averageColor.size();
    std::vector<double> mean(bands), colorSum(bands), squarePixels(bands),
                        stddevNew(bands);
    double a_current = r1->m_area;
    double a_neighbor = r2->m_area;
    double a_sum = a_current + a_neighbor;

    for (size_t b = 0; b < bands; b++)
    {
        mean[b] = ((r1->m_averageColor[b] * a_current) + (r2->m_averageColor[b] * a_neighbor)) / a_sum;
        squarePixels[b] = (r1->m_averageColorSquare[b]) + (r2->m_averageColorSquare[b]);
        colorSum[b] = r1->m_colorSum[b] + r2->m_colorSum[b];
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
        double color_f = (a_current*r1->m_stdColor[b]) + (a_neighbor*r2->m_stdColor[b]);
        color_f = (a_sum*stddev) - color_f;
        color_h += color_f;
    }
    return color_h;
}

double BaatzAlgorithm::compactnessComponentCost(EdgeSharedPtr edge)
{
    auto r1 = edge->segment1;
    auto r2 = edge->segment2;
    double spatial_h, smooth_f, compact_f;
    std::vector<double> area(3), perimeter(3), b_box_len(3); /* 0-current segment; 1-neighbor segment; 2-merged (new) segment */

    /* area */
    area[0] = r1->m_area;
    area[1] = r2->m_area;
    area[2] = area[0] + area[1];

    /* perimeter */
    perimeter[0] = r1->m_perimeter;
    perimeter[1] = r2->m_perimeter;
    perimeter[2] = r1->m_perimeter + r2->m_perimeter - 2 * r2->m_connections;

    /* bounding box lenght */
    auto mbbox = r1->m_boundingBox.merge(r2->m_boundingBox);
    b_box_len[0] = (r1->m_boundingBox.getSize(0)) * 2 + (r1->m_boundingBox.getSize(1)) * 2;
    b_box_len[1] = (r2->m_boundingBox.getSize(0)) * 2 + (r2->m_boundingBox.getSize(1)) * 2;
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
    SegmentList segments;
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
        }
    }

    // Init the edges
    for (size_t y = 0; y < height - 1; y++)
        for (size_t x = 0; x < width - 1; x++)
        {
            const unsigned int base = y * width + x;
            auto rightEdge = std::make_shared<Edge>(segments[base],segments[base+1]);
            auto bottomEdge = std::make_shared<Edge>(segments[base],segments[base+width]);

            rightEdge->cost = computeCostFunction(rightEdge);
            bottomEdge->cost = computeCostFunction(bottomEdge);

            HeapHandle rightHandle = m_heap.push(rightEdge);
            (*rightHandle)->handle = rightHandle; // store handle in node
            HeapHandle bottomHandle = m_heap.push(bottomEdge);
            (*bottomHandle)->handle = bottomHandle; // store handle in node

            m_segments[segments[base]].push_back(rightHandle);
            m_segments[segments[base+1]].push_back(rightHandle);
            m_segments[segments[base]].push_back(bottomHandle);
            m_segments[segments[base+width]].push_back(bottomHandle);
        }

    m_width = width;
    m_height = height;
}

void BaatzAlgorithm::segmentation()
{
    size_t nb_it = 0;
    double threshold = 2*m_scale*m_scale;

    std::cout << "segmentation start " << m_heap.size() << " edges" << std::endl;
    while(!m_heap.empty() && m_heap.top()->cost <= threshold) {
        std::cout << m_heap.top()->cost << " / " << threshold << std::endl;
        auto edge = m_heap.top();

        merge(edge);

        m_heap.pop();

        nb_it++;
    }
    std::cout << "segmentation end " << m_heap.size() << " edges with "<< nb_it << " iterations" << std::endl;
}

void BaatzAlgorithm::merge(EdgeSharedPtr edge)
{
    // Remove the edge
    m_segments[edge->segment1].remove(edge->handle);
    m_segments[edge->segment2].remove(edge->handle);

    // Merge segment structure
    edge->segment1->merge(*edge->segment2);

    // Replace the deleted segment by the merged one in the neighborhood list
    for(HeapHandle it : m_segments[edge->segment2]) {
        EdgeSharedPtr neihborEdge = *it;
        if(neihborEdge->segment1 == edge->segment2)
            neihborEdge->segment1 = edge->segment1;
        if(neihborEdge->segment2 == edge->segment2)
            neihborEdge->segment2 = edge->segment1;
    }

    // Merge neighbors
    m_segments[edge->segment1].insert(m_segments[edge->segment1].end(),
            std::move(m_segments[edge->segment2].begin()),
            std::move(m_segments[edge->segment2].end()));

    // Update costs
    for(HeapHandle it : m_segments[edge->segment1]) {
        EdgeSharedPtr neihborEdge = *it;
        double oldCost = neihborEdge->cost ;
        neihborEdge->cost = computeCostFunction(neihborEdge);
        m_heap.update(it);
      /*  if(neihborEdge->cost < oldCost)
            m_heap.decrease(it);
        else if(neihborEdge->cost > oldCost)
            m_heap.increase(it);*/
    }

    m_segments.erase(edge->segment2);
}

BaatzAlgorithm::SegmentList BaatzAlgorithm::getSegments() const
{
    std::vector<SegmentSharedPtr> result;

    for(auto& edge : m_segments) {
        result.push_back(edge.first);
    }

    return result;
}
