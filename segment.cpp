#include "segment.h"
#include <cmath>

void Segment::merge(const Segment& seg)
{
    size_t bands = m_averageColor.size();
    std::vector<double> mean(bands), colorSum(bands), squarePixels(bands);

    double a_sum = m_area + seg.m_area;

    for (unsigned int b = 0; b < bands; b++)
    {
        mean[b] = ((m_averageColor[b]*m_area)+(seg.m_averageColor[b]*seg.m_area))/a_sum;
        squarePixels[b] = m_averageColorSquare[b] + seg.m_averageColorSquare[b];
        colorSum[b] = m_colorSum[b] + seg.m_colorSum[b];
    }

    for(unsigned int b = 0; b < bands; b++)
    {
        m_averageColor[b] = mean[b];
        m_averageColorSquare[b] = squarePixels[b];
        m_colorSum[b] = colorSum[b];
        m_stdColor[b] = std::sqrt((squarePixels[b] - 2*mean[b]*colorSum[b] + a_sum*mean[b]*mean[b])/a_sum);
    }

    // Update spatial attributes
    m_area += seg.m_area;
    m_perimeter += seg.m_perimeter - 2*seg.m_connections;
    m_connections += seg.m_connections;

    m_pixels.insert(m_pixels.end(),std::move(seg.m_pixels.begin()),std::move(seg.m_pixels.end()));
}
