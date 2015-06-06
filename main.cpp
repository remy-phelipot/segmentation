#include <iostream>

#include <QImage>
#include "baatzalgorithm.h"

using namespace std;

int main(int argc, char **argv)
{
    QImage image("/home/remy/sf1.jpg");
    BaatzAlgorithm ba(0.9,0.9,100000);

    uchar* ptr = image.bits();
    std::vector<uchar> pixels;

    for(size_t y = 0 ; y < image.height() ; y++)
        for(size_t x = 0 ; x < image.width() ;  x++){
            auto pixel = image.pixel(x,y);
            pixels.push_back(qRed(pixel));
            pixels.push_back(qGreen(pixel));
            pixels.push_back(qBlue(pixel));
        }


    ba.loadPixelFromArray(pixels,3,image.width(),image.height());

    ba.segmentation();

    QImage result(image.width(),image.height(), QImage::Format_RGB16);

    for(auto& s : ba.getSegments()) {
        auto v = s->m_averageColor;
        auto color = qRgb(v[0],v[1],v[2]);
        for(auto pixel : s->m_pixels) {
            size_t x = pixel % image.width();
            size_t y = pixel / image.width();

            result.setPixel(x,y,color);
        }
    }

    result.save("out.png");

    return 0;
}
