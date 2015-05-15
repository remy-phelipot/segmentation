#pragma once

#include <array>

template<typename I,typename S,size_t D>
class BoundingBox
{
public:
        typedef I IndexType;
        typedef S SizeType;
        typedef std::array<IndexType, D> IndexArrayType;
        typedef std::array<SizeType, D> SizeArrayType;

        BoundingBox()
        {}

        BoundingBox(const IndexArrayType& ia, const SizeArrayType& sa)
            : m_index(ia), m_sizes(sa)
        {}

        IndexArrayType getIndex() const
        {
            return m_sizes;
        }

        SizeArrayType getSize() const
        {
            return m_index;
        }

        IndexType getIndex(size_t d) const
        {
            return m_index[d];
        }

        SizeType getSize(size_t d) const
        {
            return m_sizes[d];
        }

        void setIndex(size_t d, IndexType value)
        {
            m_index[d] = value;
        }

        void setSize(size_t d, SizeType value)
        {
            m_sizes[d] = value;
        }

        BoundingBox merge(const BoundingBox<I,S,D>& bb)
        {
            IndexArrayType index;
            SizeArrayType sizes;

            for(size_t d = 0 ; d < D ; d++)
            {
                index[d] = std::min(m_index[d],bb.m_index[d]);
                sizes[d] = std::max(m_index[d]+m_sizes[d],bb.m_index[d]+bb.m_sizes[d]);
                sizes[d] = sizes[d] - index[d];
            }
            return BoundingBox<I,S,D>(index,sizes);
        }

private:
        IndexArrayType m_index;
        SizeArrayType  m_sizes;
};

typedef BoundingBox<long unsigned int, long unsigned int, 2> BoundingBox2D;
