#ifndef ENTITY2D_HPP
#define ENTITY2D_HPP

#include "body2D.hpp"
#include "box2D.hpp"
#include "force2D.hpp"
#include "interaction2D.hpp"
#include "polygon2D.hpp"
#include "vec_ptr.hpp"
#include <unordered_set>

namespace physics
{
    class entity2D : public body2D
    {
    public:
        using body2D::body2D;

        void retrieve(const utils::vec_ptr &buffer);
        void retrieve();

        bool contains(const force2D &force) const;
        bool contains(const interaction2D &inter) const;

        const geo::box2D &bounding_box() const;
        const geo::polygon2D &shape() const;

    private:
        geo::box2D m_bbox;
        geo::polygon2D m_shape;
        utils::vec_ptr m_buffer;
        vec2 accel;

        std::unordered_set<const force2D *> m_forces;
        std::unordered_set<const interaction2D *> m_inters;

        void add(const force2D &force);
        void add(const interaction2D &inter);

        void remove(const force2D &force);
        void remove(const interaction2D &inter);

        friend class force2D;
        friend class interaction2D;
    };
}

#endif