#ifndef INTERACTION2D_HPP
#define INTERACTION2D_HPP

#include "entity2D.hpp"

namespace phys
{
    class interaction2D
    {
    public:
        interaction2D() = default;
        virtual ~interaction2D() = default;
        virtual std::pair<alg::vec2, float> force(const entity2D &e1, const entity2D &e2) const = 0;
        virtual float potential_energy(const entity2D &e1, const entity2D &e2) const { return 0.f; }

        float potential(const entity2D &e, const alg::vec2 &pos) const;
        float energy(const entity2D &e1, const entity2D &e2) const;

        bool enabled() const;
        void enabled(bool enabled);

    private:
        bool m_enabled = true;
        mutable phys::entity2D m_unit;
    };
}
#endif