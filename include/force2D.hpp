#ifndef FORCE2D_HPP
#define FORCE2D_HPP

#include "entity2D.hpp"

namespace phys
{
    class force2D
    {
    public:
        force2D() = default;
        virtual ~force2D() = default;
        virtual std::pair<alg::vec2, float> force(const entity2D &e) const = 0;

        virtual float potential_energy(const entity2D &e) const { return 0.f; }
        float energy(const entity2D &e) const;

        bool enabled() const;
        void enabled(bool enabled);

    private:
        bool m_enabled = true;
    };
}

#endif