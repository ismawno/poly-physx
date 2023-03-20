#include "interaction2D.hpp"

namespace phys
{
    float interaction2D::potential(const entity2D &e, const alg::vec2 &pos) const
    {
        m_unit.pos(pos);
        return potential_energy(m_unit, e);
    }
    float interaction2D::energy(const entity2D &e1, const entity2D &e2) const
    {
        return e1.kinetic_energy() + e2.kinetic_energy() + potential_energy(e1, e2);
    }

    bool interaction2D::enabled() const { return m_enabled; }
    void interaction2D::enabled(const bool enabled) { m_enabled = enabled; }
}