#pragma once

#include "ppx/collision/collision2D.hpp"
#include "ppx/constraints/vconstraint2D.hpp"

namespace ppx
{
class friction_constraint2D final : public vconstraint10_2D
{
  public:
    friction_constraint2D(world2D &world, const collision2D *collision, std::size_t manifold_index);

    float max_impulse;

    float constraint_velocity() const override;
    void solve_velocities() override;
    void update(const collision2D *collision, const glm::vec2 &lanchor1, const glm::vec2 &nmtv);

  private:
    float m_friction;
    glm::vec2 m_nmtv;

    glm::vec2 direction() const override;
};
} // namespace ppx