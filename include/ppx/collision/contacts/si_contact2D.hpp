#pragma once

#include "ppx/collision/contacts/si_friction2D.hpp"
#include "ppx/collision/contacts/contact_constraint2D.hpp"

namespace ppx
{
class si_contact2D final : public contact_constraint2D
{
  public:
    si_contact2D(world2D &world, const collision2D *collision, std::size_t manifold_index);

    float constraint_position() const override;
    float constraint_velocity() const override;

    void startup() override;
    void solve_velocities() override;
    void update(const collision2D *collision, std::size_t manifold_index) override;

  private:
    float m_restitution;
    float m_penetration;
    float m_pntr_correction = 0.f;
    glm::vec2 m_nmtv;
    si_friction2D m_friction;

    bool m_has_friction;
    bool m_is_adjusting_positions = false;

    float m_init_ctr_vel = 0.f;

    void update_position_data() override;
    void warmup() override;
    glm::vec2 direction() const override;
};
} // namespace ppx