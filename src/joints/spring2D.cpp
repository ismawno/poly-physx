#include "ppx/internal/pch.hpp"
#include "ppx/joints/spring2D.hpp"
#include "ppx/world2D.hpp"
#include "kit/utility/utils.hpp"

namespace ppx
{

spring2D::spring2D(world2D &world, const specs &spc)
    : joint2D(world, world.bodies.ptr(spc.bindex1), world.bodies.ptr(spc.bindex2), spc.ganchor1, spc.ganchor2),
      worldref2D(world), stiffness(spc.props.stiffness), damping(spc.props.damping), length(spc.props.length),
      non_linear_terms(spc.props.non_linear_terms), non_linear_contribution(spc.props.non_linear_contribution)
{
}

spring2D::const_ptr spring2D::as_ptr() const
{
    return world.joints.manager<spring2D>()->ptr(index);
}
spring2D::ptr spring2D::as_ptr()
{
    return world.joints.manager<spring2D>()->ptr(index);
}

glm::vec2 spring2D::non_linear_displacement(const glm::vec2 &displacement) const
{
    glm::vec2 nl_term = displacement;
    glm::vec2 nl_cummulative = displacement;
    float decay = 16.f;
    for (std::uint32_t term = 0; term < non_linear_terms; term++)
    {
        nl_term *= displacement * displacement;
        nl_cummulative += nl_term / decay;
        decay *= decay;
    }
    return nl_cummulative * non_linear_contribution;
}

glm::vec4 spring2D::force() const
{
    KIT_ASSERT_ERROR(stiffness >= 0.f, "Stiffness must be non-negative: {0}", stiffness)
    KIT_ASSERT_ERROR(damping >= 0.f, "Damping must be non-negative: {0}", damping)
    KIT_ASSERT_ERROR(length >= 0.f, "Length must be non-negative: {0}", length)
    KIT_ASSERT_ERROR(non_linear_contribution >= 0.f, "Non linear contribution must be non-negative: {0}",
                     non_linear_contribution)

    const glm::vec2 ga1 = ganchor1();
    const glm::vec2 ga2 = ganchor2();

    const glm::vec2 offset1 = ga1 - m_body1->centroid();
    const glm::vec2 offset2 = ga2 - m_body2->centroid();

    const glm::vec2 relpos = ga2 - ga1;
    const glm::vec2 direction = glm::normalize(relpos);
    const glm::vec2 relvel = direction * glm::dot(m_body2->gvelocity_at_centroid_offset(offset2) -
                                                      m_body1->gvelocity_at_centroid_offset(offset1),
                                                  direction);
    const glm::vec2 vlen = length * direction;

    const glm::vec2 displacement = relpos - vlen;
    const glm::vec2 force =
        stiffness * (non_linear_terms != 0 ? non_linear_displacement(displacement) : displacement) + damping * relvel;

    const float torque1 = kit::cross2D(offset1, force);
    const float torque2 = kit::cross2D(force, offset2);
    return {force, torque1, torque2};
}

float spring2D::kinetic_energy() const
{
    return m_body1->kinetic_energy() + m_body2->kinetic_energy();
}
float spring2D::potential_energy() const
{
    const float dist = glm::distance(ganchor1(), ganchor2()) - length;
    return 0.5f * stiffness * dist * dist;
}
float spring2D::energy() const
{
    return kinetic_energy() + potential_energy();
}

void spring2D::solve()
{
    const glm::vec4 f = force();
    m_body1->apply_simulation_force(glm::vec2(f));
    m_body2->apply_simulation_force(-glm::vec2(f));

    m_body1->apply_simulation_torque(f.z);
    m_body2->apply_simulation_torque(f.w);
}
} // namespace ppx