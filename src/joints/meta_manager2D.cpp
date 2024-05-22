#include "ppx/internal/pch.hpp"
#include "ppx/joints/joint_meta_manager2D.hpp"
#include "ppx/collision/resolution/sequential_impulses_resolution2D.hpp"
#include "ppx/world2D.hpp"

namespace ppx
{
template <IManager IM> bool joint_meta_manager2D<IM>::remove(std::size_t index)
{
    if (index >= this->m_elements.size())
        return false;
    this->m_elements.erase(this->m_elements.begin() + index);
    return true;
}

template <IManager IM> bool joint_meta_manager2D<IM>::remove(joint2D *joint)
{
    for (const auto &manager : this->m_elements)
        if (manager->remove(joint))
            return true;
    return false;
}

void actuator_meta_manager2D::solve()
{
    KIT_PERF_SCOPE("Joints solve")
    if (m_resolution)
        m_resolution->solve();
    for (const auto &manager : m_elements)
        if (manager->enabled)
            manager->solve();
}

void constraint_meta_manager2D::solve()
{
    KIT_PERF_SCOPE("Constraints solve")
    const std::size_t viters = world.constraints.velocity_iterations;
    const std::size_t piters = world.constraints.position_iterations;
    if (m_resolution)
        m_resolution->startup();

    for (const auto &manager : this->m_elements)
        if (manager->enabled)
            manager->startup();

    for (std::size_t i = 0; i < viters; i++)
    {
        if (m_resolution)
            m_resolution->solve_velocities();
        for (const auto &manager : this->m_elements)
            if (manager->enabled)
                manager->solve_velocities();
    }
    for (std::size_t i = 0; i < piters; i++)
    {
        bool solved = true;
        for (const auto &manager : this->m_elements)
            if (manager->enabled)
                solved &= manager->solve_positions();
        if (m_resolution)
            solved &= m_resolution->solve_positions();
        if (solved)
            break;
    }
}

template class joint_meta_manager2D<iactuator_manager2D>;
template class joint_meta_manager2D<iconstraint_manager2D>;

} // namespace ppx