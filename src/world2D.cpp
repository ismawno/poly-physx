#include "ppx/internal/pch.hpp"
#include "ppx/world2D.hpp"
#include "ppx/behaviours/behaviour2D.hpp"
#include "ppx/joints/distance_joint2D.hpp"

#include "ppx/collision/detection/quad_tree_detection2D.hpp"
#include "ppx/collision/resolution/spring_driven_resolution2D.hpp"
#include "ppx/collision/resolution/constraint_driven_resolution2D.hpp"

#include <cstring>
#ifdef DEBUG
#include "kit/missing/fenv.h"
#endif

namespace ppx
{
bool world2D::step()
{
    pre_step_preparation();
    const bool valid = integrator.raw_forward(*this);
    post_step_setup();
    return valid;
}

void world2D::pre_step_preparation()
{
#ifdef DEBUG
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

    m_timestep_ratio = kit::approaches_zero(integrator.ts.value) ? 1.f : m_previous_timestep / integrator.ts.value;
    collisions.detection()->clear_cached_collisions();

    bodies.send_data_to_state(integrator.state);
}
void world2D::post_step_setup()
{
    bodies.reset_impulse_forces();
    bodies.retrieve_data_from_state_variables(integrator.state.vars());
    m_previous_timestep = integrator.ts.value;

#ifdef DEBUG
    fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
}

float world2D::timestep_ratio() const
{
    return m_timestep_ratio;
}

std::vector<float> world2D::create_state_derivative() const
{
    KIT_PERF_FUNCTION()
    std::vector<float> state_derivative(6 * bodies.size(), 0.f);

    for (const body2D &body : bodies)
    {
        const std::size_t index = 6 * body.index;
        const glm::vec2 &velocity = body.velocity;
        const float angular_velocity = body.angular_velocity;
        state_derivative[index] = velocity.x;
        state_derivative[index + 1] = velocity.y;
        state_derivative[index + 2] = angular_velocity;

        const glm::vec2 &accel = body.force() * body.inv_mass();
        const float angaccel = body.torque() * body.inv_inertia();
        state_derivative[index + 3] = accel.x;
        state_derivative[index + 4] = accel.y;
        state_derivative[index + 5] = angaccel;
    }
    return state_derivative;
}

void world2D::validate()
{
    bodies.validate();
    constraints.validate();
    behaviours.validate();
    springs.validate();
}

float world2D::kinetic_energy() const
{
    float ke = 0.f;
    for (const body2D &body : bodies)
        ke += body.kinetic_energy();
    return ke;
}
float world2D::potential_energy() const
{
    float pot = 0.f;
    for (const auto &bhv : behaviours)
        if (bhv->enabled)
            pot += bhv->potential_energy();
    for (const spring2D &sp : springs)
        pot += sp.potential_energy();
    return pot;
}
float world2D::energy() const
{
    return kinetic_energy() + potential_energy();
}

std::vector<float> world2D::operator()(const float time, const float timestep, const std::vector<float> &vars)
{
    KIT_PERF_FUNCTION()
    KIT_ASSERT_CRITICAL(
        vars.size() == 6 * bodies.size(),
        "State vector size must be exactly 6 times greater than the body array size - vars: {0}, body array: {1}",
        vars.size(), bodies.size())

    bodies.reset_simulation_forces();
    bodies.retrieve_data_from_state_variables(vars);

    bodies.apply_impulse_and_persistent_forces();
    behaviours.apply_forces();
    springs.apply_forces();

    if (collisions.enabled)
        collisions.solve();

    bodies.prepare_constraint_velocities();
    constraints.solve();
    return create_state_derivative();
}

#ifdef KIT_USE_YAML_CPP
YAML::Node world2D::serializer::encode(const world2D &world) const
{
    YAML::Node node;
    node["Integrator"] = world.integrator;

    YAML::Node nc = node["Collision"];
    YAML::Node ndet = nc["Detection"];

    ndet["Method"] = (int)world.collisions.detection_method();

    if (world.collisions.detection_method() == collision_manager2D::detection_type::QUAD_TREE)
        ndet["Force square"] = world.collisions.detection<quad_tree_detection2D>()->force_square_shape;

    YAML::Node nqt = ndet["Quad tree"];
    nqt["Max bodies"] = quad_tree::max_bodies;
    nqt["Max depth"] = quad_tree::max_depth;
    nqt["Min size"] = quad_tree::min_size;

    YAML::Node nres = nc["Resolution"];
    nres["Method"] = (int)world.collisions.resolution_method();

    switch (world.collisions.resolution_method())
    {
    case collision_manager2D::resolution_type::SPRING_DRIVEN: {
        auto resol = world.collisions.resolution<spring_driven_resolution2D>();
        nres["Rigidity"] = resol->rigidity;
        nres["Normal damping"] = resol->normal_damping;
        nres["Tangent damping"] = resol->tangent_damping;
        break;
    }
    case collision_manager2D::resolution_type::CONSTRAINT_DRIVEN: {
        auto resol = world.collisions.resolution<constraint_driven_resolution2D>();
        nres["Friction"] = resol->friction;
        nres["Restitution"] = resol->restitution;
        break;
    }
    case collision_manager2D::resolution_type::CUSTOM:
        break;
    }

    YAML::Node nctrm = node["Constraint params"];
    nctrm["Iterations"] = world.constraints.iterations;
    nctrm["Warmup"] = world.constraints.warmup;
    nctrm["Position corrections"] = world.constraints.position_corrections;

    for (const body2D &body : world.bodies)
        node["Bodies"].push_back(body);

    for (const spring2D &sp : world.springs)
        node["Springs"].push_back(sp);
    for (const auto &ctr : world.constraints)
    {
        YAML::Node child;
        child[ctr->name] = *ctr;
        node["Constraints"].push_back(child);
    }

    for (const auto &bhv : world.behaviours)
        node["Behaviours"][bhv->id] = *bhv;
    return node;
}
bool world2D::serializer::decode(const YAML::Node &node, world2D &world) const
{
    if (!node.IsMap() || node.size() < 3)
        return false;

    world.bodies.clear();
    world.integrator = node["Integrator"].as<rk::integrator<float>>();
    world.integrator.state.clear();

    const YAML::Node nc = node["Collision"];
    const YAML::Node ndet = nc["Detection"];

    auto det_type = (collision_manager2D::detection_type)ndet["Method"].as<int>();
    world.collisions.detection(det_type);
    if (det_type == collision_manager2D::detection_type::QUAD_TREE)
        world.collisions.detection<quad_tree_detection2D>()->force_square_shape = ndet["Force square"].as<bool>();

    const YAML::Node nqt = ndet["Quad tree"];
    quad_tree::max_bodies = nqt["Max bodies"].as<std::size_t>();
    quad_tree::max_depth = nqt["Max depth"].as<std::uint32_t>();
    quad_tree::min_size = nqt["Min size"].as<float>();

    const YAML::Node nres = nc["Resolution"];
    auto res_type = (collision_manager2D::resolution_type)nres["Method"].as<int>();
    switch (res_type)
    {
    case collision_manager2D::resolution_type::SPRING_DRIVEN: {
        world.collisions.set_resolution<spring_driven_resolution2D>(
            nres["Rigidity"].as<float>(), nres["Normal damping"].as<float>(), nres["Tangent damping"].as<float>());
        break;
    }
    case collision_manager2D::resolution_type::CONSTRAINT_DRIVEN: {
        world.collisions.set_resolution<constraint_driven_resolution2D>(nres["Friction"].as<float>(),
                                                                        nres["Restitution"].as<float>());
        break;
    case collision_manager2D::resolution_type::CUSTOM:
        break;
    }
    }

    const YAML::Node nctrm = node["Constraint params"];
    world.constraints.iterations = nctrm["Iterations"].as<std::uint32_t>();
    world.constraints.warmup = nctrm["Warmup"].as<bool>();
    world.constraints.position_corrections = nctrm["Position corrections"].as<bool>();

    if (node["Bodies"])
        for (const YAML::Node &n : node["Bodies"])
            world.bodies.add(n.as<body2D>());

    if (node["Springs"])
        for (const YAML::Node &n : node["Springs"])
        {
            spring2D sp;
            sp.world = &world;
            n.as<spring2D>(sp);
            world.springs.add(sp);
        }

    if (node["Constraints"])
        for (const YAML::Node &n : node["Constraints"])
            if (n["Distance"])
            {
                distance_joint2D dj;
                dj.world = &world;
                n["Distance"].as<distance_joint2D>(dj);
                world.constraints.add<distance_joint2D>(dj);
            }

    if (node["Behaviours"])
        for (auto it = node["Behaviours"].begin(); it != node["Behaviours"].end(); ++it)
        {
            const auto bhv = world.behaviours.from_name<behaviour2D>(it->first.as<std::string>());
            if (bhv)
                it->second.as<behaviour2D>(*bhv);
        }

    return true;
}
#endif
} // namespace ppx