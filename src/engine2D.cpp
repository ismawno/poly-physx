#include "engine2D.hpp"
#include "ode2D.hpp"
#include "perf.hpp"
#include "rigid_bar2D.hpp"
#include <random>

namespace phys
{
    engine2D::engine2D(const rk::butcher_tableau &table,
                       const std::size_t allocations) : m_collider(&m_entities, 2 * allocations),
                                                        m_compeller(&m_entities, allocations),
                                                        m_integ(table)
    {
        m_entities.reserve(allocations);
        m_integ.state().reserve(6 * allocations);
    }
    engine2D::engine2D(const engine2D &eng) : m_entities(eng.m_entities),
                                              m_collider(&m_entities, 2 * eng.size(),
                                                         eng.m_collider.quad_tree().aabb().min(),
                                                         eng.m_collider.quad_tree().aabb().max()),
                                              m_compeller(&m_entities, eng.compeller().constraints().capacity()),
                                              m_integ(eng.m_integ),
                                              m_elapsed(eng.m_elapsed)
    {
        for (entity2D &e : m_entities)
        {
            e.m_state = &m_integ.state();
            m_collider.add_entity_intervals({&m_entities, e.index()});
        }

        m_springs.reserve(eng.m_springs.capacity());
        for (const spring2D &sp : eng.m_springs)
        {
            const const_entity2D_ptr e1 = {&m_entities, sp.e1().index()},
                                     e2 = {&m_entities, sp.e2().index()};
            if (sp.has_joints())
                m_springs.emplace_back(e1, e2, sp.joint1(), sp.joint2(), sp.length());
            else
                m_springs.emplace_back(e1, e2, sp.length());
        }

        for (const auto &ctr : eng.m_compeller.constraints())
        {
            const rigid_bar2D *rb = nullptr;
            try
            {
                rb = &dynamic_cast<const rigid_bar2D &>(*ctr);
            }
            catch (const std::bad_cast &e)
            {
                DBG_LOG("%s\n", e.what())
                continue;
            }
            const entity2D_ptr e1 = {&m_entities, rb->e1().index()},
                               e2 = {&m_entities, rb->e2().index()};
            if (rb->has_joints())
                m_compeller.add_constraint(std::make_shared<rigid_bar2D>(e1, e2, rb->joint1(), rb->joint2(), rb->length()));
            else
                m_compeller.add_constraint(std::make_shared<rigid_bar2D>(e1, e2, rb->length()));
        }
    }

    void engine2D::retrieve(const std::vector<float> &vars_buffer)
    {
        PERF_FUNCTION()
        for (std::size_t i = 0; i < m_entities.size(); i++)
            m_entities[i].retrieve(vars_buffer);
    }

    void engine2D::retrieve() { retrieve(m_integ.state().vars()); }

    bool engine2D::raw_forward(float &timestep)
    {
        const bool valid = m_integ.raw_forward(m_elapsed, timestep, *this, ode);
        register_forces_into_entities();
        reset_forces();
        retrieve();
        return valid;
    }
    bool engine2D::reiterative_forward(float &timestep, const std::size_t reiterations)
    {
        const bool valid = m_integ.reiterative_forward(m_elapsed, timestep, *this, ode, reiterations);
        register_forces_into_entities();
        reset_forces();
        retrieve();
        return valid;
    }
    bool engine2D::embedded_forward(float &timestep)
    {
        const bool valid = m_integ.embedded_forward(m_elapsed, timestep, *this, ode);
        register_forces_into_entities();
        reset_forces();
        retrieve();
        return valid;
    }

    void engine2D::load_velocities_and_added_forces(std::vector<float> &stchanges) const
    {
        PERF_FUNCTION()
        for (std::size_t i = 0; i < m_entities.size(); i++)
        {
            const std::size_t index = 6 * i;
            const alg::vec2 &vel = m_entities[i].vel();
            const float angvel = m_entities[i].angvel();
            stchanges[index] = vel.x;
            stchanges[index + 1] = vel.y;
            stchanges[index + 2] = angvel;
            if (m_entities[i].kinematic())
            {
                const alg::vec2 &force = m_entities[i].added_force();
                const float torque = m_entities[i].added_torque();
                load_force(stchanges, force, torque, index);
            }
        }
    }

    void engine2D::register_forces_into_entities() // TODO: Change name register to load wtf
    {
        const std::vector<float> step = m_integ.state().step();
        for (std::size_t i = 0; i < m_entities.size(); i++)
        {
            const std::size_t index = 6 * i;
            m_entities[i].m_force = {step[index + 3], step[index + 4]};
            m_entities[i].m_torque = step[index + 5];
        }
    }

    void engine2D::load_force(std::vector<float> &stchanges,
                              const alg::vec2 &force,
                              float torque,
                              std::size_t index)
    {
        stchanges[index + 3] += force.x;
        stchanges[index + 4] += force.y;
        stchanges[index + 5] += torque;
    }

    void engine2D::validate()
    {
        m_collider.validate();
        m_compeller.validate();
        for (auto it = m_springs.begin(); it != m_springs.end();)
            if (!it->try_validate())
                it = m_springs.erase(it);
            else
                ++it;
    }

    void engine2D::load_interactions_and_externals(std::vector<float> &stchanges) const
    {
        PERF_FUNCTION()
        for (const std::shared_ptr<const force2D> f : m_forces)
            if (f->enabled())
                for (const const_entity2D_ptr &e : f->entities())
                {
                    if (!e->kinematic())
                        continue;
                    const std::size_t index = 6 * e.index();
                    const auto [force, torque] = f->force(*e);
                    load_force(stchanges, force, torque, index);
                }
        for (const spring2D &s : m_springs)
        {
            const std::size_t index1 = 6 * s.e1().index(),
                              index2 = 6 * s.e2().index();
            const auto [force, t1, t2] = s.force();
            if (s.e1()->kinematic())
                load_force(stchanges, force, t1, index1);
            if (s.e2()->kinematic())
                load_force(stchanges, -force, t2, index2);
        }
        for (const std::shared_ptr<const interaction2D> i : m_interactions)
            for (const const_entity2D_ptr &e1 : i->entities())
            {
                if (!e1->kinematic())
                    continue;
                const std::size_t index = 6 * e1.index();
                for (const const_entity2D_ptr &e2 : i->entities())
                    if (e1 != e2)
                    {
                        const auto [force, torque] = i->force(*e1, *e2);
                        load_force(stchanges, force, torque, index);
                    }
            }
    }

    std::vector<float> engine2D::inverse_masses() const
    {
        PERF_FUNCTION()
        std::vector<float> inv_masses;
        inv_masses.reserve(3 * m_entities.size());
        for (std::size_t i = 0; i < m_entities.size(); i++)
        {
            const float inv_mass = 1.f / m_entities[i].mass(),
                        inv_inertia = 1.f / m_entities[i].inertia();
            inv_masses.insert(inv_masses.end(), {inv_mass, inv_mass, inv_inertia});
        }
        return inv_masses;
    }

    void engine2D::reset_forces()
    {
        for (entity2D &e : m_entities)
        {
            e.m_added_force = alg::vec2::zero;
            e.m_added_torque = 0.f;
        }
    }

    entity2D_ptr engine2D::add_entity(const alg::vec2 &pos,
                                      const alg::vec2 &vel,
                                      const float angpos,
                                      const float angvel,
                                      const float mass,
                                      const float charge,
                                      const std::vector<alg::vec2> &vertices,
                                      const bool kinematic)
    {
        entity2D &e = m_entities.emplace_back(pos, vel, angpos, angvel, mass, charge, vertices, kinematic);
        const entity2D_ptr e_ptr = {&m_entities, m_entities.size() - 1};

        rk::state &state = m_integ.state();
        e.m_index = m_entities.size() - 1;
        e.m_state = &state;

        state.append({pos.x, pos.y, angpos,
                      vel.x, vel.y, angvel});
        m_collider.add_entity_intervals(e_ptr);
        e.retrieve();

        DBG_LOG("Added entity with index %zu and id %zu.\n", e.m_index, e.m_id)
#ifdef DEBUG
        for (std::size_t i = 0; i < m_entities.size() - 1; i++)
            DBG_ASSERT(m_entities[i].m_id != e.m_id, "Added entity has the same id as entity with index %zu.\n", i)
#endif
        for (const add_callback &cb : m_on_entity_addition)
            cb(e_ptr);
        return e_ptr;
    }

    void engine2D::remove_entity(const std::size_t index)
    {
        DBG_ASSERT(index < m_entities.size(), "Index exceeds entity array bounds - index: %zu, size: %zu\n", index, m_entities.size())
        rk::state &state = m_integ.state();

        if (index == m_entities.size() - 1)
            m_entities.pop_back();
        else
        {
            m_entities[index] = m_entities.back();
            m_entities.pop_back();
            m_entities[index].m_index = index;
            m_entities[index].m_state = &state;
        }

        for (std::size_t i = 0; i < 6; i++)
            state[6 * index + i] = state[state.size() - 6 + i];
        state.resize(6 * m_entities.size());

        validate();
        m_collider.update_quad_tree();
        for (const remove_callback &cb : m_on_entity_removal)
            cb(index);
    }

    void engine2D::remove_entity(const entity2D &e) { remove_entity(e.index()); }

    void engine2D::add_force(const std::shared_ptr<force2D> &force) { m_forces.emplace_back(force); }
    void engine2D::add_interaction(const std::shared_ptr<interaction2D> &inter) { m_interactions.emplace_back(inter); }
    void engine2D::add_spring(const spring2D &spring) { m_springs.emplace_back(spring); }

    void engine2D::remove_force(const std::shared_ptr<force2D> &force)
    {
        m_forces.erase(std::remove(m_forces.begin(), m_forces.end(), force), m_forces.end());
    }
    void engine2D::remove_interaction(const std::shared_ptr<interaction2D> &inter)
    {
        m_interactions.erase(std::remove(m_interactions.begin(), m_interactions.end(), inter), m_interactions.end());
    }
    void engine2D::remove_spring(std::size_t index)
    {
        if (index >= m_springs.size())
            return;
        m_springs.erase(m_springs.begin() + index);
    }

    void engine2D::clear_entities()
    {
        for (std::size_t i = m_entities.size() - 1; i < m_entities.size(); i--)
            remove_entity(i);
    }
    void engine2D::clear_forces() { m_forces.clear(); }
    void engine2D::clear_interactions() { m_interactions.clear(); }
    void engine2D::clear_springs() { m_springs.clear(); }
    void engine2D::clear()
    {
        m_forces.clear();
        m_interactions.clear();
        m_springs.clear();
        m_compeller.clear_constraints();
        clear_entities();
    }

    void engine2D::write(ini::output &out) const
    {
        out.write("elapsed", m_elapsed);
        out.begin_section("tableau");
        m_integ.tableau().write(out);
        out.end_section();
        out.begin_section("collider");
        m_collider.write(out);
        out.end_section();

        std::string section = "entity";
        for (const entity2D &e : m_entities)
        {
            out.begin_section(section + std::to_string(e.index()));
            e.write(out);
            out.end_section();
        }

        section = "spring";
        std::size_t index = 0;
        for (const spring2D &sp : m_springs)
        {
            out.begin_section(section + std::to_string(index++));
            sp.write(out);
            out.end_section();
        }

        section = "rigid_bar";
        index = 0;
        for (const auto &ctr : m_compeller.constraints())
        {
            const rigid_bar2D *rb = nullptr;
            try
            {
                rb = &dynamic_cast<const rigid_bar2D &>(*ctr);
            }
            catch (const std::bad_cast &e)
            {
                DBG_LOG("%s\n", e.what())
                continue;
            }
            out.begin_section(section + std::to_string(index++));
            rb->write(out);
            out.end_section();
        }
    }

    void engine2D::read(ini::input &in)
    {
        clear_entities();
        m_elapsed = in.readf("elapsed");
        in.begin_section("tableau");
        rk::butcher_tableau tb;
        tb.read(in);
        m_integ.tableau(tb);
        in.end_section();
        in.begin_section("collider");
        m_collider.read(in);
        in.end_section();

        std::string section = "entity";
        std::size_t index = 0;
        while (true)
        {
            in.begin_section(section + std::to_string(index++));
            if (!in.contains_section())
            {
                in.end_section();
                break;
            }
            add_entity()->read(in);
            in.end_section();
        }

        section = "spring";
        index = 0;
        while (true)
        {
            in.begin_section(section + std::to_string(index++));
            if (!in.contains_section())
            {
                in.end_section();
                break;
            }
            const bool has_joints = (bool)in.readi("has_joints");
            const std::size_t idx1 = in.readi("e1"), idx2 = in.readi("e2");
            const entity2D_ptr e1 = (*this)[idx1], e2 = (*this)[idx2];

            if (has_joints)
            {
                alg::vec2 joint1, joint2;
                in.begin_section("joint1");
                joint1.read(in);
                in.end_section();
                in.begin_section("joint2");
                joint2.read(in);
                in.end_section();

                spring2D sp(e1, e2, joint1, joint2);
                sp.read(in);
                add_spring(sp);
            }
            else
            {
                spring2D sp(e1, e2);
                sp.read(in);
                add_spring(sp);
            }
            in.end_section();
        }

        section = "rigid_bar";
        index = 0;
        while (true)
        {
            in.begin_section(section + std::to_string(index++));
            if (!in.contains_section())
            {
                in.end_section();
                break;
            }
            const bool has_joints = (bool)in.readi("has_joints");
            const std::size_t idx1 = in.readi("e1"), idx2 = in.readi("e2");
            const entity2D_ptr e1 = (*this)[idx1], e2 = (*this)[idx2];

            if (has_joints)
            {
                alg::vec2 joint1, joint2;
                in.begin_section("joint1");
                joint1.read(in);
                in.end_section();
                in.begin_section("joint2");
                joint2.read(in);
                in.end_section();

                const auto rb = std::make_shared<rigid_bar2D>(e1, e2, joint1, joint2);
                rb->read(in);
                m_compeller.add_constraint(rb);
            }
            else
            {
                const auto rb = std::make_shared<rigid_bar2D>(e1, e2);
                rb->read(in);
                m_compeller.add_constraint(rb);
            }
            in.end_section();
        }
    }

    void engine2D::checkpoint() { m_checkpoint = std::make_tuple(m_elapsed, m_integ.state().vars(), m_entities); }
    void engine2D::revert()
    {
        const auto &[elapsed, vars, entities] = m_checkpoint;
        DBG_ASSERT(m_integ.state().vars().size() == vars.size() &&
                       m_entities.size() == entities.size(),
                   "Cannot revert to a checkpoint where the number of entities differ. Entities now: %zu, entities before: %zu.\n", m_entities.size(), entities.size())

        m_elapsed = elapsed;
        m_integ.state().vars(vars);
        m_entities = entities;
    }

    float engine2D::kinetic_energy() const
    {
        float ke = 0.f;
        for (const entity2D &e : m_entities)
            ke += e.kinetic_energy();
        return ke;
    }
    float engine2D::potential_energy() const
    {
        float pot = 0.f;
        for (const auto &force : m_forces)
            pot += force->potential_energy();
        for (const auto &inter : m_interactions)
            pot += inter->potential_energy();
        for (const spring2D &sp : m_springs)
            pot += sp.potential_energy();
        return pot;
    }
    float engine2D::energy() const { return kinetic_energy() + potential_energy(); }

    void engine2D::on_entity_addition(const add_callback &on_add) { m_on_entity_addition.emplace_back(on_add); }
    void engine2D::on_entity_removal(const remove_callback &on_remove) { m_on_entity_removal.emplace_back(on_remove); }

    const_entity2D_ptr engine2D::operator[](std::size_t index) const { return {&m_entities, index}; }
    entity2D_ptr engine2D::operator[](std::size_t index) { return {&m_entities, index}; }

    std::vector<const_entity2D_ptr> engine2D::operator[](const geo::aabb2D &aabb) const
    {
        std::vector<const_entity2D_ptr> in_area;
        in_area.reserve(m_entities.size() / 2);
        for (const entity2D &e : m_entities)
            if (e.aabb().overlaps(aabb))
                in_area.emplace_back(&m_entities, e.index());
        return in_area;
    }
    std::vector<entity2D_ptr> engine2D::operator[](const geo::aabb2D &aabb)
    {
        std::vector<entity2D_ptr> in_area;
        in_area.reserve(m_entities.size() / 2);
        for (const entity2D &e : m_entities)
            if (e.aabb().overlaps(aabb))
                in_area.emplace_back(&m_entities, e.index());
        return in_area;
    }

    const std::vector<std::shared_ptr<force2D>> &engine2D::forces() const { return m_forces; }
    const std::vector<std::shared_ptr<interaction2D>> &engine2D::interactions() const { return m_interactions; }
    const std::vector<spring2D> &engine2D::springs() const { return m_springs; }

    utils::vector_view<std::shared_ptr<force2D>> engine2D::forces() { return m_forces; }
    utils::vector_view<std::shared_ptr<interaction2D>> engine2D::interactions() { return m_interactions; }
    utils::vector_view<spring2D> engine2D::springs() { return m_springs; }

    const_entity2D_ptr engine2D::operator[](const alg::vec2 &point) const
    {
        const geo::aabb2D aabb = point;
        for (const entity2D &e : m_entities)
            if (e.aabb().overlaps(aabb))
                return {&m_entities, e.index()};
        return nullptr;
    }
    entity2D_ptr engine2D::operator[](const alg::vec2 &point)
    {
        const geo::aabb2D aabb = point;
        for (const entity2D &e : m_entities)
            if (e.aabb().overlaps(aabb))
                return {&m_entities, e.index()};
        return nullptr;
    }

    const std::vector<entity2D> &engine2D::entities() const { return m_entities; }
    utils::vector_view<entity2D> engine2D::entities() { return m_entities; }
    std::size_t engine2D::size() const { return m_entities.size(); }

    const rk::integrator &engine2D::integrator() const { return m_integ; }
    rk::integrator &engine2D::integrator() { return m_integ; }

    const compeller2D &engine2D::compeller() const { return m_compeller; }
    compeller2D &engine2D::compeller() { return m_compeller; }

    const collider2D &engine2D::collider() const { return m_collider; }
    collider2D &engine2D::collider() { return m_collider; }

    float engine2D::elapsed() const { return m_elapsed; }
}