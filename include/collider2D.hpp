#ifndef COLLIDER_HPP
#define COLLIDER_HPP

#include "entity2D.hpp"
#include "entity_ptr.hpp"
#include "constrain2D.hpp"
#include "quad_tree2D.hpp"
#include <vector>
#include <utility>

namespace phys
{
    class collider2D
    {
    public:
        collider2D(const std::vector<entity2D> &entities,
                   float stiffness = 1000.f,
                   float dampening = 10.f,
                   std::size_t allocations = 40);

        void add_entity_intervals(const const_entity_ptr &e); // TODO: Implement remove
        void solve_and_load_collisions(std::vector<float> &stchanges);

        float stiffness() const;
        float dampening() const;

        void stiffness(float stiffness);
        void dampening(float dampening);

        void update_quad_tree();

    private:
        struct interval
        {
        public:
            enum end
            {
                LOWER,
                HIGHER
            };

            interval(const const_entity_ptr &e, end end_type);

            const const_entity_ptr &entity() const;
            float value() const;
            end type() const;

        private:
            const_entity_ptr m_entity;
            end m_end;
        };

        struct collision
        {
            const_entity_ptr e1, e2;
            alg::vec2 touch1, touch2;
        };

        const std::vector<entity2D> &m_entities;
        std::vector<interval> m_intervals;
        quad_tree2D m_quad_tree;
        float m_stiffness, m_dampening;

        void sort_intervals();
        static bool collide(const const_entity_ptr &e1, const const_entity_ptr &e2, collision &c);
        std::vector<collision> brute_force() const;
        std::vector<collision> sort_and_sweep();
        std::vector<collision> quad_tree();

        void load_collisions(const std::vector<collision> &collisions,
                             std::vector<float> &stchanges) const;
        std::array<float, 6> forces_upon_collision(const collision &c) const;

        static bool gjk_epa(const const_entity_ptr &e1, const const_entity_ptr &e2, collision &col);

        static bool gjk(const geo::polygon2D &poly1, const geo::polygon2D &poly2, std::vector<alg::vec2> &simplex);
        static void line_case(const std::vector<alg::vec2> &simplex, alg::vec2 &dir);
        static bool triangle_case(std::vector<alg::vec2> &simplex, alg::vec2 &dir);

        static alg::vec2 epa(const geo::polygon2D &poly1, const geo::polygon2D &poly2, std::vector<alg::vec2> &simplex);
        static std::pair<alg::vec2, alg::vec2> touch_points(const geo::polygon2D &poly1,
                                                            const geo::polygon2D &poly2,
                                                            const alg::vec2 &mtv);
    };
}

#endif