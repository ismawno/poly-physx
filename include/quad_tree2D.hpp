#ifndef QUAD_TREE2D_HPP
#define QUAD_TREE2D_HPP

#include "vec2.hpp"
#include "entity_ptr.hpp"
#include <memory>

namespace phys
{
    class quad_tree2D
    {
    public:
        quad_tree2D() = delete;
        quad_tree2D(const alg::vec2 &pos, const alg::vec2 &dim, std::size_t max_entities = 5);

        void add_if_inside(const const_entity_ptr &e);
        void partitions(std::vector<const std::vector<const_entity_ptr> *> &partitions) const;
        void clear();

        bool full() const;
        const alg::vec2 &pos() const;
        const alg::vec2 &dim() const;
        std::size_t max_entities() const;
        bool partitioned() const;
        const std::vector<const_entity_ptr> &entities() const;

    private:
        std::unique_ptr<quad_tree2D> m_top_left, m_top_right, m_bottom_left, m_bottom_right;
        const alg::vec2 m_pos, m_dim;
        const std::size_t m_max_entities;
        bool m_partitioned, m_has_children;
        std::vector<const_entity_ptr> m_entities;

        void create_children();
        void partition();
        bool contains(const alg::vec2 &p) const;
        bool contains(const const_entity_ptr &e) const;
        void add_to_children(const const_entity_ptr &e);
    };
}

#endif