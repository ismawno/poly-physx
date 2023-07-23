#ifndef PPX_QUAD_TREE2D_HPP
#define PPX_QUAD_TREE2D_HPP

#include "ppx/body2D.hpp"
#include "kit/memory/scope.hpp"
#include <memory>
#include <array>

namespace ppx
{
class quad_tree2D final
{
  public:
    quad_tree2D(const glm::vec2 &min, const glm::vec2 &max, std::size_t max_entities = 12, std::uint32_t depth = 0);

    void partitions(std::vector<const std::vector<const body2D *> *> &partitions) const;
    void insert(const body2D *bd);
    void clear();

    const geo::aabb2D &aabb() const;
    void aabb(const geo::aabb2D &aabb);

    std::size_t max_entities() const;
    void max_entities(std::size_t max_entities);

    bool partitioned() const;
    const std::vector<const body2D *> &entities() const;

    const std::array<kit::scope<quad_tree2D>, 4> &children() const;
    const quad_tree2D &child(std::size_t index) const;
    const quad_tree2D &operator[](std::size_t index) const;

    static std::uint32_t max_depth();
    static void max_depth(std::uint32_t max_depth);

    static float min_size();
    static void min_size(float min_size);

  private:
    std::array<kit::scope<quad_tree2D>, 4> m_children = {nullptr, nullptr, nullptr, nullptr}; // TL, TR, BL, BR
    geo::aabb2D m_aabb;
    std::size_t m_max_entities;
    std::uint32_t m_depth;
    bool m_partitioned = false, m_has_children = false;
    std::vector<const body2D *> m_entities;

    static std::uint32_t s_max_depth;
    static float s_min_size;

    bool full() const;
    bool rock_bottom() const;
    void create_children();
    void reset_children();
    void partition();
    void insert_to_children(const body2D *bd);

    quad_tree2D(quad_tree2D &&) = default;
    quad_tree2D &operator=(quad_tree2D &&) = default;
};
} // namespace ppx

#endif