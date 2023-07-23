#ifndef PPX_BODY2D_HPP
#define PPX_BODY2D_HPP

#include "geo/aabb2D.hpp"
#include "geo/polygon.hpp"
#include "geo/circle.hpp"
#include "rk/state.hpp"
#include "ppx/events/body_events.hpp"
#include "kit/interface/identifiable.hpp"
#include "kit/interface/indexable.hpp"
#include "kit/memory/track_ptr.hpp"
#include <variant>

namespace ppx
{
class body2D : public kit::identifiable<>, public kit::indexable
{
  public:
    using ptr = kit::track_ptr<body2D>;
    using const_ptr = kit::const_track_ptr<body2D>;

    enum class shape_type
    {
        POLYGON = 0,
        CIRCLE = 1
    };
    struct specs
    {
        glm::vec2 position{0.f}, velocity{0.f};
        float rotation = 0.f, angular_velocity = 0.f, mass = 1.f, charge = 1.f;
        kit::block_vector<glm::vec2> vertices = geo::polygon::box(5.f);
        float radius = 2.5f;
        bool kinematic = true;
        shape_type shape = shape_type::POLYGON;
        static specs from_body(const body2D &body);
    };

    body2D(const kit::block_vector<glm::vec2> &vertices, const glm::vec2 &position = glm::vec2(0.f),
           const glm::vec2 &velocity = glm::vec2(0.f), float rotation = 0.f, float angular_velocity = 0.f,
           float mass = 1.f, float charge = 1.f, bool kinematic = true);
    body2D(float radius, const glm::vec2 &position = glm::vec2(0.f), const glm::vec2 &velocity = glm::vec2(0.f),
           float rotation = 0.f, float angular_velocity = 0.f, float mass = 1.f, float charge = 1.f,
           bool kinematic = true);
    body2D(const glm::vec2 &position = glm::vec2(0.f), const glm::vec2 &velocity = glm::vec2(0.f), float rotation = 0.f,
           float angular_velocity = 0.f, float mass = 1.f, float charge = 1.f, bool kinematic = true);
    body2D(const specs &spc);

    void retrieve();
    void dispatch() const;
    float kinetic_energy() const;

    void add_force(const glm::vec2 &force);
    void add_torque(float torque);

    const glm::vec2 &added_force() const;
    float added_torque() const;

    const geo::shape2D &shape() const;

    template <typename T> const T &shape() const;

    template <typename T> const T *shape_if() const;

    void shape(const kit::block_vector<glm::vec2> &vertices);
    void shape(float radius);
    void shape(const geo::polygon &poly);
    void shape(const geo::circle &c);

    shape_type type() const;

    float inertia() const;
    float inverse_inertia() const;

    bool kinematic() const;
    void kinematic(bool kinematic);

    void translate(const glm::vec2 &dpos);
    void rotate(float dangle);

    const body_events &events() const;
    body_events &events();

    const glm::vec2 &position() const;
    const glm::vec2 &velocity() const;
    const glm::vec2 vel_at(const glm::vec2 &at) const;
    float rotation() const;
    float angular_velocity() const;
    float mass() const;
    float inverse_mass() const;
    float charge() const;

    void position(const glm::vec2 &position);
    void velocity(const glm::vec2 &velocity);
    void rotation(float rotation);
    void angular_velocity(float angular_velocity);
    void mass(float mass);
    void charge(float charge);

  private:
    std::variant<geo::polygon, geo::circle> m_shape;
    rk::state *m_state = nullptr;
    glm::vec2 m_vel{0.f}, m_added_force{0.f};
    body_events m_events;
    float m_angvel, m_added_torque = 0.f, m_mass, m_inv_mass, m_inertia, m_inv_inertia, m_charge;
    bool m_kinematic;

    geo::shape2D &get_shape();
    void retrieve(const std::vector<float> &vars_buffer);
    void compute_inertia(const geo::shape2D &sh);

    friend class world2D;
};
#ifdef KIT_USE_YAML_CPP
YAML::Emitter &operator<<(YAML::Emitter &out, const body2D &body);
#endif
} // namespace ppx

#ifdef KIT_USE_YAML_CPP
namespace YAML
{
template <> struct convert<ppx::body2D>
{
    static Node encode(const ppx::body2D &body);
    static bool decode(const Node &node, ppx::body2D &body);
};
} // namespace YAML
#endif

#endif