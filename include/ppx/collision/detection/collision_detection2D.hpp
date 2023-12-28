#pragma once

#include "ppx/body2D.hpp"
#include "ppx/collision/collision2D.hpp"

#include "ppx/collision/manifold/manifold_algorithms2D.hpp"
#include "kit/memory/scope.hpp"
#include "kit/utility/utils.hpp"

#ifndef PPX_THREAD_COUNT
#define PPX_THREAD_COUNT 8
#endif

namespace ppx
{
class world2D;

class collision_detection2D
{
  public:
    virtual ~collision_detection2D() = default;

    world2D *world = nullptr;
    float epa_threshold = 1.e-3f;

    template <typename T = cc_manifold_algorithm2D> const T *cc_manifold_algorithm() const
    {
        return kit::const_get_casted_raw_ptr<T>(m_cc_manifold);
    }
    template <typename T = cc_manifold_algorithm2D> T *cc_manifold_algorithm()
    {
        return kit::get_casted_raw_ptr<T>(m_cc_manifold);
    }
    template <typename T = cp_manifold_algorithm2D> const T *cp_manifold_algorithm() const
    {
        return kit::const_get_casted_raw_ptr<T>(m_cp_manifold);
    }
    template <typename T = cp_manifold_algorithm2D> T *cp_manifold_algorithm()
    {
        return kit::get_casted_raw_ptr<T>(m_cp_manifold);
    }
    template <typename T = pp_manifold_algorithm2D> const T *pp_manifold_algorithm() const
    {
        return kit::const_get_casted_raw_ptr<T>(m_pp_manifold);
    }
    template <typename T = pp_manifold_algorithm2D> T *pp_manifold_algorithm()
    {
        return kit::get_casted_raw_ptr<T>(m_pp_manifold);
    }

    template <typename T, class... ManifArgs> const T *set_cc_manifold_algorithm(ManifArgs &&...args)
    {
        static_assert(std::is_base_of_v<cc_manifold_algorithm2D, T>,
                      "Manifold algorithm must inherit from cc_manifold_algorithm2D");

        auto manalg = kit::make_scope<T>(std::forward<ManifArgs>(args)...);
        T *ptr = manalg.get();

        m_cc_manifold = std::move(manalg);
        return ptr;
    }
    template <typename T, class... ManifArgs> const T *set_cp_manifold_algorithm(ManifArgs &&...args)
    {
        static_assert(std::is_base_of_v<cp_manifold_algorithm2D, T>,
                      "Manifold algorithm must inherit from cp_manifold_algorithm2D");
        auto manalg = kit::make_scope<T>(std::forward<ManifArgs>(args)...);
        T *ptr = manalg.get();

        m_cp_manifold = std::move(manalg);
        return ptr;
    }
    template <typename T, class... ManifArgs> const T *set_pp_manifold_algorithm(ManifArgs &&...args)
    {
        static_assert(std::is_base_of_v<pp_manifold_algorithm2D, T>,
                      "Manifold algorithm must inherit from pp_manifold_algorithm2D");
        auto manalg = kit::make_scope<T>(std::forward<ManifArgs>(args)...);
        T *ptr = manalg.get();

        m_pp_manifold = std::move(manalg);
        return ptr;
    }

    const std::vector<collision2D> &detect_collisions_cached();
    void clear_cached_collisions();

    const std::vector<collision2D> &collisions() const;

    virtual void on_attach()
    {
    }

    void inherit(collision_detection2D &coldet);

  protected:
    std::vector<collision2D> m_collisions;
#ifdef PPX_MULTITHREADED
    std::array<std::vector<collision2D>, PPX_THREAD_COUNT> m_mt_collisions;
#endif

    collision2D generate_collision(body2D &body1, body2D &body2) const;
    void cc_narrow_collision_check(body2D &body1, body2D &body2, collision2D &collision) const;
    void cp_narrow_collision_check(body2D &body1, body2D &body2, collision2D &collision) const;
    void pp_narrow_collision_check(body2D &body1, body2D &body2, collision2D &collision) const;

    void try_enter_or_stay_callback(const collision2D &c) const;
    void try_exit_callback(body2D &body1, body2D &body2) const;

  private:
    kit::scope<cc_manifold_algorithm2D> m_cc_manifold;
    kit::scope<cp_manifold_algorithm2D> m_cp_manifold;
    kit::scope<pp_manifold_algorithm2D> m_pp_manifold;

    virtual void detect_collisions() = 0;
};

} // namespace ppx
