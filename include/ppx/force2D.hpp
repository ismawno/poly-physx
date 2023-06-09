#ifndef FORCE2D_HPP
#define FORCE2D_HPP

#include "ppx/entity2D_set.hpp"

namespace ppx
{
    class force2D : public entity2D_set
    {
    public:
        using entity2D_set::entity2D_set;
        virtual ~force2D() = default;
        virtual std::pair<glm::vec2, float> force(const entity2D &e) const = 0;

        virtual float potential_energy(const entity2D &e) const { return 0.f; }
        float potential_energy() const;

        float energy(const entity2D &e) const;
        float energy() const;

    private:
        force2D(const force2D &) = delete;
        force2D &operator=(const force2D &) = delete;
    };
}

#endif