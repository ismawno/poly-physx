#include "ppx/internal/pch.hpp"
#include "ppx/utility/identifiable.hpp"

namespace ppx
{
identifiable::identifiable(const uuid id) : m_uuid(id)
{
}

uuid identifiable::id() const
{
    return m_uuid;
}
void identifiable::id(const uuid id)
{
    m_uuid = id;
}
} // namespace ppx