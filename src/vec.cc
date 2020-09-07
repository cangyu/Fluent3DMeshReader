#include "vec.h"

/// From A to B
void delta(const VEC& na, const VEC& nb, VEC& dst)
{
    dst.x() = nb.x() - na.x();
    dst.y() = nb.y() - na.y();
    dst.z() = nb.z() - na.z();
}

double line_length(const VEC& na, const VEC& nb)
{
    double ret = 0.0;
    ret += std::pow(na.x() - nb.x(), 2);
    ret += std::pow(na.y() - nb.y(), 2);
    ret += std::pow(na.z() - nb.z(), 2);
    return std::sqrt(ret);
}

void line_center(const VEC& na, const VEC& nb, VEC& dst)
{
    dst.x() = 0.5 * (na.x() + nb.x());
    dst.y() = 0.5 * (na.y() + nb.y());
    dst.z() = 0.5 * (na.z() + nb.z());
}

double triangle_area(const VEC& na, const VEC& nb, const VEC& nc)
{
    const auto c = line_length(na, nb);
    const auto a = line_length(nb, nc);
    const auto b = line_length(nc, na);
    const auto p = 0.5 * (a + b + c);
    return std::sqrt(p * (p - a) * (p - b) * (p - c)); /// Heron's formula
}

void triangle_center(const VEC& na, const VEC& nb, const VEC& nc, VEC& dst)
{
    dst.x() = (na.x() + nb.x() + nc.x()) / 3.0;
    dst.y() = (na.y() + nb.y() + nc.y()) / 3.0;
    dst.z() = (na.z() + nb.z() + nc.z()) / 3.0;
}

void triangle_normal(const VEC& na, const VEC& nb, const VEC& nc, VEC& dst)
{
    VEC rab, rac;
    delta(na, nb, rab);
    delta(na, nc, rac);
    rab.cross(rac, dst); /// Take cross product to find normal direction
    dst.normalize();
}

double quadrilateral_area(const VEC& n1, const VEC& n2, const VEC& n3, const VEC& n4)
{
    return triangle_area(n1, n2, n3) + triangle_area(n1, n3, n4);
}

void quadrilateral_center(const VEC& n1, const VEC& n2, const VEC& n3, const VEC& n4, VEC& dst)
{
    const auto S123 = triangle_area(n1, n2, n3);
    const auto S134 = triangle_area(n1, n3, n4);

    VEC rc123, rc134;
    triangle_center(n1, n2, n3, rc123);
    triangle_center(n1, n3, n4, rc134);

    const auto alpha = S123 / (S123 + S134);
    const auto beta = 1.0 - alpha;

    rc123 *= alpha;
    rc134 *= beta;

    dst.x() = rc123.x() + rc134.x();
    dst.y() = rc123.y() + rc134.y();
    dst.z() = rc123.z() + rc134.z();
}

void quadrilateral_normal(const VEC& n1, const VEC& n2, const VEC& n3, const VEC& n4, VEC& dst)
{
    VEC ra, rb;
    delta(n2, n4, ra);
    delta(n1, n3, rb);
    ra.cross(rb, dst);
    dst.normalize();
}
