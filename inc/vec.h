#ifndef VEC_H
#define VEC_H

#include <cstddef>
#include <cmath>

class VEC
{
private:
    double m_comp[3];

public:
    VEC() : m_comp{0.0, 0.0, 0.0} {}

    explicit VEC(double val) : m_comp{val, val, val} {}

    VEC(double v1, double v2, double v3) : m_comp{v1, v2, v3} {}

    VEC(const VEC &rhs) : m_comp{rhs.m_comp[0], rhs.m_comp[1], rhs.m_comp[2]} {}

    [[nodiscard]] double x() const { return m_comp[0]; }
    [[nodiscard]] double y() const { return m_comp[1]; }
    [[nodiscard]] double z() const { return m_comp[2]; }

    double &x() { return m_comp[0]; }
    double &y() { return m_comp[1]; }
    double &z() { return m_comp[2]; }

    double &at(size_t loc_idx)
    {
        return m_comp[loc_idx];
    }

    [[nodiscard]] double at(size_t loc_idx) const
    {
        return m_comp[loc_idx];
    }

    VEC &operator+=(const VEC &rhs)
    {
        x() += rhs.x();
        y() += rhs.y();
        z() += rhs.z();
        return *this;
    }

    VEC &operator-=(const VEC &rhs)
    {
        x() -= rhs.x();
        y() -= rhs.y();
        z() -= rhs.z();
        return *this;
    }

    VEC &operator*=(double a)
    {
        x() *= a;
        y() *= a;
        z() *= a;
        return *this;
    }

    VEC &operator/=(double a)
    {
        x() /= a;
        y() /= a;
        z() /= a;
        return *this;
    }

    VEC &operator=(const VEC &rhs)
    {
        if(this != &rhs)
        {
            x() = rhs.x();
            y() = rhs.y();
            z() = rhs.z();
        }
        return *this;
    }

    [[nodiscard]] double dot(const VEC &b) const
    {
        double ret = 0.0;
        ret += x() * b.x();
        ret += y() * b.y();
        ret += z() * b.z();
        return ret;
    }

    void cross(const VEC &b, VEC &ret) const
    {
        ret.x() = y() * b.z() - z() * b.y();
        ret.y() = z() * b.x() - x() * b.z();
        ret.z() = x() * b.y() - y() * b.x();
    }

    [[nodiscard]] double norm() const
    {
        double ret = 0.0;
        ret += std::pow(x(), 2);
        ret += std::pow(y(), 2);
        ret += std::pow(z(), 2);
        return std::sqrt(ret);
    }

    void normalize()
    {
        this->operator/=(norm());
    }
};

void delta(const VEC &na, const VEC &nb, VEC &dst);
double line_length(const VEC &na, const VEC &nb);
void line_center(const VEC &na, const VEC &nb, VEC &dst);
double triangle_area(const VEC &na, const VEC &nb, const VEC &nc);
void triangle_center(const VEC &na, const VEC &nb, const VEC &nc, VEC &dst);
void triangle_normal(const VEC &na, const VEC &nb, const VEC &nc, VEC &dst);
double quadrilateral_area(const VEC &n1, const VEC &n2, const VEC &n3, const VEC &n4);
void quadrilateral_center(const VEC &n1, const VEC &n2, const VEC &n3, const VEC &n4, VEC &dst);
void quadrilateral_normal(const VEC &n1, const VEC &n2, const VEC &n3, const VEC &n4, VEC &dst);

#endif
