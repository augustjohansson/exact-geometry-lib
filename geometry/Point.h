#ifndef __POINT_H
#define __POINT_H

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <string>
#include <ostream>
#include <sstream>

namespace simpex
{
  // Optional: aggressive inlining hint for hot builds
#if defined(_MSC_VER)
#define DOLFIN_FORCEINLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define DOLFIN_FORCEINLINE inline __attribute__((always_inline))
#else
#define DOLFIN_FORCEINLINE inline
#endif

  class Point
  {
  public:
    // Construction / assignment: keep trivial where possible
    constexpr explicit Point(double x = 0.0, double y = 0.0, double z = 0.0) noexcept
      : _x{{x, y, z}} {}

    // Copy the first dim entries (dim expected <= 3). Keep asserts in debug.
    explicit Point(std::size_t dim, const double* x) noexcept : _x{{0.0, 0.0, 0.0}}
    {
      assert(dim <= 3);
      assert(x != nullptr || dim == 0);
      for (std::size_t i = 0; i < dim; ++i) _x[i] = x[i];
    }

    constexpr Point(const Point&) noexcept = default;
    constexpr Point(Point&&) noexcept = default;
    constexpr Point& operator=(const Point&) noexcept = default;
    constexpr Point& operator=(Point&&) noexcept = default;
    ~Point() = default;

    // Element access
    DOLFIN_FORCEINLINE double& operator[](std::size_t i) noexcept
    { assert(i < 3); return _x[i]; }

    DOLFIN_FORCEINLINE const double& operator[](std::size_t i) const noexcept
    { assert(i < 3); return _x[i]; }

    // Named coordinates
    DOLFIN_FORCEINLINE constexpr double x() const noexcept { return _x[0]; }
    DOLFIN_FORCEINLINE constexpr double y() const noexcept { return _x[1]; }
    DOLFIN_FORCEINLINE constexpr double z() const noexcept { return _x[2]; }

    // Raw pointer access
    DOLFIN_FORCEINLINE double* coordinates() noexcept { return _x.data(); }
    DOLFIN_FORCEINLINE const double* coordinates() const noexcept { return _x.data(); }

    DOLFIN_FORCEINLINE constexpr std::array<double, 3> array() const noexcept { return _x; }

    // Arithmetic
    DOLFIN_FORCEINLINE constexpr Point operator+(const Point& p) const noexcept
    { return Point(_x[0] + p._x[0], _x[1] + p._x[1], _x[2] + p._x[2]); }

    DOLFIN_FORCEINLINE constexpr Point operator-(const Point& p) const noexcept
    { return Point(_x[0] - p._x[0], _x[1] - p._x[1], _x[2] - p._x[2]); }

    DOLFIN_FORCEINLINE Point& operator+=(const Point& p) noexcept
    { _x[0] += p._x[0]; _x[1] += p._x[1]; _x[2] += p._x[2]; return *this; }

    DOLFIN_FORCEINLINE Point& operator-=(const Point& p) noexcept
    { _x[0] -= p._x[0]; _x[1] -= p._x[1]; _x[2] -= p._x[2]; return *this; }

    DOLFIN_FORCEINLINE constexpr Point operator-() const noexcept
    { return Point(-_x[0], -_x[1], -_x[2]); }

    DOLFIN_FORCEINLINE constexpr Point operator*(double a) const noexcept
    { return Point(a*_x[0], a*_x[1], a*_x[2]); }

    DOLFIN_FORCEINLINE Point& operator*=(double a) noexcept
    { _x[0] *= a; _x[1] *= a; _x[2] *= a; return *this; }

    DOLFIN_FORCEINLINE Point operator/(double a) const noexcept
    { return Point(_x[0]/a, _x[1]/a, _x[2]/a); }

    DOLFIN_FORCEINLINE Point& operator/=(double a) noexcept
    { _x[0] /= a; _x[1] /= a; _x[2] /= a; return *this; }

    DOLFIN_FORCEINLINE constexpr bool operator==(const Point& p) const noexcept
    { return _x[0] == p._x[0] && _x[1] == p._x[1] && _x[2] == p._x[2]; }

    DOLFIN_FORCEINLINE constexpr bool operator!=(const Point& p) const noexcept
    { return !(*this == p); }

    // Hot geometry ops: inline them
    DOLFIN_FORCEINLINE constexpr double squared_distance(const Point& p) const noexcept
    {
      const double dx = _x[0] - p._x[0];
      const double dy = _x[1] - p._x[1];
      const double dz = _x[2] - p._x[2];
      return dx*dx + dy*dy + dz*dz;
    }

    DOLFIN_FORCEINLINE double distance(const Point& p) const noexcept
    { return std::sqrt(squared_distance(p)); }

    DOLFIN_FORCEINLINE double norm() const noexcept
    { return std::sqrt(_x[0]*_x[0] + _x[1]*_x[1] + _x[2]*_x[2]); }

    DOLFIN_FORCEINLINE constexpr double squared_norm() const noexcept
    { return _x[0]*_x[0] + _x[1]*_x[1] + _x[2]*_x[2]; }

    DOLFIN_FORCEINLINE constexpr Point cross(const Point& p) const noexcept
    {
      return Point(
		   _x[1]*p._x[2] - _x[2]*p._x[1],
		   _x[2]*p._x[0] - _x[0]*p._x[2],
		   _x[0]*p._x[1] - _x[1]*p._x[0]
		   );
    }

    DOLFIN_FORCEINLINE constexpr double dot(const Point& p) const noexcept
    { return _x[0]*p._x[0] + _x[1]*p._x[1] + _x[2]*p._x[2]; }

    Point rotate(const Point& a, double theta) const
    {
      assert(std::abs(a.norm() - 1.0) < 3.0e-16);

      const double s = std::sin(theta);
      const double c = std::cos(theta);
      const double one_c = 1.0 - c;

      const double ax = a[0], ay = a[1], az = a[2];
      const double x  = _x[0], y  = _x[1], z  = _x[2];

      const double dot = ax*x + ay*y + az*z;

      const double cx = ay*z - az*y;
      const double cy = az*x - ax*z;
      const double cz = ax*y - ay*x;

      return Point(x*c + cx*s + ax*dot*one_c,
		   y*c + cy*s + ay*dot*one_c,
		   z*c + cz*s + az*dot*one_c);
    }
      
    std::string str(bool verbose=false) const
    {
      std::stringstream ss;
      ss << x() << ' ' << y() << ' ' << z();
      return ss.str();
    }

  private:
    std::array<double, 3> _x;
  };

  DOLFIN_FORCEINLINE constexpr Point operator*(double a, const Point& p) noexcept
  { return p * a; }

  inline std::ostream& operator<<(std::ostream& stream, const Point& point)
  { stream << point.str(false); return stream; }

} 

#endif
