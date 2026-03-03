// Python bindings for exact-geometry-lib using pybind11

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../geometry/Point.h"
#include "../geometry/CollisionPredicates.h"
#include "../geometry/IntersectionConstruction.h"
#include "../geometry/ConvexTriangulation.h"
#include "../geometry/SimplexQuadrature.h"

namespace py = pybind11;
using namespace dolfin;

PYBIND11_MODULE(geometry, m)
{
  m.doc() = "Python bindings for exact-geometry-lib";

  // -------------------------------------------------------------------------
  // Point class
  // -------------------------------------------------------------------------
  py::class_<Point>(m, "Point")
    .def(py::init<double, double, double>(),
         py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0)
    .def("x", &Point::x)
    .def("y", &Point::y)
    .def("z", &Point::z)
    .def("norm", &Point::norm)
    .def("squared_norm", &Point::squared_norm)
    .def("distance", &Point::distance)
    .def("squared_distance", &Point::squared_distance)
    .def("dot", &Point::dot)
    .def("cross", &Point::cross)
    .def("__getitem__", [](const Point& p, std::size_t i) { return p[i]; })
    .def("__setitem__", [](Point& p, std::size_t i, double v) { p[i] = v; })
    .def("__repr__", [](const Point& p) {
      return "Point(" + std::to_string(p.x()) + ", "
                      + std::to_string(p.y()) + ", "
                      + std::to_string(p.z()) + ")";
    })
    .def("__eq__", &Point::operator==)
    .def("__ne__", &Point::operator!=)
    .def("__add__", &Point::operator+)
    .def("__sub__",
         static_cast<Point (Point::*)(const Point&) const>(&Point::operator-))
    .def("__mul__", static_cast<Point (Point::*)(double) const>(&Point::operator*))
    .def("__truediv__", &Point::operator/);

  // -------------------------------------------------------------------------
  // CollisionPredicates
  // -------------------------------------------------------------------------
  py::class_<CollisionPredicates>(m, "CollisionPredicates")
    // segment-point
    .def_static("collides_segment_point_1d",
                &CollisionPredicates::collides_segment_point_1d)
    .def_static("collides_segment_point_2d",
                &CollisionPredicates::collides_segment_point_2d)
    .def_static("collides_segment_point_3d",
                &CollisionPredicates::collides_segment_point_3d)
    // segment-segment
    .def_static("collides_segment_segment_1d",
                &CollisionPredicates::collides_segment_segment_1d)
    .def_static("collides_segment_segment_2d",
                &CollisionPredicates::collides_segment_segment_2d)
    .def_static("collides_segment_segment_3d",
                &CollisionPredicates::collides_segment_segment_3d)
    // triangle-point
    .def_static("collides_triangle_point_2d",
                &CollisionPredicates::collides_triangle_point_2d)
    .def_static("collides_triangle_point_3d",
                &CollisionPredicates::collides_triangle_point_3d)
    // triangle-segment
    .def_static("collides_triangle_segment_2d",
                &CollisionPredicates::collides_triangle_segment_2d)
    .def_static("collides_triangle_segment_3d",
                &CollisionPredicates::collides_triangle_segment_3d)
    // triangle-triangle
    .def_static("collides_triangle_triangle_2d",
                &CollisionPredicates::collides_triangle_triangle_2d)
    .def_static("collides_triangle_triangle_3d",
                &CollisionPredicates::collides_triangle_triangle_3d)
    // tetrahedron-point
    .def_static("collides_tetrahedron_point_3d",
                &CollisionPredicates::collides_tetrahedron_point_3d)
    // tetrahedron-segment
    .def_static("collides_tetrahedron_segment_3d",
                &CollisionPredicates::collides_tetrahedron_segment_3d)
    // tetrahedron-triangle
    .def_static("collides_tetrahedron_triangle_3d",
                &CollisionPredicates::collides_tetrahedron_triangle_3d)
    // tetrahedron-tetrahedron
    .def_static("collides_tetrahedron_tetrahedron_3d",
                &CollisionPredicates::collides_tetrahedron_tetrahedron_3d);

  // -------------------------------------------------------------------------
  // IntersectionConstruction
  // -------------------------------------------------------------------------
  py::class_<IntersectionConstruction>(m, "IntersectionConstruction")
    // point-point
    .def_static("intersection_point_point_1d",
                &IntersectionConstruction::intersection_point_point_1d)
    .def_static("intersection_point_point_2d",
                &IntersectionConstruction::intersection_point_point_2d)
    .def_static("intersection_point_point_3d",
                &IntersectionConstruction::intersection_point_point_3d)
    // segment-point
    .def_static("intersection_segment_point_1d",
                &IntersectionConstruction::intersection_segment_point_1d)
    .def_static("intersection_segment_point_2d",
                &IntersectionConstruction::intersection_segment_point_2d)
    .def_static("intersection_segment_point_3d",
                &IntersectionConstruction::intersection_segment_point_3d)
    // triangle-point
    .def_static("intersection_triangle_point_2d",
                &IntersectionConstruction::intersection_triangle_point_2d)
    .def_static("intersection_triangle_point_3d",
                &IntersectionConstruction::intersection_triangle_point_3d)
    // tetrahedron-point
    .def_static("intersection_tetrahedron_point_3d",
                &IntersectionConstruction::intersection_tetrahedron_point_3d)
    // segment-segment
    .def_static("intersection_segment_segment_1d",
                &IntersectionConstruction::intersection_segment_segment_1d)
    .def_static("intersection_segment_segment_2d",
                &IntersectionConstruction::intersection_segment_segment_2d)
    .def_static("intersection_segment_segment_3d",
                &IntersectionConstruction::intersection_segment_segment_3d)
    // triangle-segment
    .def_static("intersection_triangle_segment_2d",
                &IntersectionConstruction::intersection_triangle_segment_2d)
    .def_static("intersection_triangle_segment_3d",
                &IntersectionConstruction::intersection_triangle_segment_3d)
    // tetrahedron-segment
    .def_static("intersection_tetrahedron_segment_3d",
                &IntersectionConstruction::intersection_tetrahedron_segment_3d)
    // triangle-triangle
    .def_static("intersection_triangle_triangle_2d",
                &IntersectionConstruction::intersection_triangle_triangle_2d)
    .def_static("intersection_triangle_triangle_3d",
                &IntersectionConstruction::intersection_triangle_triangle_3d)
    // tetrahedron-triangle
    .def_static("intersection_tetrahedron_triangle_3d",
                &IntersectionConstruction::intersection_tetrahedron_triangle_3d)
    // tetrahedron-tetrahedron
    .def_static("intersection_tetrahedron_tetrahedron_3d",
                &IntersectionConstruction::intersection_tetrahedron_tetrahedron_3d);

  // -------------------------------------------------------------------------
  // ConvexTriangulation
  // -------------------------------------------------------------------------
  py::class_<ConvexTriangulation>(m, "ConvexTriangulation")
    .def_static("triangulate",
                &ConvexTriangulation::triangulate)
    .def_static("triangulate_1d",
                &ConvexTriangulation::triangulate_1d)
    .def_static("triangulate_graham_scan_2d",
                &ConvexTriangulation::triangulate_graham_scan_2d)
    .def_static("triangulate_graham_scan_3d",
                &ConvexTriangulation::triangulate_graham_scan_3d)
    .def_static("selfintersects",
                &ConvexTriangulation::selfintersects);

  // -------------------------------------------------------------------------
  // SimplexQuadrature
  // -------------------------------------------------------------------------
  py::class_<SimplexQuadrature>(m, "SimplexQuadrature")
    .def(py::init<std::size_t, std::size_t>(),
         py::arg("tdim"), py::arg("order"))
    .def("compute_quadrature_rule",
         &SimplexQuadrature::compute_quadrature_rule)
    .def("compute_quadrature_rule_interval",
         &SimplexQuadrature::compute_quadrature_rule_interval)
    .def("compute_quadrature_rule_triangle",
         &SimplexQuadrature::compute_quadrature_rule_triangle)
    .def("compute_quadrature_rule_tetrahedron",
         &SimplexQuadrature::compute_quadrature_rule_tetrahedron)
    .def_static("compress",
                [](std::pair<std::vector<double>, std::vector<double>> qr,
                   std::size_t gdim, std::size_t order) {
                  return SimplexQuadrature::compress(qr, gdim, order);
                });
}
