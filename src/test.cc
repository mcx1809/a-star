#include <array>
#include <cmath>
#include <iostream>

#include "a_star.h"

struct CostScalar {
  CostScalar() = default;

  CostScalar(double v) { this->v = v; }

  bool operator<(const CostScalar &o) const { return v < o.v; }

  CostScalar operator+(const CostScalar &o) const { return v + o.v; }

  double v;
};

struct Vertex3D {
  Vertex3D() = default;

  Vertex3D(double x, double y, double z) { xyz = {{x, y, z}}; }

  bool operator<(const Vertex3D &o) const {
    if (Distance2(o) < 1e-2) {
      return false;
    } else {
      return xyz < o.xyz;
    }
  }

  double Distance2(const Vertex3D &o) const {
    return std::pow(X() - o.X(), 2) + std::pow(Y() - o.Y(), 2) +
           std::pow(Z() - o.Z(), 2);
  }

  double Distance(const Vertex3D &o) const { return std::sqrt(Distance2(o)); }

  double X() const { return xyz[0]; }

  double Y() const { return xyz[1]; }

  double Z() const { return xyz[2]; }

  std::array<double, 3> xyz;
};

int main(int, char **) {
  auto start = Vertex3D(11.0, 0.0, 0.0);
  auto goal = Vertex3D(-15.0, 0.0, 0.0);

  AStar<Vertex3D, CostScalar> a_star;

  a_star.SetDistanceFn(
      [](const Vertex3D &v1, const Vertex3D &v2) { return v1.Distance(v2); });

  a_star.SetHeuristicFn([&](const Vertex3D &v) { return v.Distance(goal); });

  a_star.SetNeighborsFn(
      [](const Vertex3D &v,
         const std::function<void(const Vertex3D &, bool)> &add) {
        auto check_add = [&](double x, double y, double z) {
          auto new_v = Vertex3D(x, y, z);
          if (new_v < v || v < new_v) {
            auto r2 = new_v.Distance2(Vertex3D(0.0, 0.0, 0.0));
            if (r2 >= std::pow(10.0, 2) && r2 <= std::pow(20.0, 2)) {
              add(new_v, false);
            }
          }
        };

        for (auto i = -1; i <= 1; ++i) {
          for (auto j = -1; j <= 1; ++j) {
            for (auto k = -1; k <= 1; ++k) {
              check_add(v.X() + i, v.Y() + j, v.Z() + k);
            }
          }
        }
      });

  a_star.SetMaxNodes(1024 * 1024);

  a_star.Search(start, goal);

  for (const auto &i : a_star) {
    const auto &v = *std::get<0>(i);
    std::cout << v.X() << " " << v.Y() << " " << v.Z() << std::endl;
  }

  return 0;
}
