#pragma once

#include <vector>
#include <algorithm>

#include "options.hpp"

extern shared_ptr<Options> _options;

struct Int3D {

  int64_t x, y, z;

  Int3D() : x(0), y(0), z(0) {}

  Int3D(int64_t x, int64_t y, int64_t z) : x(x), y(y), z(z) {}

  Int3D & operator=(const Int3D &val) {
    x = val.x;
    y = val.y;
    z = val.z;
    return *this;
  }

};

struct Double3D {

  double x, y, z;

  Double3D() : x(0.0), y(0.0), z(0.0) {}

  Double3D(double x, double y, double z) : x(x), y(y), z(z) {}

  Double3D & operator=(const Double3D &val) {
    x = val.x;
    y = val.y;
    z = val.z;
    return *this;
  }

};

class Alveolar {

 public:

   Alveolar() {}

   ~Alveolar() {}

   const std::vector<Int3D> &getEpiLocations() { return epiCellPositions3D; }

   const std::set<int64_t> &getEpiCellIds() { return epiCellPositions1D; }

 private:

   std::vector<Int3D> epiCellPositions3D;
   std::set<int64_t> epiCellPositions1D;

   void construct(const Int3D & pos) {
     // Cube dimensions for dim = 5 epi cells across are x,y,z = [-2 2]
      int dim = 5, idim = 2;
      epiCellPositions3D.clear();
      epiCellPositions1D.clear();
      for (int x = -idim; x <= idim; x++) {
       for (int y = -idim; y <= idim; y++) {
         for (int z = 0; z < dim; z++) {
           // Cells in the two x planes
           if (x == -idim || x == idim) {
             epiCellPositions3D.push_back({.x = x + pos.x,
               .y = y + pos.y,
               .z = z + pos.z});
             epiCellPositions1D.insert(pos.x
               + pos.y * gridSize.x
               + pos.z * gridSize.x * gridSize.y);
           }
           // Cells in the two y planes
           if (y == -idim || y == idim) {
             epiCellPositions3D.push_back({.x = x + pos.x,
               .y = y + pos.y,
               .z = z + pos.z});
             epiCellPositions1D.insert(pos.x
               + pos.y * gridSize.x
               + pos.z * gridSize.x * gridSize.y);
           }
           // Cells in the one z plane at bottom of alveoli
           if (z == dim - 1) {
             epiCellPositions3D.push_back({.x = x + pos.x,
               .y = y + pos.y,
               .z = z + pos.z});
             epiCellPositions1D.insert(pos.x
               + pos.y * gridSize.x
               + pos.z * gridSize.x * gridSize.y);
           }
         }
       }
      }
      //TODO rotate
   }

};
