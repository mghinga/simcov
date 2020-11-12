#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "upcxx_utils.hpp"
#include "utils.hpp"
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

struct Level {

  int64_t L = 0, d = 0, count = 0;
  double bAngle = 0.0, gAngle = 0.0;
  Int3D centroid;
  Double3D direction;

};

class Lung {

 public:

  Lung(shared_ptr<Random> rnd_gen) {
    this->rnd_gen = rnd_gen;
    gridSize.x = _options->dimensions[0];
    gridSize.y = _options->dimensions[1];
    gridSize.z = _options->dimensions[2];
    if (_options->lung_model_type == "algorithmic") { // Algorithmic
      loadEstimatedParameters();
      // Draw root
      Level lvl = levels.at(0);
      Int3D child = constructSegment({gridSize.x/2 + 3*lvl.d, gridSize.y/2, 0},
        lvl,
        0.0);
      // Recursively build tree
      construct(child, 1, lvl.bAngle);
      // // Build alveolus structures
      // for (Int3D leaf : leafs) {
      // //TODO Int3D leaf(0, 0, 0);
      //   Int3D left(leaf.x + 2, leaf.y + 2, leaf.z);
      //   constructAlveoli(left);
      //   Int3D right(leaf.x + 7, leaf.y + 7, leaf.z);
      //   constructAlveoli(right);
      // }
      SLOG("Number of alveoli ", 2 * leafs.size(), "\n");
    } else if (_options->lung_model_type == "empirical") { // Empirical
      loadEmpiricalData();
      // Draw all segments
      for (int i = 0; i < levels.size(); i++) {
        if (levels.at(i).L < 1 || levels.at(i).d < 1) {
          skipped++;
          continue;
        }
        constructSegment(levels.at(i));
      }
      SLOG(skipped, " skipped airway segments\n");
    } else {
      SDIE("Invalid lung model type ", _options->lung_model_type,
           " should be 'empirical' or 'algorithmic'");
    }
  }

  ~Lung() {}

  const std::set<int64_t> &getAirwayEpiCellIds() { return airwayEpiCellPositions1D; }

  const std::set<int64_t> &getAlveoliEpiCellIds() { return alveoliEpiCellPositions1D; }

  const std::vector<Int3D> &getEpiLocations() { return positions; }

 private:

   int skipped = 0; // Bauer et al 2019
   double scale = 10;// 1 / 5e-3; // Convert cm to um
   std::vector<Int3D> leafs;
   std::vector<Level> levels;
   std::set<int64_t> alveoliEpiCellPositions1D;
   std::set<int64_t> airwayEpiCellPositions1D;
   std::vector<Int3D> positions;
   Int3D gridSize;
   shared_ptr<Random> rnd_gen;

   Int3D rotate(const Int3D & vec,
     const Double3D & axis,
     const double & angle) {
    /**
    * http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
    */
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);
    double x = vec.x, y = vec.y, z = vec.z;
    //double u = 0, v = 1, w = 0;
    double u = axis.x, v = axis.y, w = axis.z;
    double dotProduct = (x * u) + (y * v) + (z * w);
    double newX = u * dotProduct * (1 - cosAngle) + x * cosAngle
      + (-w * y + v * z) * sinAngle;
    double newY = v * dotProduct * (1 - cosAngle) + y * cosAngle
      + (w * x - u * z) * sinAngle;
    double newZ = w * dotProduct * (1 - cosAngle) + z * cosAngle
      + (-v * x + u * y) * sinAngle;
    Int3D rval((int64_t) round(newX),
      (int64_t) round(newY),
      (int64_t) round(newZ));
    return rval;
   }

   void construct(const Int3D & root,
     int generation,
     double previousBranchAngle) {
       if (generation > (levels.size() - 1)) {
         leafs.push_back({.x = root.x, .y = root.y, .z = root.z});
         return;
       }
       double rotateZ = (generation < 7) ? 0.0 : rnd_gen->get(0,180)*M_PI/180;
       // Draw left child bronchiole
       Level lvl = levels.at(generation);
       lvl.bAngle = previousBranchAngle - lvl.bAngle;
       Int3D child = constructSegment(root, lvl, rotateZ);
       construct(child, generation + 1, lvl.bAngle);
       // Draw right child bronchiole
       lvl = levels.at(generation);
       lvl.bAngle = previousBranchAngle + lvl.bAngle;
       if (lvl.count > 1) {
         lvl.gAngle = -lvl.gAngle;
         child = constructSegment(root, lvl, rotateZ);
         construct(child, generation + 1, lvl.bAngle);
       }
     }

   Int3D constructSegment(const Int3D &root,
     const Level &level,
     double rotateZ) {
     std::vector<Int3D> positions;
     // Build cylinder at origin along y-axis
     int64_t radius = level.d / 2;
     double az = 0;
     double inc = M_PI / 2;
     for (int64_t z = 0; z <= level.L; z++) {
       if (_options->lung_model_is_skeleton) {
         positions.push_back({.x = 0, .y = 0, .z = z});
       } else {
         for (az = 0; az < 2 * M_PI; az += M_PI / 180) {
           int64_t x = (int64_t) round(radius * sin(inc) * cos(az));
           int64_t y = (int64_t) round(radius * sin(inc) * sin(az));
           positions.push_back({.x = x, .y = y, .z = z});
         }
       }
     }
     // Treat as positional vectors and apply rotations and translate
     for (int64_t i = 0; i < positions.size(); i++) {
       // Rotate y
       Int3D newPosition0 = rotate(positions.at(i),
        {.x = 0.0, .y = 1.0, .z = 0.0},
        level.bAngle);
       // Rotate z
       Int3D newPosition1 = rotate(newPosition0,
        {.x = 0.0, .y = 0.0, .z = 1.0},
        rotateZ);
       newPosition1.x += root.x;
       newPosition1.y += root.y;
       newPosition1.z += root.z;
       positions.at(i) = newPosition1;
       // Verify new location is within grid
       if ((0 <= newPosition1.x && newPosition1.x < gridSize.x)
        && (0 <= newPosition1.y && newPosition1.y < gridSize.y)
        && (0 <= newPosition1.z && newPosition1.z < gridSize.z)) {
          airwayEpiCellPositions1D.insert(newPosition1.x
            + newPosition1.y * gridSize.x
            + newPosition1.z * gridSize.x * gridSize.y);
       }
     }
     // Return root for next generation
     Int3D base0 = rotate({.x = 0, .y = 0, .z = level.L},
       {.x = 0.0, .y = 1.0, .z = 0.0},
       level.bAngle);
     Int3D base1 = rotate(base0,
       {.x = 0.0, .y = 0.0, .z = 1.0},
       rotateZ);
     Int3D rval(base1.x + root.x,
       base1.y + root.y,
       base1.z + root.z);
     return rval;
   }

   void constructSegment(const Level &level) {
     std::vector<Int3D> positions;
     // Build cylinder at origin along y-axis
     int64_t radius = level.d;
     double az = 0, inc = M_PI / 2;
     for (int64_t z = 0; z <= level.L; z++) {
       if (_options->lung_model_is_skeleton) {
         positions.push_back({.x = 0, .y = 0, .z = z});
       } else {
         for (az = 0; az < 2 * M_PI; az += M_PI / 180) {
           int64_t x = (int64_t) round(radius * sin(inc) * cos(az));
           int64_t y = (int64_t) round(radius * sin(inc) * sin(az));
           positions.push_back({.x = x, .y = y, .z = z});
         }
       }
     }
     // Treat as positional vectors and apply rotations and translate
     az = (level.direction.x == 0) ? 0 : atan(level.direction.y / level.direction.x);
     inc = acos(level.direction.z);
     for (int64_t i = 0; i < positions.size(); i++) {
         Int3D newPosition0 = rotate(positions.at(i),
             {.x = 0.0, .y = 1.0, .z = 0.0},
             inc);
         newPosition0 = rotate(newPosition0,
             {.x = 0.0, .y = 0.0, .z = 1.0},
             az);
         newPosition0.x += level.centroid.x;
         newPosition0.y += level.centroid.y;
         newPosition0.z += level.centroid.z;
         positions.at(i) = newPosition0;
         // Verify new location is within grid
         if ((0 <= newPosition0.x && newPosition0.x < gridSize.x) &&
             (0 <= newPosition0.y && newPosition0.y < gridSize.y) &&
             (0 <= newPosition0.z && newPosition0.z < gridSize.z)) {
           airwayEpiCellPositions1D.insert(newPosition0.x + newPosition0.y * gridSize.x +
                                     newPosition0.z * gridSize.x * gridSize.y);
         }
       }
   }

   void constructAlveoli(const Int3D & pos) {
     // Cube dimensions for dim = 5 epi cells across are x,y,z = [-2 2]
     int dim = 5, idim = 2;
     for (int x = -idim; x <= idim; x++) {
       for (int y = -idim; y <= idim; y++) {
         for (int z = 0; z < dim; z++) {
           // Cells in the two x planes
           if (x == -idim || x == idim) {
             addPosition(x, y, z, pos);
           }
           // Cells in the two y planes
           if (y == -idim || y == idim) {
             addPosition(x, y, z, pos);
           }
           // Cells in the one z plane at bottom of alveoli
           if (z == dim - 1) {
             addPosition(x, y, z, pos);
           }
         }
       }
     }
     // TODO rotate
   }

   void addPosition(int x, int y, int z, const Int3D & pos) {
     positions.push_back({.x = x + pos.x,
       .y = y + pos.y,
       .z = z + pos.z});
     alveoliEpiCellPositions1D.insert((x + pos.x)
      + (y + pos.y) * gridSize.x
      + (z + pos.z) * gridSize.x * gridSize.y);
   }

   void loadEstimatedParameters() {
     // Yeh et al 1980
     std::ifstream file;
     file.open(_options->lung_model_dir + "/table.txt");
     std::string line;
     if (file.is_open()) {
       while (std::getline(file, line)) {
         std::stringstream lstream(line);
         Level e;
         int tempI;
         double tempD;
         lstream >> tempI;
         lstream >> e.count;
         lstream >> tempD;
         e.L = (int64_t) round(scale * tempD);
         lstream >> tempD;
         e.d = (int64_t) round(scale * tempD);
         lstream >> tempI;
         e.bAngle = tempI * M_PI / 180;
         lstream >> tempI;
         e.gAngle = tempI * M_PI / 180;
         levels.push_back(e);
       }
       SLOG("Loaded ", levels.size(), " estimated levels\n");
       file.close();
     } else {
       SDIE("Failed to open file ", _options->lung_model_dir + "/table.txt");
     }
   }

   void loadEmpiricalData() {
     // Bauer et al 2019
     std::ifstream file;
     file.open(_options->lung_model_dir + "/table.csv");
     std::string line;
     if (file.is_open()) {
       while (std::getline(file, line)) {
         std::stringstream lstream(line);
         Level e;
         int tempI;
         double tempD, x, y, z;
         std::string tempS, comma;
         while(lstream.good()) {
           string str;
           getline(lstream, str, ','); // 0
           getline(lstream, str, ','); // 1
           getline(lstream, str, ','); // 2 length
           e.L = (int64_t) round(scale * std::stod(str));
           getline(lstream, str, ','); // 3 diameter
           e.d = (int64_t) round(scale * std::stod(str));
           getline(lstream, str, ','); // 4 empty?
           // 5-7 centroid
           getline(lstream, str, ',');
           e.centroid.x = (int64_t) round(scale * std::stod(str));
           getline(lstream, str, ',');
           e.centroid.y = (int64_t) round(scale * std::stod(str));
           getline(lstream, str, ',');
           e.centroid.z = (int64_t) round(scale * std::stod(str));
           // 8-10 Direction
           getline(lstream, str, ',');
           e.direction.x = std::stod(str);
           getline(lstream, str, ',');
           e.direction.y = std::stod(str);
           getline(lstream, str);
           e.direction.z = std::stod(str);
           levels.push_back(e);
         }
       }
       SLOG("Loaded ", levels.size(), " airway segments\n");
       file.close();
     } else {
       SDIE("Failed to open file ", _options->lung_model_dir + "/table.csv");
     }
   }

};
