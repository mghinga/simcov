#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

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

 Lung() {
   std::cout <<"Data folder: " << dataFolder << std::endl;
   std::cout << "Executable working directory: ";
   char pBuf[256];
   ssize_t len = sizeof(pBuf);
   int bytes = min(readlink("/proc/self/exe", pBuf, len), len - 1);
   if (bytes >= 0) {
     pBuf[bytes] = '\0';
     std::cout << pBuf << std::endl;
   } else {
     std::cout << "FAILED" << std::endl;
   }
 }

 ~Lung() {}

 const std::set<int64_t>& getEpiCellIds() {
   return epiCellPositions1D;
 }

 const std::vector<Int3D>& getEpiLocations() {
   return epiCellPositions3D;
 }

 const std::vector<Level>& getLevels() {
   return levels;
 }

 int getSkippedAirwayCount() {
   return skipped;
 }

 void init(int64_t grid_x, int64_t grid_y, int64_t grid_z) {
   gridSize.x = grid_x;
   gridSize.y = grid_y;
   gridSize.z = grid_z;
   switch (model) {
     case 0: {
       loadEstimatedParameters();
       // Draw root
       Level lvl = levels.at(startIndex);
       Int3D child = buildSegment({gridSize.x/2 - lvl.d, gridSize.y/2, 0}, lvl);
       // Recursively build tree
       build(child, (startIndex + 1), lvl.bAngle);
       break;
     }
     case 1: {
       loadEmpiricalData();
       // Draw all segments
       for (int i = 0; i < airwaysLimit && i < levels.size(); i++) {
         if (levels.at(i).L < 1 || levels.at(i).d < 1) {
           skipped++;
           continue;
         }
         buildSegment(levels.at(i));
       }
       break;
     }
     default: break;
   }
 }

 private:

   string dataFolder = "/mnt/c/Users/3n571w2/Downloads/simcov/data/";
   bool isSkeleton = true;//false;//
   int model = 0;//1;//
   int startIndex = 0;
   int maxDepth = 0; // Yeh et al 1980
   int airwaysLimit = 9999, skipped = 0; // Bauer et al 2019
   //TODO std::vector<Int3D> leafs;
   std::vector<Level> levels;
   std::vector<Int3D> epiCellPositions3D;
   std::set<int64_t> epiCellPositions1D;
   Int3D gridSize;

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

   void build(const Int3D & root,
     int generation,
     double previousBranchAngle) {
       if (generation > (maxDepth - 1)) {
         //TODO leafs.add(root);
         return;
       }
       // Draw left child bronchiole
       Level lvl = levels.at(generation);
       lvl.bAngle = previousBranchAngle - lvl.bAngle;
       Int3D child = buildSegment(root, lvl);
       build(child, generation + 1, lvl.bAngle);
       // Draw right child bronchiole
       lvl = levels.at(generation);
       lvl.bAngle = previousBranchAngle + lvl.bAngle;
       if (lvl.count > 1) {
         lvl.gAngle = -lvl.gAngle;
         child = buildSegment(root, lvl);
         build(child, generation + 1, lvl.bAngle);
       }
     }

   Int3D buildSegment(const Int3D &root, const Level &level) {
     std::vector<Int3D> epiCellPositions3D;
     // Build cylinder at origin along y-axis
     int64_t radius = level.d / 2;
     double az = 0;
     double inc = M_PI / 2;
     for (int64_t z = 0; z <= level.L; z++) {
       if (isSkeleton) {
         epiCellPositions3D.push_back({.x = 0, .y = 0, .z = z});
       } else {
         for (az = 0; az < 2 * M_PI; az += M_PI / 180) {
           int64_t x = (int64_t) round(radius * sin(inc) * cos(az));
           int64_t y = (int64_t) round(radius * sin(inc) * sin(az));
           epiCellPositions3D.push_back({.x = x, .y = y, .z = z});
         }
       }
     }
     // Treat as positional vectors and apply rotations and translate
     for (int64_t i = 0; i < epiCellPositions3D.size(); i++) {
       Int3D newPosition0 = rotate(epiCellPositions3D.at(i),
        {.x = 0.0, .y = 1.0, .z = 0.0},
        level.bAngle);
       newPosition0.x += root.x;
       newPosition0.y += root.y;
       newPosition0.z += root.z;
       epiCellPositions3D.at(i) = newPosition0;
       // Verify new location is within grid
       if ((0 <= newPosition0.x && newPosition0.x < gridSize.x)
        && (0 <= newPosition0.y && newPosition0.y < gridSize.y)
        && (0 <= newPosition0.z && newPosition0.z < gridSize.z)) {
          epiCellPositions1D.insert(newPosition0.x
            + newPosition0.y * gridSize.x
            + newPosition0.z * gridSize.x * gridSize.y);
       }
     }
     // Return root for next generation
     Int3D base = rotate({.x = 0, .y = 0, .z = level.L},
       {.x = 0.0, .y = 1.0, .z = 0.0},
       level.bAngle);
     Int3D rval(base.x + root.x,
       base.y + root.y,
       base.z + root.z);
     return rval;
   }

   void buildSegment(const Level &level) {
     std::vector<Int3D> epiCellPositions3D;
     // Build cylinder at origin along y-axis
     int64_t radius = level.d;
     double az = 0, inc = M_PI / 2;
     for (int64_t z = 0; z <= level.L; z++) {
       if (isSkeleton) {
         epiCellPositions3D.push_back({.x = 0, .y = 0, .z = z});
       } else {
         for (az = 0; az < 2 * M_PI; az += M_PI / 180) {
           int64_t x = (int64_t) round(radius * sin(inc) * cos(az));
           int64_t y = (int64_t) round(radius * sin(inc) * sin(az));
           epiCellPositions3D.push_back({.x = x, .y = y, .z = z});
         }
       }
     }
     // Treat as positional vectors and apply rotations and translate
     az = (level.direction.x == 0) ? 0 : atan(level.direction.y
         / level.direction.x);
     inc = acos(level.direction.z);
     for (int64_t i = 0; i < epiCellPositions3D.size(); i++) {
         Int3D newPosition0 = rotate(epiCellPositions3D.at(i),
             {.x = 0.0, .y = 1.0, .z = 0.0},
             inc);
         newPosition0 = rotate(newPosition0,
             {.x = 0.0, .y = 0.0, .z = 1.0},
             az);
         newPosition0.x += level.centroid.x;
         newPosition0.y += level.centroid.y;
         newPosition0.z += level.centroid.z;
         epiCellPositions3D.at(i) = newPosition0;
         // Verify new location is within grid
         if ((0 <= newPosition0.x && newPosition0.x < gridSize.x)
          && (0 <= newPosition0.y && newPosition0.y < gridSize.y)
          && (0 <= newPosition0.z && newPosition0.z < gridSize.z)) {
            epiCellPositions1D.insert(newPosition0.x
              + newPosition0.y * gridSize.x
              + newPosition0.z * gridSize.x * gridSize.y);
         }
       }
   }

   void loadEstimatedParameters() {
     // Yeh et al 1980
     std::ifstream file;
     file.open(dataFolder + "table.txt");
     maxDepth = 0;
     std::string line;
     double m = 10; // Convert cm to mm
     if (file.is_open()) {
       while (std::getline(file, line)) {
         std::stringstream lstream(line);
         Level e;
         int tempI;
         double tempD;
         lstream >> tempI;
         lstream >> e.count;
         lstream >> tempD;
         e.L = (int64_t) round(m * tempD);
         lstream >> tempD;
         e.d = (int64_t) round(m * tempD);
         lstream >> tempI;
         e.bAngle = tempI * M_PI / 180;
         lstream >> tempI;
         e.gAngle = tempI * M_PI / 180;
         levels.push_back(e);
         maxDepth++;
       }
       std::cout << "Loaded " << levels.size() <<" estimated levels"<<std::endl;
       file.close();
     } else {
       std::cout << "Failed open file: " << dataFolder+"table.txt" << std::endl;
       std::exit(1);
     }
   }

   void loadEmpiricalData() {
     // Bauer et al 2019
     std::ifstream file;
     file.open(dataFolder + "table.csv");
     std::string line;
     double m = 10; // Convert 10^0 mm to 10^1 mm
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
           e.L = (int64_t) round(m * std::stod(str));
           getline(lstream, str, ','); // 3 diameter
           e.d = (int64_t) round(m * std::stod(str));
           getline(lstream, str, ','); // 4 empty?
           // 5-7 centroid
           getline(lstream, str, ',');
           e.centroid.x = (int64_t) round(m * std::stod(str));
           getline(lstream, str, ',');
           e.centroid.y = (int64_t) round(m * std::stod(str));
           getline(lstream, str, ',');
           e.centroid.z = (int64_t) round(m * std::stod(str));
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
       std::cout << "Loaded " << levels.size() << " airway segments"<<std::endl;
       file.close();
     } else {
       std::cout << "Failed open file: " << dataFolder+"table.csv" << std::endl;
       std::exit(1);
     }
   }

};
