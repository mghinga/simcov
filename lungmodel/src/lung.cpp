
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <memory>
#include <fcntl.h>
#include <inttypes.h>

class Random {
private:
    std::mt19937_64 generator;
    std::uniform_int_distribution<int32_t> distribution;
public:
    static constexpr double deg2rad = M_PI / 180;
    Random(unsigned seed, int32_t minv, int32_t maxv)
        : generator(seed), distribution(minv, maxv) {}
    double get() {
        return distribution(generator) * deg2rad;
    }
};

struct Level {
    int32_t L = 0, r = 0, count = 0;
    double bAngle = 0.0, gAngle = 0.0;
};

// Model variables
const double twoPi = M_PI * 2;
const double cylinderRadialIncrement = M_PI / 2;
const int idim = 20;
const int idim2 = 2 * idim;
const int idim3 = idim2 - 1;
double yaxis[3] = { 0.0, 1.0, 0.0 };
double zaxis[3] = { 0.0, 0.0, 1.0 };
std::vector<Level> levels;
std::shared_ptr<Random> rnd_gen = std::make_shared<Random>(753695190, 0, 179);
time_t startProgram, endProgram, startFunction, endFunction;
double timeProgram = 0.0, timeFunction = 0.0;

// Input parameters
int32_t gridSize[3] = {0}; // Simulation space size
int32_t gridOffset[3] = {0}; // Position of simulation space in model

// Output variables
int64_t numAirways = 0, numAirwayCells = 0;
int64_t numAlveoli = 0, numAlveoliCells = 0;
int64_t numIntersectCells = 0;
int32_t minx = 0, maxx = 0, miny = 0, maxy = 0, minz = 0, maxz = 0;
std::unordered_set<int64_t> epiCellPositions1D;

void rotate(int32_t (&vec)[3], const double (&axis)[3], double angle) {
    /**
    * http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
    */
    double cosAngle = cos(angle);
    double sinAngle = sin(angle);
    double x = vec[0], y = vec[1], z = vec[2];
    // double u = 0, v = 1, w = 0;
    double u = axis[0], v = axis[1], w = axis[2];
    double dotProduct = (x * u) + (y * v) + (z * w);
    double newX = u * dotProduct * (1 - cosAngle) + x
        * cosAngle + (-w * y + v * z) * sinAngle;
    double newY = v * dotProduct * (1 - cosAngle) + y
        * cosAngle + (w * x - u * z) * sinAngle;
    double newZ = w * dotProduct * (1 - cosAngle) + z
        * cosAngle + (-v * x + u * y) * sinAngle;
    vec[0] = (int32_t)round(newX);
    vec[1] = (int32_t)round(newY);
    vec[2] = (int32_t)round(newZ);
}

void setModelBounds(const int32_t(&pos)[3]) {
    if (pos[0] > maxx) {
        maxx = pos[0];
    } else if (pos[0] < minx) {
        minx = pos[0];
    }
    if (pos[1] > maxy) {
        maxy = pos[1];
    } else if (pos[1] < miny) {
        miny = pos[1];
    }
    if (pos[2] > maxz) {
        maxz = pos[2];
    } else if (pos[2] < minz) {
        minz = pos[2];
    }
}

void addPosition(int32_t x,
    int32_t y,
    int32_t z,
    const int32_t(&pos)[3],
    double bAngle,
    double rotateZ,
    int64_t& numCells) {
    // Treat as positional vectors and apply rotations and translate y
    int32_t position[3] = { x, y, z };
    rotate(position, yaxis, bAngle);
    rotate(position, zaxis, rotateZ);
    position[0] += pos[0];
    position[1] += pos[1];
    position[2] += pos[2];
#ifdef SLM_WRITE_TO_FILE
    // Verify new location is within grid
    if (position[0] < gridOffset[0] ||//TODO optimize for production
        position[0] > (gridSize[0] + gridOffset[0]) ||
        position[1] < gridOffset[1] ||
        position[1] > (gridSize[1] + gridOffset[1]) ||
        position[2] < gridOffset[2] ||
        position[2] > (gridSize[2] + gridOffset[2])) {
        return;
    }
#endif
    // Compute model min max dimension boundaries
    setModelBounds(position);
    // Record new location and if the cell intersects another
    auto success = epiCellPositions1D.insert(position[0] + position[1]
        * gridSize[0] + position[2] * gridSize[0] * gridSize[1]);
    if (!success.second) {
        numIntersectCells++;
    }
    // Increment total cell count
    numCells++;
}

void constructAlveoli(const int32_t (&pos)[3], double bAngle, double rotateZ) {
    // Single alveolar volume 200x200x200 um, 40x40x40 units, [-20, 20]
    for (int32_t x = -idim; x <= idim; x++) {
        for (int32_t y = -idim; y <= idim; y++) {
            for (int32_t z = 0; z < idim2; z++) {
                if (x == -idim // Cells in the x planes
                    || x == idim
                    || y == -idim // Cells in the y planes
                    || y == idim
                    || z == idim3) { // Cells in z alveolar bottom
                    addPosition(x,
                        y,
                        z,
                        pos,
                        bAngle,
                        rotateZ,
                        numAlveoliCells);
                }
            }
        }
    }
    numAlveoli++;
}

void constructSegment(const int32_t (&root)[3],
    const Level &level,
    double rotateZ,
    bool isTerminal,
    int32_t (&newRoot)[3]) {
    // Build cylinder at origin along y-axis
    time(&startFunction);
    for (int32_t z = 0; z <= level.L; z++) {
        for (double az = 0; az < twoPi; az += Random::deg2rad) {
            int32_t x = (int32_t)round(level.r
                * sin(cylinderRadialIncrement)
                * cos(az));
            int32_t y = (int32_t)round(level.r
                * sin(cylinderRadialIncrement)
                * sin(az));
            addPosition(x, y, z, root, level.bAngle, rotateZ, numAirwayCells);
        }
    }
    // Create root for next generation
    int32_t base[3] = {0, 0, level.L};
    rotate(base, yaxis, level.bAngle);
    rotate(base, zaxis, rotateZ);
    newRoot[0] = root[0] + base[0];
    newRoot[1] = root[1] + base[1];
    newRoot[2] = root[2] + base[2];
    numAirways++;
    // Draw alveolus at each terminal airway
    if (isTerminal) {
        constructAlveoli(newRoot, level.bAngle, rotateZ);
    }
    time(&endFunction);
    timeFunction += double(endFunction - startFunction);
}

void construct(const int32_t (&root)[3],
    int iteration,
    int index,
    int end,
    double previousBranchAngle,
    double previousRotAngle) {
    if (iteration > end) {
        return;
    }
    bool isTerminal = (iteration == end) ? true : false;
    // Uniform randomly rotate each branch
    double rotateZ = (iteration >= 2) ? rnd_gen->get() : 0.0;
    rotateZ = previousRotAngle - rotateZ;
    // Draw left child bronchiole
    Level lvl = levels.at(index);
    lvl.bAngle = previousBranchAngle - lvl.bAngle;
    int32_t lchild[3] = { 0, 0, 0 };
    constructSegment(root, lvl, rotateZ, isTerminal, lchild);
    construct(lchild, iteration + 1, index + 1, end, lvl.bAngle, rotateZ);
    // Uniform randomly rotate each branch
    rotateZ = (iteration >= 2) ? rnd_gen->get() : 0.0;
    rotateZ = previousRotAngle + rotateZ;
    // Draw right child bronchiole
    lvl = levels.at(index);
    lvl.bAngle = previousBranchAngle + lvl.bAngle;
    lvl.gAngle = -lvl.gAngle;
    int32_t rchild[3] = { 0, 0, 0 };
    constructSegment(root, lvl, rotateZ, isTerminal, rchild);
    construct(rchild, iteration + 1, index + 1, end, lvl.bAngle, rotateZ);
}

void loadEstimatedParameters() {
    /* Yeh et al 1980
    scale = 2000  => 1 unit = 5um, 1cm = 10^-2m = 10^4 um, 10^4/5 = 2000 units
    */
    double scale = 2000; //10; // for reduced scale full lung model 300x300x300
    std::ifstream file;
    file.open("table.txt");
    std::string line;
    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::stringstream lstream(line);
            Level e;
            int32_t tempI;
            double tempD;
            lstream >> tempI;
            lstream >> e.count;
            lstream >> tempD;
            e.L = (int32_t)round(scale * tempD);
            lstream >> tempD;
            e.r = (int32_t)(round(scale * tempD) / 2);
            lstream >> tempI;
            e.bAngle = tempI * Random::deg2rad;
            lstream >> tempI;
            e.gAngle = tempI * Random::deg2rad;
            if (levels.size() > 73) {  // Negate for left lobe data
                e.bAngle *= -1.0;
            }
            levels.push_back(e);
        }
        std::fprintf(stderr, "Loaded estimated levels %ld\n", levels.size());
        file.close();
    } else {
        std::fprintf(stderr, "Failed to open table.txt\n");
    }
}

void reset() {
    minx = maxx = miny = maxy = minz = maxz = 0;
    numIntersectCells = 0;
    numAlveoli = numAlveoliCells = 0;
    numAirways = numAirwayCells = 0;
    epiCellPositions1D.clear();
}

int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("Usage: %s <dim_x> <dim_y> <dim_z> <offset_x> <offset_y> <offset_z>\n", argv[0]);
        return -1;
    }
    // Set input parameters
    gridSize[0] = atoi(argv[1]);
    gridSize[1] = atoi(argv[2]);
    gridSize[2] = atoi(argv[3]);
    gridOffset[0] = atoi(argv[4]);
    gridOffset[1] = atoi(argv[5]);
    gridOffset[2] = atoi(argv[6]);
    time(&startProgram);
    loadEstimatedParameters();
    /**
    * Starting at root and recursively build tree
    * 
    * lung lobe, num gen to model, starting row in table.txt
    * ******************************************************
    * upper right lobe, 24, 0
    * middle right, 24, 24
    * lower right, 26, 48
    * upper left, 24, 74
    * lower left, 25, 98
    */
    std::ofstream ofs;
    int generations[] = { 16 };//TODO 24, 24, 26, 24, 25};
    int startIndex[] = { 0 };//TODO 24, 48, 74, 98};
    int32_t base[] = { 12628, 10516, 0 }; // Base of btree at roundUp(bounds/2)
    for (int i = 0; i < 1; i++) {//TODO 5; i++) {
        std::fprintf(stderr,
            "Processing lobe %d generations %d\n",
            i,
            generations[i]);
        Level lvl = levels.at(startIndex[i]);
        int32_t root[3] = { 0, 0, 0 };
        constructSegment(base, lvl, 0.0, false, root);
        construct(root,
            1,
            startIndex[i] + 1,
            generations[i] - 1,
            lvl.bAngle,
            0.0);
        // Print stats
        std::fprintf(stderr,
            "Alveolars " "%" PRId64 " cells " "%" PRId64 "\n",
            numAlveoli,
            numAlveoliCells);
        std::fprintf(stderr,
            "Airways " "%" PRId64 " cells " "%" PRId64 "\n",
            numAirways,
            numAirwayCells);
        std::fprintf(stderr,
            "Total cells " "%" PRId64 "\n",
            numAlveoliCells + numAirwayCells);
#ifdef SLM_WRITE_TO_FILE
        ofs.open("airway.csv", std::ofstream::out | std::ofstream::app);
        if (!ofs.is_open()) {
            std::fprintf(stderr, "Could not create file");
            exit(1);
        }
        for (const int64_t& a : epiCellPositions1D) {
            ofs << a << std::endl;
        }
        std::fprintf(stderr, "Bytes written %ld\n", (long)ofs.tellp());
        ofs.close();
#endif
        std::fprintf(stderr,
            "%d %d %d %d %d %d\n",
            minx, maxx, miny, maxy, minz, maxz);
        std::fprintf(stderr,
            "Cells intersecting " "%" PRId64 "\n",
            numIntersectCells);
        std::fprintf(stderr,
            "Recorded cells " "%" PRId64 "\n",
            epiCellPositions1D.size());
        reset();
    }
    time(&endProgram);
    // Print timing
    timeProgram = double(endProgram - startProgram);
    std::fprintf(stderr,
        "Total time %g Time in function %g\n",
        timeProgram,
        timeFunction);
    return 0;
}
