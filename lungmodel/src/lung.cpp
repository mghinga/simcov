// LungModel
//
// Akil Andrews, UNM May 2020

#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <vector>
#include <stack>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <memory>
#include <fcntl.h>
#include <inttypes.h>
#include <chrono>
#include <omp.h>

#define NOW std::chrono::high_resolution_clock::now

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

struct Branch {
    int32_t root[3];
    int iteration, index, end;
    double previousBranchAngle, previousRotAngle;
    Branch(int32_t *root,
        int iteration,
        int index,
        int previousBranchAngle,
        int previousRotAngle) :
        iteration(iteration),
        index(index),
        previousBranchAngle(previousBranchAngle),
        previousRotAngle(previousRotAngle) {
        std::copy(root, root + 3, this->root);
    }
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
std::stack<Branch> branches;
std::shared_ptr<Random> rnd_gen = std::make_shared<Random>(753695190, 0, 179);

// Input parameters
int32_t gridSize[3] = {0}; // Simulation space size
int32_t gridOffset[3] = {0}; // Position of simulation space in model

// Output variables
int64_t numAirways = 0, numAlveoli = 0;
int64_t numCells = 0, numIntersectCells = 0;
int32_t minx = 0, maxx = 0, miny = 0, maxy = 0, minz = 0, maxz = 0;
std::unordered_set<int64_t> results;
std::vector<int64_t> epiCellPositions1D;

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

int32_t addPosition(int32_t x,
    int32_t y,
    int32_t z,
    const int32_t(&pos)[3],
    double bAngle,
    double rotateZ) {
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
#pragma omp critical
    {
        // Compute model min max dimension boundaries
        setModelBounds(position);//TODO maybe reduce?
    }
    // Return new location in 1D
    return (position[0] + position[1]
        * gridSize[0] + position[2] * gridSize[0] * gridSize[1]);
}

void constructAlveoli(const int32_t (&pos)[3],
    double bAngle,
    double rotateZ,
    std::vector<int64_t>& positions) {
    // Single alveolar volume 200x200x200 um, 40x40x40 units, [-20, 20]
    int64_t newPos;
    for (int32_t x = -idim; x <= idim; x++) {
        for (int32_t y = -idim; y <= idim; y++) {
            for (int32_t z = 0; z < idim2; z++) {
                if (x == -idim // Cells in the x planes
                    || x == idim
                    || y == -idim // Cells in the y planes
                    || y == idim
                    || z == idim3) { // Cells in z alveolar bottom
                    newPos = addPosition(x,
                        y,
                        z,
                        pos,
                        bAngle,
                        rotateZ);
                    positions.push_back(newPos);
                }
            }
        }
    }
#pragma omp critical
    {
        numAlveoli++;
    }
}

void constructSegment(const int32_t (&root)[3],
    const int32_t(&newRoot)[3],
    const Level &level,
    double rotateZ,
    bool isTerminal,
    std::vector<int64_t> &positions) {
    // Build cylinder at origin along y-axis
    int32_t x, y, newPos;
    for (int32_t z = 0; z <= level.L; z++) {
        for (double az = 0; az < twoPi; az += Random::deg2rad) {
            x = (int32_t)round(level.r
                * sin(cylinderRadialIncrement)
                * cos(az));
            y = (int32_t)round(level.r
                * sin(cylinderRadialIncrement)
                * sin(az));
            newPos = addPosition(x, y, z, root, level.bAngle, rotateZ);
            positions.push_back(newPos);
        }
    }
#pragma omp critical
    {
        numAirways++;
    }
    // Draw alveolus at each terminal airway
    if (isTerminal) {
        constructAlveoli(newRoot, level.bAngle, rotateZ, positions);
    }
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
    /*TODO
    Computer the average number of alveoli per sac
    alveolusPerSac
    24, 0
    24, 24
    26, 48
    24, 74
    25, 98
    -- , 123
     */

}

void print() {
    std::fprintf(stderr,
        "Airways " "%" PRId64 " Alveolus " "%" PRId64 "\n",
        numAirways,
        numAlveoli);
#ifdef SLM_WRITE_TO_FILE
    std::ofstream ofs;
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
        "Cells intersecting " "%" PRId64 " total cells " "%" PRId64 "\n",
        numIntersectCells,
        numCells);
}

void reduce() {
    // Merge all positions into unordered set
    for (int i = 0; i < epiCellPositions1D.size(); i++) {
        // Record location and if the cell intersects another
        auto success = results.insert(epiCellPositions1D[i]);
        if (!success.second) {
            numIntersectCells++;
        }
        numCells++;
    }
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
    loadEstimatedParameters();
    /**
    * Starting at root and preorder iteratively build tree
    * 
    * lung lobe, num gen to model, starting row in table.txt
    * ******************************************************
    * upper right lobe, 24, 0
    * middle right, 24, 24
    * lower right, 26, 48
    * upper left, 24, 74
    * lower left, 25, 98
    * 
    * Note* last generations are alveolus
    */
    int generations[] = { 24, 24, 26, 24, 25 };
    int startIndex[] = { 0, 24, 48, 74, 98 };
    int32_t base[] = { 12628, 10516, 0 }; // Base of btree at roundUp(bounds/2)
    for (int i = 0; i < 1; i++) {//TODO 5; i++) {
        std::fprintf(stderr,
            "Processing lobe %d generations %d\n",
            i,
            generations[i]);
        auto start = NOW();
#pragma omp parallel
        {
#pragma omp single
            {
                // Create root for next generation
                Level lvl = levels.at(startIndex[i]);
                rotate(base, yaxis, lvl.bAngle);
                rotate(base, zaxis, 0.0);
                int32_t root[3] = { 0, 0, 0 };
                root[0] = root[0] + base[0];
                root[1] = root[1] + base[1];
                root[2] = root[2] + base[2];
                {
                    std::vector<int64_t> positions;
                    constructSegment(base, root, lvl, 0.0, false, positions);
                    epiCellPositions1D.insert(epiCellPositions1D.end(),
                        positions.begin(),
                        positions.end());
                }
                branches.push(Branch(root,
                    1,
                    startIndex[i] + 1,
                    lvl.bAngle,
                    0.0));
                int lastGeneration = generations[i] - 2; // -1 for array indexing, -1 since last gen=alveolus
#pragma omp task untied // Thread switching for starvation
                {
                    while (!branches.empty()) {
                        Branch branch = branches.top();
                        branches.pop();
                        if (branch.iteration <= lastGeneration) {
                            // Determine if this is a terminal bronchiole
                            bool isTerminal = (branch.iteration == lastGeneration) ? true : false;
                            // Uniform randomly rotate branch
                            double rotateZ = (branch.iteration >= 2) ? rnd_gen->get() : 0.0;
                            rotateZ = branch.previousRotAngle + rotateZ;
                            // Draw right child bronchiole
                            Level lvl = levels.at(branch.index);
                            lvl.bAngle = branch.previousBranchAngle + lvl.bAngle;
                            lvl.gAngle = -lvl.gAngle;
                            // Create root for next generation
                            int32_t base[3] = { 0, 0, lvl.L };
                            rotate(base, yaxis, lvl.bAngle);
                            rotate(base, zaxis, rotateZ);
                            int32_t rchild[3] = { 0, 0, 0 };
                            rchild[0] = branch.root[0] + base[0];
                            rchild[1] = branch.root[1] + base[1];
                            rchild[2] = branch.root[2] + base[2];
#pragma omp task
                            {
                                std::vector<int64_t> positions;
                                constructSegment(branch.root,
                                    rchild,
                                    lvl,
                                    rotateZ,
                                    isTerminal,
                                    positions);
#pragma omp critical
                                {
                                    epiCellPositions1D.insert(epiCellPositions1D.end(),
                                        positions.begin(),
                                        positions.end());
                                }
                            }
                            // Push right child to stack first for preorder
                            branches.push(Branch(rchild,
                                branch.iteration + 1,
                                branch.index + 1,
                                lvl.bAngle,
                                rotateZ));
                            // Uniform randomly rotate branch
                            rotateZ = (branch.iteration >= 2) ? rnd_gen->get() : 0.0;
                            rotateZ = branch.previousRotAngle - rotateZ;
                            // Draw left child bronchiole
                            lvl = levels.at(branch.index);
                            lvl.bAngle = branch.previousBranchAngle - lvl.bAngle;
                            // Create root for next generation
                            rotate(base, yaxis, lvl.bAngle);
                            rotate(base, zaxis, rotateZ);
                            int32_t lchild[3] = { 0, 0, 0 };
                            lchild[0] = branch.root[0] + base[0];
                            lchild[1] = branch.root[1] + base[1];
                            lchild[2] = branch.root[2] + base[2];
#pragma omp task
                            {
                                std::vector<int64_t> positions;
                                constructSegment(branch.root,
                                    lchild,
                                    lvl,
                                    rotateZ,
                                    isTerminal,
                                    positions);
#pragma omp critical
                                {
                                    epiCellPositions1D.insert(epiCellPositions1D.end(),
                                        positions.begin(),
                                        positions.end());
                                }
                            }
                            // Push left child to stack
                            branches.push(Branch(lchild,
                                branch.iteration + 1,
                                branch.index + 1,
                                lvl.bAngle,
                                rotateZ));
                        }
                    } // end while
                } // end main task
#pragma omp taskwait
            } // end single
        } // end parallel
        std::chrono::duration<double> t_elapsed = NOW() - start;
        std::fprintf(stderr, "%d total time %g\n", i, t_elapsed.count());
        reduce();
        print();
    }
    return 0;
}
