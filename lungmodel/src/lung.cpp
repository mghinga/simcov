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
#include <atomic>
#include <omp.h>

#include "prime.hpp"

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
    int iteration, index;//, end;
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
int32_t minx = 0, maxx = 0, miny = 0, maxy = 0, minz = 0, maxz = 0;
std::vector<int64_t> epiCellPositions1D;

std::string get_cur_time() {
    std::time_t result = std::time(nullptr);
    char buffer[64];
    buffer[0] = '\0';
    size_t sz = strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S ", std::localtime(&result));
    return std::string(buffer);
}
  
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

void setLocalModelBounds(const int32_t(&pos)[3], int32_t &_minx, int32_t &_maxx, int32_t &_miny, int32_t &_maxy, int32_t &_minz, int32_t &_maxz) {
    if (pos[0] > _maxx) {
        _maxx = pos[0];
    } else if (pos[0] < _minx) {
        _minx = pos[0];
    }
    if (pos[1] > _maxy) {
        _maxy = pos[1];
    } else if (pos[1] < _miny) {
        _miny = pos[1];
    }
    if (pos[2] > _maxz) {
        _maxz = pos[2];
    } else if (pos[2] < _minz) {
        _minz = pos[2];
    }
}

void setModelBounds(int32_t _minx, int32_t _maxx, int32_t _miny, int32_t _maxy, int32_t _minz, int32_t _maxz) {
    minx = std::min(minx, _minx);
    maxx = std::max(maxx, _maxx);
    miny = std::min(miny, _miny);
    maxy = std::max(maxy, _maxy);
    minz = std::min(minz, _minz);
    maxz = std::max(maxz, _maxz);
}

int32_t addPosition(int32_t x,
    int32_t y,
    int32_t z,
    const int32_t(&pos)[3],
    double bAngle,
    double rotateZ,
    int32_t &_minx, int32_t &_maxx, int32_t &_miny, int32_t &_maxy, int32_t &_minz, int32_t &_maxz
    ) {
    // Treat as positional vectors and apply rotations and translate y
    int32_t position[3] = { x, y, z };
    rotate(position, yaxis, bAngle);
    rotate(position, zaxis, rotateZ);
    position[0] += pos[0];
    position[1] += pos[1];
    position[2] += pos[2];
    /*
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
    */
    // Compute model min max dimension boundaries
    setLocalModelBounds(position, _minx, _maxx, _miny, _maxy, _minz, _maxz);//TODO maybe reduce?
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
    int32_t _minx = 0, _maxx = 0, _miny = 0, _maxy = 0, _minz = 0, _maxz = 0;
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
                        rotateZ,
                        _minx, _maxx, _miny, _maxy, _minz, _maxz);
                    positions.push_back(newPos);
                }
            }
        }
    }
#pragma omp critical
    {
        setModelBounds(_minx, _maxx, _miny, _maxy, _minz, _maxz);
        numAlveoli++;
    }
}

void constructSegment(const int32_t (&root)[3],
    const int32_t(&newRoot)[3],
    const Level &level,
    double rotateZ,
    bool isTerminal) {
    std::vector<int64_t> positions;
    int32_t _minx = 0, _maxx = 0, _miny = 0, _maxy = 0, _minz = 0, _maxz = 0;
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
            newPos = addPosition(x, y, z, root, level.bAngle, rotateZ,
                                 _minx, _maxx, _miny, _maxy, _minz, _maxz);
            positions.push_back(newPos);
        }
    }
#pragma omp critical
    {
        setModelBounds(_minx, _maxx, _miny, _maxy, _minz, _maxz);
        numAirways++;
    }
    // Draw alveolus at each terminal airway
    if (isTerminal) {
        constructAlveoli(newRoot, level.bAngle, rotateZ, positions);
    }
#pragma omp critical
    {
        epiCellPositions1D.insert(epiCellPositions1D.end(),
                                  positions.begin(),
                                  positions.end());
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

void print(bool write_to_file, int lobe) {
    std::fprintf(stderr,
        "Airways " "%" PRId64 " Alveolus " "%" PRId64 "\n",
        numAirways,
        numAlveoli);
    if (write_to_file) {
        // sort to get same values even after parallel creation - only needed for debugging
        //std::sort(epiCellPositions1D.begin(), epiCellPositions1D.end());
        std::ofstream ofs;
        /*
          ofs.open("airway.csv", std::ofstream::out);
          if (!ofs.is_open()) {
          std::fprintf(stderr, "Could not create file");
          exit(1);
          }
          for (const int64_t& a : epiCellPositions1D) {
          ofs << a << std::endl;
          }
          std::fprintf(stderr, "Bytes written %ld\n", (long)ofs.tellp());
          ofs.close();
        */
        std::string fname = std::string("airway-") + std::to_string(lobe) + std::string(".bin");
        ofs.open(fname, std::ios::out | std::ios::binary);
        if (!ofs.is_open()) {
            std::fprintf(stderr, "Could not create file");
            exit(1);
        }
        ofs.write((char*)&epiCellPositions1D[0], epiCellPositions1D.size() * sizeof(int64_t));
        double gbs_written = (double)ofs.tellp() / (1024 * 1024 * 1024);
        ofs.close();
        std::cerr << "Wrote " << gbs_written << " GBs to " << fname << "\n";
    }
    std::fprintf(stderr,
        "%d %d %d %d %d %d\n",
        minx, maxx, miny, maxy, minz, maxz);
}

void reduce() {
    auto t = NOW();
    int MAX_PROBE = 100;
    size_t n = epiCellPositions1D.size();
    primes::Prime prime;
    prime.set(n * 1.5, true);
    size_t max_elems = prime.get();
    // choose max int64_t as empty marker
    const int64_t EMPTY_SLOT = std::numeric_limits<int64_t>::max();
    std::vector<std::atomic<int64_t>> elems(max_elems);
    std::cerr << "Calculating intersections for " << n << " epicells\n";
    std::cerr << "sizeof size_t " << sizeof(size_t) << "\n";
#pragma omp parallel for
    for (size_t i = 0; i < max_elems; i++) {
        elems[i] = EMPTY_SLOT;
    }
    size_t numIntersectCells = 0;
    size_t num_dropped_inserts = 0;
#pragma omp parallel
    {
        size_t tenth = n / (10 * omp_get_max_threads());
        size_t my_intersects = 0;
#pragma omp for
        // Merge all positions into unordered set
        for (size_t i = 0; i < n; i++) {
            if (omp_get_thread_num() == 0 && i > 0 && i % tenth == 0) {
                int my_count = i * omp_get_max_threads();
                std::cerr << get_cur_time() << " epicells processed: " << my_count
                          << " out of " << n << " (" << round(100.0 * my_count / n) << " %)\n";
            }
            // Record location and if the cell intersects another
            auto slot = std::hash<int64_t>{}(epiCellPositions1D[i]) % max_elems;
            auto start_slot = slot;
            for (int j = 0; j < MAX_PROBE; j++) {
                int64_t old_elem = EMPTY_SLOT;
                elems[slot].compare_exchange_strong(old_elem, epiCellPositions1D[i], std::memory_order_release, std::memory_order_relaxed);
                if (old_elem == EMPTY_SLOT) {
                    break;
                } else if (old_elem == epiCellPositions1D[i]) {
                    my_intersects++;
                    break;
                }
                // quadratic probe
                slot = (start_slot + (j + 1) * (j + 1)) % max_elems;
                if (j == MAX_PROBE - 1) num_dropped_inserts++;
            }
        }
#pragma omp atomic update
        numIntersectCells += my_intersects;
    }
    size_t numCells = n;
    std::chrono::duration<double> t_elapsed = NOW() - t;
    std::cerr << "Calculating intersections took " << t_elapsed.count() << " s\n"
              << "Total cells " << numCells << " num intersections " << numIntersectCells
              << " (" << 100.0 * (double)numIntersectCells / numCells << " %)\n";
    if (num_dropped_inserts) std::cerr << "dropped " << num_dropped_inserts << std::endl;
}

int main(int argc, char *argv[]) {
    if (argc != 9) {
        printf("Usage: %s <dim_x> <dim_y> <dim_z> <offset_x> <offset_y> <offset_z> <generations> <dump>\n", argv[0]);
        return -1;
    }
    // Set input parameters
    gridSize[0] = atoi(argv[1]);
    gridSize[1] = atoi(argv[2]);
    gridSize[2] = atoi(argv[3]);
    gridOffset[0] = atoi(argv[4]);
    gridOffset[1] = atoi(argv[5]);
    gridOffset[2] = atoi(argv[6]);
    int lobe0_gens = atoi(argv[7]);
    std::string write_str(argv[8]);
    bool write_to_file = (write_str == "dump");
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

    //epiCellPositions.1d.resize(50000000);
    //size_t num_epicells = 0;
    //int generations[] = { 24, 24, 26, 24, 25 };
    int generations[] = { lobe0_gens, 24, 26, 24, 25 };
    int startIndex[] = { 0, 24, 48, 74, 98 };
    int32_t base[] = { 12628, 10516, 0 }; // Base of btree at roundUp(bounds/2)

    double t_construct = 0;

    for (int i = 0; i < 1; i++) {//TODO 5; i++) {
        std::fprintf(stderr,
            "Processing lobe %d for %d generations\n",
            i,
            generations[i]);
        auto max_branches = pow(2.0, generations[i] - 1);
        auto one_tenth = max_branches / 10;
        auto start = NOW();
        auto est_max_epicells = 2500000.0 * pow(max_branches, 0.704);
        std::cerr << "Estimated max epicells " << est_max_epicells << " max memory estimate "
                  << ((2.5 * est_max_epicells * sizeof(int64_t)) / 1024 / 1024 / 1024) << " GBs\n";
        epiCellPositions1D.reserve(est_max_epicells);
#pragma omp parallel
        {
            int num_branches_processed = 0;
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
                constructSegment(base, root, lvl, 0.0, false);
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
                        int prev_tenth_branches_processed = num_branches_processed / one_tenth;
                        num_branches_processed++;
                        int tenth_branches_processed = num_branches_processed / one_tenth;
                        if (tenth_branches_processed > prev_tenth_branches_processed) {
                            std::cerr << get_cur_time() << " branches processed: " << num_branches_processed
                                      << " out of " << max_branches
                                      << " (" << (int)(100 * num_branches_processed / max_branches) << " %)\n";;
                        }
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
                                auto tc = NOW();
                                constructSegment(branch.root,
                                    rchild,
                                    lvl,
                                    rotateZ,
                                    isTerminal);
                                std::chrono::duration<double> t_elapsed = NOW() - tc;
#pragma omp atomic update
                                t_construct += t_elapsed.count();
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
                                auto tc = NOW();
                                constructSegment(branch.root,
                                    lchild,
                                    lvl,
                                    rotateZ,
                                    isTerminal);
                                std::chrono::duration<double> t_elapsed = NOW() - tc;
#pragma omp atomic update
                                t_construct += t_elapsed.count();
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
        t_construct /= omp_get_max_threads();
        std::chrono::duration<double> t_elapsed = NOW() - start;
        std::cerr << "Construction completed in " << t_elapsed.count()
                  << " s, with " << epiCellPositions1D.size() << " epicells\n";
        auto t_reduce = NOW();
        reduce();
        std::chrono::duration<double> t_reduce_elapsed = NOW() - t_reduce;
        auto t_print = NOW();
        print(write_to_file, i);
        std::chrono::duration<double> t_print_elapsed = NOW() - t_print;
        t_elapsed = NOW() - start;
        std::fprintf(stderr, "-------\nconstructSegment time %.4f s, %.2f %%\n", t_construct, 100.0 * t_construct / t_elapsed.count());
        std::fprintf(stderr, "reduce time %.4f s, %.2f %%\n", t_reduce_elapsed.count(), 100.0 * t_reduce_elapsed.count() / t_elapsed.count());
        std::fprintf(stderr, "print time %.4f s, %.2f %%\n", t_print_elapsed.count(), 100.0 * t_print_elapsed.count() / t_elapsed.count());
        std::fprintf(stderr, "Lobe %d total time %g\n\n", i, t_elapsed.count());
    }
    return 0;
}
