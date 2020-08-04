// SimCov
//
// Steven Hofmeyr, LBNL May 2020

#include <fcntl.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <upcxx/upcxx.hpp>

using namespace std;

#include "options.hpp"
#include "tissue.hpp"
#include "upcxx_utils.hpp"
#include "utils.hpp"

using namespace upcxx;
using namespace upcxx_utils;

#define NOW chrono::high_resolution_clock::now

struct SimStats {
  int64_t num_infected = 0;
  int64_t num_dead_epicells = 0;
  int64_t num_tcells = 0;
  int64_t num_tcells_extravasated = 0;
  int64_t tot_num_viral_kills = 0;

  void clear() {
    num_infected = 0;
    num_dead_epicells = 0;
    num_tcells = 0;
    num_tcells_extravasated = 0;
  }

  int64_t get_all_num_infected() {
    return upcxx::reduce_one(num_infected, upcxx::op_fast_add, 0).wait();
  }

  int64_t get_all_num_tcells() {
    return upcxx::reduce_one(num_tcells, upcxx::op_fast_add, 0).wait();
  }

  int64_t get_all_num_tcells_extravasated() {
    return upcxx::reduce_one(num_tcells_extravasated, upcxx::op_fast_add, 0).wait();
  }

  int64_t get_all_num_dead_epicells() {
    return upcxx::reduce_one(num_dead_epicells, upcxx::op_fast_add, 0).wait();
  }

  int64_t get_all_num_viral_kills() {
    return upcxx::reduce_one(tot_num_viral_kills, upcxx::op_fast_add, 0).wait();
  }

};

ofstream _logstream;
bool _verbose = false;
shared_ptr<Options> _options;
shared_ptr<Random> _rnd_gen;
SimStats _sim_stats;

IntermittentTimer _generate_tcell_timer(__FILENAME__ + string(":") + "generate_tcells");
IntermittentTimer _update_tcell_timer(__FILENAME__ + string(":") + "update_tcell");
IntermittentTimer _update_virus_timer(__FILENAME__ + string(":") + "update_virus");
IntermittentTimer _update_cytokines_timer(__FILENAME__ + string(":") + "update_cytokines");
IntermittentTimer _sample_timer(__FILENAME__ + string(":") + "sample");

void initial_infection(Tissue &tissue) {
  BarrierTimer timer(__FILEFUNC__);
  // each rank generates a block of infected cells
  // for large dimensions, the chance of repeat sampling for a few points is very small
  int local_num_infections = ceil((double)_options->num_infections / rank_n());
  int64_t min_i = rank_me() * local_num_infections;
  int64_t max_i = (rank_me() + 1) * local_num_infections;
  if (max_i > _options->num_infections) max_i = _options->num_infections;
  if (min_i > _options->num_infections) min_i = _options->num_infections;
  int64_t num_infected = 0;
  barrier();
  DBG("Infecting ", (max_i - min_i), " cells from ", min_i, " to ", max_i, "\n");
  ProgressBar progbar(max_i - min_i, "Setting initial infections");
  for (int i = min_i; i < max_i; i++) {
    progbar.update();
    GridCoords coords(_rnd_gen, Tissue::grid_size);
    DBG("infection: ", coords.str() + "\n");
    tissue.inc_incoming_virus(coords, 1.0);
    num_infected++;
    upcxx::progress();
  }
  progbar.done();
  barrier();
  SLOG("Initially infected ", reduce_one(num_infected, op_fast_add, 0).wait(), " epicells\n");
}

void generate_tcells(Tissue &tissue) {
  _generate_tcell_timer.start();
  // each rank could generate some t-cells
  int local_num_tcells = _options->tcell_generation / rank_n();
  int remaining_tcells = _options->tcell_generation - local_num_tcells * rank_n();
  if (rank_me() < remaining_tcells) local_num_tcells++;
  if (local_num_tcells) {
    for (int i = 0; i < local_num_tcells; i++) {
      GridCoords coords(_rnd_gen, Tissue::grid_size);
      string tcell_id = to_string(rank_me()) + "-" + to_string(tissue.tcells_generated);
      tissue.tcells_generated++;
      DBG("tcell ", tcell_id, " starting at ", coords.str(), "\n");
      tissue.add_tcell(coords, {.id = tcell_id, .in_vasculature = true});
      upcxx::progress();
    }
  }
  _generate_tcell_timer.stop();
  barrier();
  SLOG_VERBOSE("Generated ", reduce_one(local_num_tcells, op_fast_add, 0).wait(), " t-cells\n");
}

int64_t get_rnd_coord(int64_t x, int64_t max_x) {
  int64_t new_x = x + _rnd_gen->get(0, 3) - 1;
  if (new_x < 0) new_x = 0;
  if (new_x >= max_x) new_x = max_x - 1;
  return new_x;
}

void update_tcell(int time_step, Tissue &tissue, GridPoint *grid_point, TCell &tcell) {
  _update_tcell_timer.start();
  if (tcell.in_vasculature && grid_point->icytokine > 0) tcell.in_vasculature = false;
  if (tcell.in_vasculature) {
    auto rnd_nb_i = _rnd_gen->get(0, (int64_t)grid_point->neighbors.size());
    GridCoords selected_coords = grid_point->neighbors[rnd_nb_i];
    tissue.add_tcell(selected_coords, tcell);
  } else {
    if (grid_point->virus > 0) {
      // induce apoptosis
      DBG(time_step, ": tcell ", tcell.id, " at ", grid_point->coords.str(),
          " induced apoptosis\n");
      if (grid_point->epicell->status != EpiCellStatus::APOPTOTIC) {
        grid_point->epicell->status = EpiCellStatus::APOPTOTIC;
        grid_point->epicell->start_apoptosis = time_step;
      }
      // need to add here to ensure it gets to the next time step in the vector swap
      tissue.add_tcell(grid_point->coords, tcell);
    } else {
      GridCoords selected_coords;
      double highest_chemokine = 0;
      for (auto &nb_coords : grid_point->neighbors) {
        if (nb_coords == grid_point->coords) continue;
        double chemokine = tissue.get_chemokine(nb_coords);
        if (chemokine > highest_chemokine) {
          highest_chemokine = chemokine;
          selected_coords = nb_coords;
        }
      }
      if (highest_chemokine == 0) {
        auto rnd_nb_i = _rnd_gen->get(0, (int64_t)grid_point->neighbors.size());
        selected_coords = grid_point->neighbors[rnd_nb_i];
      } else {
        DBG(time_step, ": highest nb chemokine for tcell ", tcell.id, " at ",
            grid_point->coords.str(), " is at ", selected_coords.str(), " with ", highest_chemokine,
            "\n");
      }
      DBG(time_step, ": tcell ", tcell.id, " at ", grid_point->coords.str(), " moving to ",
          selected_coords.str(), "\n");
      tissue.add_tcell(selected_coords, tcell);
    }
  }
  _update_tcell_timer.stop();
}

void update_virus(int time_step, Tissue &tissue, GridPoint *grid_point) {
  if (grid_point->virus == 0) return;
  _update_virus_timer.start();
  grid_point->epicell->num_steps_infected++;
  switch (grid_point->epicell->status) {
    case EpiCellStatus::INCUBATING:
      if (grid_point->epicell->num_steps_infected > _options->incubation_period) {
        grid_point->epicell->status = EpiCellStatus::EXPRESSING;
        grid_point->virus = 1.0;
        grid_point->chemokine = 1.0;
        grid_point->icytokine = 1.0;
      }
      break;
    case EpiCellStatus::APOPTOTIC:
      if (time_step - grid_point->epicell->start_apoptosis > _options->apoptosis_period) {
        grid_point->epicell->status = EpiCellStatus::DEAD;
        grid_point->virus = 0;
        _sim_stats.tot_num_viral_kills++;
        break;
      }
      // otherwise it's the same as expressing, so no break here
    case EpiCellStatus::EXPRESSING:
      if (grid_point->epicell->num_steps_infected >
          _options->incubation_period + _options->expressing_period) {
        grid_point->epicell->status = EpiCellStatus::DEAD;
        grid_point->virus = 0;
        break;
      }
      for (auto &nb_coords : grid_point->neighbors) {
        if (nb_coords == grid_point->coords) continue;
        tissue.inc_incoming_virus(nb_coords, grid_point->virus);
      }
      break;
    default: DIE("Invalid state ", (int)grid_point->epicell->status, " for infected cell");
  }
  _update_virus_timer.stop();
}

void update_cytokines(int time_step, Tissue &tissue, GridPoint *grid_point) {
  _update_cytokines_timer.start();
  if (grid_point->chemokine > 0) {
    if (grid_point->epicell->status != EpiCellStatus::EXPRESSING) {
      grid_point->chemokine -= _options->chemokine_decay_rate;
      if (grid_point->chemokine < 0) grid_point->chemokine = 0;
    }
    if (grid_point->chemokine > 0) {
      double chemokine_to_nb = grid_point->chemokine * _options->chemokine_diffusion_rate /
                               grid_point->neighbors.size();
      for (auto &nb_coords : grid_point->neighbors) {
        if (nb_coords == grid_point->coords) continue;
        tissue.inc_incoming_chemokines(nb_coords, chemokine_to_nb);
      }
    }
  }
  if (grid_point->icytokine > 0) {
    if (grid_point->epicell->status != EpiCellStatus::EXPRESSING) {
      grid_point->icytokine -= _options->icytokine_decay_rate;
      if (grid_point->icytokine < 0) grid_point->icytokine = 0;
    }
    if (grid_point->incoming_icytokine > 0) {
      double icytokine_to_nb = grid_point->incoming_virus * _options->icytokine_diffusion_rate /
                               grid_point->neighbors.size();
      for (auto &nb_coords : grid_point->neighbors) {
        if (nb_coords == grid_point->coords) continue;
        tissue.inc_incoming_icytokines(nb_coords, icytokine_to_nb);
      }
    }
  }
  _update_cytokines_timer.stop();
}

void finish_round(Tissue &tissue) {
  barrier();
  _sim_stats.clear();
  for (auto grid_point = tissue.get_first_local_grid_point(); grid_point;
       grid_point = tissue.get_next_local_grid_point()) {
    if (grid_point->tcells) {
      grid_point->switch_tcells_vector();
      _sim_stats.num_tcells += grid_point->tcells->size();
    }
    if (grid_point->incoming_virus > 0 && grid_point->epicell->status == EpiCellStatus::HEALTHY) {
      double infection_prob = grid_point->incoming_virus * _options->virus_infection_prob;
      if (_rnd_gen->get_prob() >= infection_prob) {
        grid_point->incoming_virus = 0;
      } else {
        grid_point->virus = grid_point->incoming_virus;
        grid_point->incoming_virus = 0;
        grid_point->epicell->num_steps_infected = 0;
        grid_point->epicell->status = EpiCellStatus::INCUBATING;
      }
    }
    if (grid_point->incoming_chemokine) {
      grid_point->chemokine = grid_point->incoming_chemokine;
      grid_point->incoming_chemokine = 0;
    }
    if (grid_point->incoming_icytokine) {
      grid_point->icytokine = grid_point->incoming_icytokine;
      grid_point->incoming_icytokine = 0;
    }
    switch (grid_point->epicell->status) {
      case EpiCellStatus::INCUBATING:
      case EpiCellStatus::EXPRESSING:
      case EpiCellStatus::APOPTOTIC: _sim_stats.num_infected++; break;
      case EpiCellStatus::DEAD: _sim_stats.num_dead_epicells++; break;
      case EpiCellStatus::HEALTHY: break;
    }
  }
  barrier();
}

void sample(int time_step, Tissue &tissue, ViewObject view_object) {
  // each rank writes its own blocks into the file at the appropriate locations. We compute the
  // location knowing that each position takes exactly 4 characters (up to 3 numbers and a space),
  // so we can work out how much space a block takes up, in order to position the data correctly.
  // Each point is represented by a scalar char. From 1-127 are used to indicate viral load, with a
  // color of red, going from faint to strong. From -127 to -1 are used to indicate tcell count,
  // with a color of blue, going from faint to strong. The last color (0) is used to indicate a dead
  // epicell (how do we make this a separate color?)
  _sample_timer.start();
  // each grid point takes up a single char
  size_t tot_sz = tissue.get_num_grid_points() * 1;
  string fname = "samples/sample_" + view_object_str(view_object) + "_" + to_string(time_step) +
                 ".vtk";
  int x_dim = _options->dimensions[0];
  int y_dim = _options->dimensions[1];
  int z_dim = _options->dimensions[2];
  ostringstream header_oss;
  header_oss << "# vtk DataFile Version 4.2\n"
             << "SimCov sample " << basename(_options->output_dir.c_str()) << time_step
             << "\n"
             //<< "ASCII\n"
             << "BINARY\n"
             << "DATASET STRUCTURED_POINTS\n"
             // we add one in each dimension because these are for drawing the
             // visualization points, and our visualization entities are cells
             << "DIMENSIONS " << (x_dim + 1) << " " << (y_dim + 1) << " " << (z_dim + 1) << "\n"
             << "SPACING 1 1 1\n"
             << "ORIGIN 0 0 0\n"
             << "CELL_DATA " << (x_dim * y_dim * z_dim) << "\n";
  switch (view_object) {
    case ViewObject::VIRUS: header_oss << "SCALARS virus unsigned_char 1\n"; break;
    case ViewObject::TCELL_VAS: header_oss << "SCALARS t-cell-vas unsigned_char 1\n"; break;
    case ViewObject::TCELL_XVAS: header_oss << "SCALARS t-cell-exvas unsigned_char 1\n"; break;
    case ViewObject::EPICELL: header_oss << "SCALARS epicell unsigned_char 1\n"; break;
    case ViewObject::ICYTOKINE: header_oss << "SCALARS icytokine unsigned_char 1\n"; break;
    case ViewObject::CHEMOKINE: header_oss << "SCALARS chemokine unsigned_char 1\n"; break;
    default: SDIE("unknown view object");
  }
  header_oss << "LOOKUP_TABLE default\n";
  if (!rank_me()) {
    // all ranks have the same size since we zero pad ASCII chars to ensure we
    // always write out 3 chars per point
    tot_sz += header_oss.str().size();
    // rank 0 creates the file and truncates it to the correct length
    auto fileno = open(fname.c_str(), O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (fileno == -1) SDIE("Cannot create file ", fname, ": ", strerror(errno), "\n");
    if (ftruncate(fileno, tot_sz) == -1)
      DIE("Could not truncate ", fname, " to ", tot_sz, " bytes\n");
    close(fileno);
    DBG("Truncated sample file ", fname, " to ", tot_sz, " bytes\n");
  }
  upcxx::barrier();
  // wait until rank 0 has finished setting up the file
  auto [bytes_written, grid_points_written] = tissue.dump_blocks(fname, header_oss.str(),
                                                                 view_object);
  upcxx::barrier();
  auto tot_bytes_written = reduce_one(bytes_written, op_fast_add, 0).wait();
  auto tot_grid_points_written = reduce_one(grid_points_written, op_fast_add, 0).wait();
  if (!rank_me()) assert(tot_grid_points_written == tissue.get_num_grid_points());
  SLOG_VERBOSE("Successfully wrote ", tot_grid_points_written, " grid points, with size ",
               get_size_str(tot_bytes_written), " to ", fname, "\n");
  _sample_timer.stop();
}

void run_sim(Tissue &tissue) {
  BarrierTimer timer(__FILEFUNC__);
  double time_step_ticks = (double)_options->num_timesteps / 20;
  auto start_t = NOW();
  auto curr_t = start_t;
  for (int time_step = 0; time_step < _options->num_timesteps; time_step++) {
    generate_tcells(tissue);
    // iterate through all local grid points
    // FIXME: this should be iteration through the local active grid points only
    for (auto grid_point = tissue.get_first_local_grid_point(); grid_point;
         grid_point = tissue.get_next_local_grid_point()) {
      upcxx::progress();
      // the tcells are moved (added to the new list, but only cleared out at the end of all
      // updates)
      if (grid_point->tcells) {
        for (auto &tcell : *grid_point->tcells) {
          update_tcell(time_step, tissue, grid_point, tcell);
        }
      }
      update_virus(time_step, tissue, grid_point);
      update_cytokines(time_step, tissue, grid_point);
    }
    finish_round(tissue);
    // print every 5% of the iterations
    if (time_step >= time_step_ticks || time_step == _options->num_timesteps - 1) {
      chrono::duration<double> t_elapsed = NOW() - curr_t;
      curr_t = NOW();
      // memory doesn't really change so don't report it every iteration
      SLOG(setw(5), left, time_step, " [", get_current_time(true), " ", setprecision(2), fixed,
           setw(5), right, t_elapsed.count(), "s", "]: ", " infections ",
           _sim_stats.get_all_num_infected(), ", dead epicells ",
           perc_str(_sim_stats.get_all_num_dead_epicells(), tissue.get_num_grid_points()),
           ", viral kills ", _sim_stats.get_all_num_viral_kills(),
           ", t-cells: ", _sim_stats.get_all_num_tcells(), "\n");
      time_step_ticks += (double)_options->num_timesteps / 20;
    }
    // sample
    if (time_step % _options->sample_period == 0) {
      sample(time_step, tissue, ViewObject::TCELL_VAS);
      sample(time_step, tissue, ViewObject::TCELL_XVAS);
      sample(time_step, tissue, ViewObject::VIRUS);
      sample(time_step, tissue, ViewObject::EPICELL);
      sample(time_step, tissue, ViewObject::ICYTOKINE);
      sample(time_step, tissue, ViewObject::CHEMOKINE);
    }
  }
  _update_tcell_timer.done_all();
  _update_virus_timer.done_all();
  _sample_timer.done_all();
  chrono::duration<double> t_elapsed = NOW() - start_t;
  SLOG("Finished ", _options->num_timesteps, " time steps in ", setprecision(4), fixed,
       t_elapsed.count(), " s (", (double)t_elapsed.count() / _options->num_timesteps,
       " s per step)\n");
}

int main(int argc, char **argv) {
  upcxx::init();
  auto start_t = NOW();
  _options = make_shared<Options>();
  if (!_options->load(argc, argv)) return 0;
  ProgressBar::SHOW_PROGRESS = _options->show_progress;

  _rnd_gen = make_shared<Random>(_options->rnd_seed + rank_me());

  if (pin_thread(getpid(), local_team().rank_me()) == -1)
    WARN("Could not pin process ", getpid(), " to core ", rank_me());
  else
    SLOG_VERBOSE("Pinned processes, with process 0 (pid ", getpid(), ") pinned to core ",
                 local_team().rank_me(), "\n");

  MemoryTrackerThread memory_tracker;
  memory_tracker.start();
  auto start_free_mem = get_free_mem();
  SLOG(KBLUE, "Starting with ", get_size_str(start_free_mem), " free on node 0", KNORM, "\n");
  Tissue tissue;
  tissue.construct({_options->dimensions[0], _options->dimensions[1], _options->dimensions[2]});
  initial_infection(tissue);
  finish_round(tissue);
  SLOG(KBLUE, "Memory used on node 0 after initialization is  ",
       get_size_str(start_free_mem - get_free_mem()), KNORM, "\n");

  run_sim(tissue);

  memory_tracker.stop();
  chrono::duration<double> t_elapsed = NOW() - start_t;
  SLOG("Finished in ", setprecision(2), fixed, t_elapsed.count(), " s at ", get_current_time(),
       " for SimCov version ", SIMCOV_VERSION, "\n");
  barrier();
  upcxx::finalize();
  return 0;
}
