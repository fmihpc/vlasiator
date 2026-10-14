#include "sysboundary/outflow.h"
#include "object_wrapper.h"
#include "readparameters.h"
#include "velocity_mesh_parameters.h"
#include "grid.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#ifdef USE_GPU
#error This regression requires a CPU build.
#endif

class TestOutflow : public SBC::Outflow {
 public:
   const SBC::OutflowSpeciesParameters &species(std::size_t i) const { return speciesParams.at(i); }
   static uint copyScheme() { return COPY; }
   static uint noneScheme() { return NONE; }
   void setInterior(CellID target, CellID source) { allClosestNonsysboundaryCells[target] = {source}; }
   void setNonperiodic() { periodic.fill(false); }
   void setLayerOne(CellID target, CellID source) { allClosestL1OutflowCells[target] = {source}; }
};

bool configured(const std::string &name, const std::string &scenario, std::size_t face) {
   if (scenario == "schemes")
      return true;
   if (name == "proton" || scenario == "identical")
      return face == 1;
   return scenario != "empty" && (face == 4 || face == 5);
}

bool copies(const std::string &name, const std::string &scenario, std::size_t face) {
   if (scenario == "schemes")
      return (face % 2 == 0) == (name == "proton");
   return !(name == "proton" && face == 1);
}

std::array<bool, 6> geometricFaces(const std::array<double, 3> &center) {
   std::array<bool, 6> faces{};
   for (std::size_t d = 0; d < 3; ++d) {
      faces[2 * d] = center[d] > 6;
      faces[2 * d + 1] = center[d] < 2;
   }
   return faces;
}

bool copiesOnCell(const std::string &name, const std::string &scenario, const std::array<bool, 6> &faces) {
   for (std::size_t face = 0; face < 6; ++face) {
      if (faces[face] && configured(name, scenario, face) && copies(name, scenario, face))
         return true;
   }
   return false;
}

template <typename Check>
void distributedCheck(TestOutflow &boundary, const std::vector<std::string> &names, int rank, Check check) {
   dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> grid;
   grid.set_initial_length({{8, 8, 8}})
       .set_neighborhood_length(2)
       .set_maximum_refinement_level(0)
       .set_load_balancing_method("NONE")
       .set_periodic(false, false, false)
       .initialize(MPI_COMM_WORLD);
   dccrg::Cartesian_Geometry::Parameters geometry;
   for (int d = 0; d < 3; ++d) {
      geometry.start[d] = 0;
      geometry.level_0_cell_length[d] = 1;
   }
   grid.set_geometry(geometry);
   std::vector<std::array<int, 3>> axial, extended;
   for (int z = -2; z <= 2; ++z)
      for (int y = -2; y <= 2; ++y)
         for (int x = -2; x <= 2; ++x) {
            if (x || y || z) {
               extended.push_back({{x, y, z}});
               if ((x != 0) + (y != 0) + (z != 0) == 1)
                  axial.push_back({{x, y, z}});
            }
         }
   check(grid.add_neighborhood(Neighborhoods::VLASOV_SOLVER, axial), "distributed axial stencil");
   check(grid.add_neighborhood(Neighborhoods::SYSBOUNDARIES_EXTENDED, extended), "distributed transfer stencil");
   const auto cells = grid.get_cells();
   auto value = [](CellID id, std::size_t pop, std::size_t v) -> Realf { return 10000 + 100 * pop + 2 * id + v; };
   for (const auto id : cells) {
      const auto faces = geometricFaces(grid.geometry.get_center(id));
      grid[id]->sysBoundaryFlag = std::any_of(faces.begin(), faces.end(), [](bool b) { return b; })
                                      ? boundary.getIndex()
                                      : sysboundarytype::NOT_SYSBOUNDARY;
      grid[id]->sysBoundaryLayer = 1;
      for (std::size_t pop = 0; pop < names.size(); ++pop) {
         check(grid[id]->add_velocity_block(0, pop), "distributed block allocation");
         for (std::size_t v = 0; v < WID3; ++v)
            grid[id]->get_data(pop)[v] = value(id, pop, v);
         grid[id]->get_population(pop).RHO = 1000 + 100 * pop + id;
         grid[id]->get_population(pop).N_blocks = grid[id]->size(pop);
      }
   }
   // Transfer through production block-list preparation and VDF communication,
   // rather than populating remote cells directly in the test.
   SpatialCell::set_mpi_transfer_type(Transfer::ALL_SPATIAL_DATA);
   grid.update_copies_of_remote_neighbors(Neighborhoods::SYSBOUNDARIES_EXTENDED);
   for (std::size_t pop = 0; pop < names.size(); ++pop) {
      updateRemoteVelocityBlockLists(grid, pop, Neighborhoods::SYSBOUNDARIES_EXTENDED);
      SpatialCell::set_mpi_transfer_type(Transfer::VEL_BLOCK_DATA);
      grid.update_copies_of_remote_neighbors(Neighborhoods::SYSBOUNDARIES_EXTENDED);
   }
   int remoteTargets = 0;
   for (const auto target : cells) {
      if (grid[target]->sysBoundaryFlag != boundary.getIndex())
         continue;
      CellID source = INVALID_CELLID;
      for (const auto &neighbor : *grid.get_neighbors_of(target, Neighborhoods::SYSBOUNDARIES_EXTENDED)) {
         if (neighbor.first && grid.get_process(neighbor.first) != rank && grid[neighbor.first] &&
             grid[neighbor.first]->sysBoundaryFlag == sysboundarytype::NOT_SYSBOUNDARY) {
            source = neighbor.first;
            break;
         }
      }
      if (source == INVALID_CELLID)
         continue;
      ++remoteTargets;
      check(grid.get_process(source) != rank, "source must be owned by another rank");
      boundary.setInterior(target, source);
      for (std::size_t calledPop = 0; calledPop < names.size(); ++calledPop) {
         for (std::size_t pop = 0; pop < names.size(); ++pop) {
            check(grid[source]->size(pop) == 1, "remote block list received");
            for (std::size_t v = 0; v < WID3; ++v) {
               check(grid[source]->get_data(pop)[v] == value(source, pop, v), "remote VDF received");
               grid[target]->get_data(pop)[v] = value(target, pop, v);
            }
            grid[target]->get_population(pop).RHO = 1000 + 100 * pop + target;
         }
         boundary.vlasovBoundaryCondition(grid, target, calledPop, false);
         for (std::size_t pop = 0; pop < names.size(); ++pop) {
            const bool copied = pop == calledPop &&
                                copiesOnCell(names[pop], "schemes", geometricFaces(grid.geometry.get_center(target)));
            for (std::size_t v = 0; v < WID3; ++v) {
               check(grid[target]->get_data(pop)[v] == value(copied ? source : target, pop, v),
                     "distributed population VDF");
               check(grid[source]->get_data(pop)[v] == value(source, pop, v), "distributed source unchanged");
            }
            check(grid[target]->get_population(pop).RHO == 1000 + 100 * pop + (copied ? source : target),
                  "distributed density");
         }
      }
   }
   check(remoteTargets > 0, "at least one genuinely remote source on each rank");
   std::cout << "OUTFLOW_REMOTE rank=" << rank << " targets=" << remoteTargets << '\n';
}

int run(int argc, char **argv, int rank) {
   const std::vector<std::string> scenarios{"single", "different", "reversed", "identical", "empty", "schemes"};
   if (argc < 3 || std::find(scenarios.begin(), scenarios.end(), argv[1]) == scenarios.end())
      return 1;
   const std::string scenario = argv[1];
   std::vector<char *> arguments{argv[0]};
   for (int i = 2; i < argc; ++i)
      arguments.push_back(argv[i]);
   const std::vector<std::string> names = scenario == "single"     ? std::vector<std::string>{"proton"}
                                          : scenario == "reversed" ? std::vector<std::string>{"helium", "proton"}
                                                                   : std::vector<std::string>{"proton", "helium"};
   auto &populations = getObjectWrapper().particleSpecies;
   populations.resize(names.size());
   for (std::size_t i = 0; i < names.size(); ++i) {
      populations[i].name = names[i];
      populations[i].velocityMesh = 0;
      populations[i].sparseMinValue = 0;
   }
   Readparameters parameters(arguments.size(), arguments.data());
   TestOutflow boundary;
   boundary.addParameters();
   std::vector<std::string> invalid, files;
   parameters.parse(invalid, files, false);
   boundary.getParameters();
   boundary.setNonperiodic();
   int failures = 0;
   int assertions = 0;
   auto check = [&](bool condition, const std::string &label) {
      ++assertions;
      if (!condition) {
         ++failures;
         std::cerr << "rank " << rank << " scenario " << scenario << ": " << label << '\n';
      }
   };
   check(invalid.empty(), "parser rejected fixture");
   std::array<bool, 6> expectedGlobal{}, actualGlobal{};
   for (std::size_t pop = 0; pop < names.size(); ++pop) {
      for (std::size_t face = 0; face < 6; ++face) {
         const bool active = configured(names[pop], scenario, face);
         expectedGlobal[face] = expectedGlobal[face] || active;
         check(boundary.species(pop).facesToSkipVlasov[face] == !active, "face mask");
         check(boundary.species(pop).faceVlasovScheme[face] ==
                   (copies(names[pop], scenario, face) ? TestOutflow::copyScheme() : TestOutflow::noneScheme()),
               "face scheme");
      }
   }
   boundary.getFaces(actualGlobal.data());
   check(actualGlobal == expectedGlobal, "geometric face union");

   vmesh::allocateMeshWrapper();
   vmesh::MeshParameters mesh;
   for (int d = 0; d < 3; ++d) {
      mesh.meshLimits[2 * d] = -1;
      mesh.meshLimits[2 * d + 1] = 1;
      mesh.gridLength[d] = 1;
      mesh.blockLength[d] = WID;
   }
   vmesh::getMeshWrapper()->velocityMeshesCreation->push_back(&mesh);
   vmesh::getMeshWrapper()->initVelocityMeshes(1);
   Parameters::xcells_ini = Parameters::ycells_ini = Parameters::zcells_ini = 8;
   {
      dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> grid;
      grid.set_initial_length({{8, 8, 8}})
          .set_neighborhood_length(2)
          .set_maximum_refinement_level(0)
          .set_load_balancing_method("NONE")
          .set_periodic(false, false, false)
          .initialize(MPI_COMM_SELF);
      dccrg::Cartesian_Geometry::Parameters geometry;
      for (int d = 0; d < 3; ++d) {
         geometry.start[d] = 0;
         geometry.level_0_cell_length[d] = 1;
      }
      grid.set_geometry(geometry);
      std::vector<std::array<int, 3>> stencil;
      for (int d = -2; d <= 2; ++d) {
         if (d) {
            stencil.push_back({{d, 0, 0}});
            stencil.push_back({{0, d, 0}});
            stencil.push_back({{0, 0, d}});
         }
      }
      check(grid.add_neighborhood(Neighborhoods::VLASOV_SOLVER, stencil), "register Vlasov stencil");
      std::vector<std::array<int, 3>> nearest, extended;
      for (int z = -2; z <= 2; ++z)
         for (int y = -2; y <= 2; ++y)
            for (int x = -2; x <= 2; ++x) {
               if (x || y || z) {
                  extended.push_back({{x, y, z}});
                  if (std::abs(x) <= 1 && std::abs(y) <= 1 && std::abs(z) <= 1)
                     nearest.push_back({{x, y, z}});
               }
            }
      check(grid.add_neighborhood(Neighborhoods::SYSBOUNDARIES, nearest), "register boundary stencil");
      check(grid.add_neighborhood(Neighborhoods::SYSBOUNDARIES_EXTENDED, extended), "register extended stencil");
      const auto cells = grid.get_cells();
      for (const auto id : cells)
         grid[id]->sysBoundaryFlag = sysboundarytype::NOT_SYSBOUNDARY;
      auto findCell = [&](const std::array<double, 3> &center) -> CellID {
         for (const auto id : cells)
            if (grid.geometry.get_center(id) == center)
               return id;
         return INVALID_CELLID;
      };
      const CellID source = findCell({{3.5, 3.5, 3.5}});
      check(source != INVALID_CELLID, "interior cell exists");
      if (source == INVALID_CELLID)
         MPI_Abort(MPI_COMM_WORLD, 1);
      grid[source]->sysBoundaryFlag = sysboundarytype::NOT_SYSBOUNDARY;
      // Six face interiors, twelve edges and eight corners; None is not a veto
      // when another configured face on the same cell requests Copy.
      std::vector<std::array<double, 3>> centers;
      for (double z : {0.5, 3.5, 7.5})
         for (double y : {0.5, 3.5, 7.5})
            for (double x : {0.5, 3.5, 7.5}) {
               if (x != 3.5 || y != 3.5 || z != 3.5)
                  centers.push_back({{x, y, z}});
            }
      for (const auto &center : centers) {
         const CellID target = findCell(center);
         check(target != INVALID_CELLID, "boundary cell exists");
         if (target == INVALID_CELLID)
            MPI_Abort(MPI_COMM_WORLD, 1);
         grid[target]->sysBoundaryFlag = boundary.getIndex();
         boundary.setInterior(target, source);
         for (uint layer : {1u, 2u})
            for (bool velocityMoments : {false, true}) {
               grid[target]->sysBoundaryLayer = layer;
               for (std::size_t calledPop = 0; calledPop < names.size(); ++calledPop) {
                  for (std::size_t pop = 0; pop < names.size(); ++pop) {
                     if (grid[source]->size(pop) == 0)
                        check(grid[source]->add_velocity_block(0, pop), "source block allocation");
                     if (grid[target]->size(pop) == 0)
                        check(grid[target]->add_velocity_block(0, pop), "target block allocation");
                     for (std::size_t v = 0; v < WID3; ++v) {
                        grid[source]->get_data(pop)[v] = 1000 + 100 * pop + v;
                        grid[target]->get_data(pop)[v] = 10 + 100 * pop + v;
                     }
                     grid[source]->get_population(pop).RHO = 1000 + pop;
                     grid[target]->get_population(pop).RHO = 10 + pop;
                     grid[source]->get_population(pop).RHO_R = 2000 + pop;
                     grid[target]->get_population(pop).RHO_R = 20 + pop;
                     grid[source]->get_population(pop).RHO_V = 3000 + pop;
                     grid[target]->get_population(pop).RHO_V = 30 + pop;
                     grid[source]->get_population(pop).N_blocks = grid[source]->size(pop);
                     grid[target]->get_population(pop).N_blocks = grid[target]->size(pop);
                  }
                  boundary.vlasovBoundaryCondition(grid, target, calledPop, velocityMoments);
                  for (std::size_t pop = 0; pop < names.size(); ++pop) {
                     const bool copied = pop == calledPop && copiesOnCell(names[pop], scenario, geometricFaces(center));
                     const std::string label = names[pop] + " cell " + std::to_string(target) + " layer " +
                                               std::to_string(layer) + " after call " + std::to_string(calledPop);
                     check(grid[target]->size(pop) == 1, label + " block count");
                     for (std::size_t v = 0; v < WID3; ++v) {
                        check(grid[target]->get_data(pop)[v] == (copied && layer == 1 ? 1000 : 10) + 100 * pop + v,
                              label + " VDF");
                        check(grid[source]->get_data(pop)[v] == 1000 + 100 * pop + v, label + " source VDF");
                     }
                     check(grid[target]->get_population(pop).RHO == (copied ? 1000 : 10) + pop, label + " density");
                     check(grid[target]->get_population(pop).RHO_R ==
                               (copied && (layer == 1 || !velocityMoments) ? 2000 : 20) + pop,
                           label + " R density");
                     check(grid[target]->get_population(pop).RHO_V ==
                               (copied && (layer == 1 || velocityMoments) ? 3000 : 30) + pop,
                           label + " V density");
                     check(grid[source]->get_population(pop).RHO == 1000 + pop, label + " source density");
                  }
               }
            }
      }

      // Exercise the production restart helper, including both R/V moment passes.
      recalculateLocalCellsCache(grid);
      for (std::size_t face = 0; face < 6; ++face) {
         for (const auto id : cells)
            grid[id]->sysBoundaryFlag = sysboundarytype::NOT_SYSBOUNDARY;
         std::array<double, 3> center{{3.5, 3.5, 3.5}};
         center[face / 2] = face % 2 == 0 ? 6.5 : 1.5;
         const CellID l1 = findCell(center);
         center[face / 2] = face % 2 == 0 ? 7.5 : 0.5;
         const CellID l2 = findCell(center);
         grid[l1]->sysBoundaryFlag = grid[l2]->sysBoundaryFlag = boundary.getIndex();
         grid[l1]->sysBoundaryLayer = 1;
         grid[l2]->sysBoundaryLayer = 2;
         boundary.setInterior(l1, source);
         boundary.setLayerOne(l2, l1);
         for (std::size_t pop = 0; pop < names.size(); ++pop) {
            for (const CellID id : {source, l1, l2}) {
               if (!grid[id]->size(pop))
                  check(grid[id]->add_velocity_block(0, pop), "restart block allocation");
               std::fill_n(grid[id]->get_data(pop), WID3, id == source ? 1000 + pop : 10 + pop);
               auto &p = grid[id]->get_population(pop);
               p.RHO = id == source ? 1000 + pop : 10 + pop;
               p.RHO_R = id == source ? 2000 + pop : 20 + pop;
               p.RHO_V = id == source ? 3000 + pop : 30 + pop;
               p.N_blocks = grid[id]->size(pop);
            }
         }
         boundary.setupL2OutflowAtRestart(grid);
         for (std::size_t pop = 0; pop < names.size(); ++pop) {
            const bool copied = configured(names[pop], scenario, face) && copies(names[pop], scenario, face);
            for (std::size_t v = 0; v < WID3; ++v) {
               check(grid[l1]->get_data(pop)[v] == (copied ? 1000 : 10) + pop, "restart L1 VDF");
               check(grid[l2]->get_data(pop)[v] == 10 + pop, "restart L2 VDF unchanged");
            }
            for (const CellID id : {l1, l2}) {
               check(grid[id]->get_population(pop).RHO == (copied ? 1000 : 10) + pop, "restart density");
               check(grid[id]->get_population(pop).RHO_R == (copied ? 2000 : 20) + pop, "restart R density");
               check(grid[id]->get_population(pop).RHO_V == (copied ? 3000 : 30) + pop, "restart V density");
            }
         }
      }
   }
   int ranks = 0;
   MPI_Comm_size(MPI_COMM_WORLD, &ranks);
   if (scenario == "schemes" && ranks > 1)
      distributedCheck(boundary, names, rank, check);
   vmesh::getMeshWrapper()->velocityMeshesCreation->clear();
   delete vmesh::host_getMeshWrapper();
   std::cout << "OUTFLOW scenario=" << scenario << " rank=" << rank << " assertions=" << assertions
             << " failures=" << failures << '\n';
   return failures;
}

int main(int argc, char **argv) {
   MPI_Init(&argc, &argv);
   phiprof::initialize();
   int rank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   const int failures = run(argc, argv, rank);
   int total = 0;
   MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Finalize();
   return total ? 1 : 0;
}
