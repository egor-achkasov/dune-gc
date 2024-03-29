// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>   // for std::cout
#include <vector>     // for std::vector
#include <cstdint>    // for uint, uint64_t
#include <chrono>     // for std::chrono::[high_resolution_clock, duration]
#include <set>        // for std::set
#include <unordered_map>  // for std::unordered_map

#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/operators.hh> // for MatrixAdapter
#include <dune/istl/preconditioners.hh> // for SeqILU
#include <dune/istl/solvers.hh> // for CGSolver
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>

// Global typedefs
typedef std::vector<uint> Colors;

/// A wrapper for a pair of a Dune Grid and colored DOFs within it
class ColoredGrid
{ /* empty for now */ };

// GC Algorithms
//=========================================================================

enum DOF_LOC { VERT, EDGE, BOTH };

template <class GridView>
class Topology
{
public:
  typedef GridView::IndexSet::IndexType IndexType;

  /// O(gv.size(0))
  Topology(const DOF_LOC df, const Basis& b)
  {
    // Init
    e2dof.reserve(gv.size(0));
    dof2e.reserve(gv.size(???));
    auto index_set = gv.indexSet();
    constexpr int DOF_CODIM = df == EDGE ? 2 : ???;

    // Iter over elements and their subentities
    for (const auto e& : Dune::elements(gv));
      for (const auto dof& : Dune::subEntities<DOF_CODIM>(e)){
        e2dof[index_set.index(e)] = index_set.index(dof);
        dof2e[index_set.index(dof)] = index_set.index(e);
      }
  }

  std::unordered_map<IndexType, IndexType> e2dof;
  std::unordered_map<IndexType, IndexType> dof2e;
};

/** Greedy 2D graph coloring. https://doi.org/10.1145/3458744.3473362
* \param g 2D grid
* \return {c}_i^N, where c -- uint color, i -- vertex index, N -- number of vertices
*/
template <class GridView>
Colors gc_greedy(const GridView& gv)
{
  // Init
  size_t num_verts = gv.size(GridView::dimension);
  Topology<GridView> t = Topology<GridView>(gv, ???);

  // A1 Iterative Parallel Graph Coloring
  // 1.1
  Colors C;
  C.resize(num_verts, 0);
  // 1.2
  std::set<uint> conf;
  auto hint = conf.begin();
  for (size_t v = 0; v != num_verts; ++v)
    hint = conf.insert(hint, v);
  // 1.3
  while (!conf.empty()){
    // 1.4
    // A2 Assign Colors
    // 2.1
    std::vector<uint64_t> forbidden;
    forbidden.resize(num_verts / 64 + 1);
    // 2.2
    for (uint v : conf){ // TODO parallel
      // 2.3
      for (uint64_t& f : forbidden) f = 0;
      // 2.4
      for (uint v_adj_id : g.vertices[v].adj)
        forbidden[v_adj_id / 64] |= 1 << (v_adj_id % 64);
      // 2.5
      for (uint i = 0; i != num_verts; ++i){
        if ((forbidden[i / 64] >> (i % 64)) % 2) continue;
        C[v] = i;
        break;
      }
    }

    // 1.5
    // A3 DetectConflicts
    // 3.1
    std::set<uint> new_conf;
    // 3.2
    for (uint v : conf) // TODO parallel
      for (uint u : g.adj(v))
        if (u < v && C[v] == C[u])
          new_conf.insert(v); // TODO atomic
    // 3.9
    conf = new_conf;
  }

  return C;
}

//=========================================================================
// end GC Algorithms

// Benchmarking
//=========================================================================

using CHRONO_CLOCK = std::chrono::high_resolution_clock;
#define CHRONO_DURCAST_MILS std::chrono::duration_cast<std::chrono::milliseconds>
static CHRONO_CLOCK::time_point g_chrono_tp1;
static int64_t g_chrono_dur_mils;
/// Measures and prints execution time of a function call
#define TIME(fcall) \
  std::cout << "Run \"" << #fcall << "\":\t"; \
  g_chrono_tp1 = CHRONO_CLOCK::now(); \
  fcall; \
  g_chrono_dur_mils = CHRONO_DURCAST_MILS(CHRONO_CLOCK::now() - g_chrono_tp1).count(); \
  std::cout << g_chrono_dur_mils << " ms" << std::endl;

// TODO memory
// TODO valgrid(massif) / perf

//=========================================================================
// end Benchmarking


// FEM algorithms
//=========================================================================

template<class LocalView, class Matrix>
extern void
assembleElementStiffnessMatrix(
  const LocalView& localView,
  Matrix& elementMatrix);

template<class LocalView>
extern void
assembleElementVolumeTerm(
  const LocalView& localView,
  Dune::BlockVector<double>& localB,
  const std::function<
    double(
      Dune::FieldVector<
        double, 
        LocalView::Element::dimension
      >
    )
  > volumeTerm);

template<class Basis>
extern void
getOccupationPattern(
  const Basis& basis,
  Dune::MatrixIndexSet& nb);

template<class Basis>
extern void
assemblePoissonProblem(
  const Basis& basis,
  Dune::BCRSMatrix<double>& matrix,
  Dune::BlockVector<double>& b,
  const std::function<
    double(
      Dune::FieldVector<
        double,
        Basis::GridView::dimension
      >
    )
  > volumeTerm);

//=========================================================================
// end FEM algorithms

int main(int argc, char** argv)
{
  // Init
  //=======================================================================
  // init main consts
  const int dim = 2;
  // init main typedefs
  typedef Dune::YaspGrid<dim>       Grid;
  typedef Grid::LeafGridView        GridView;
  typedef Topology<GridView>        Topology;
  typedef Dune::BCRSMatrix<double>  Matrix;
  typedef Dune::BlockVector<double> Vector;
  // init grid
  Dune::FieldVector<double,dim> len; for (auto& l : len) l=1.0;
  std::array<int,dim> cells; for (auto& c : cells) c=4;
  Grid grid(len,cells);
  grid.globalRefine(2);
  GridView gv = grid.leafGridView();
  // init FEM 1
  //-----------------------------------------------------------------------
  // Poisson Equation problem:
  // - \Delta u = -5 on (0,1)^2 \ [0.5,1)^2 with
  //   u = 0   on {0}x[0,1] \cup [0,1]x{0}
  //   u = 0.5 on {0.5}x[0.5,1] \cup [0.5,1]x{0.5}
  //   u' = 0  otherwise
  Dune::Functions::LagrangeBasis<GridView,1> basis(gv);
  std::vector<bool> dirichletNodes;
  auto predicate = [](auto x)
  {
    return x[0] < 1e-8
        || x[1] < 1e-8
        || (x[0] > 0.4999 && x[1] > 0.4999);
  };
  Dune::Functions::interpolate(basis, dirichletNodes, predicate);
  // FEM1 -- init Stiffness Matrix
  Matrix stiffnessMatrix1;
  for (size_t i=0; i<stiffnessMatrix1.N(); i++){
    if (!dirichletNodes[i]) continue;
    auto cIt    = stiffnessMatrix1[i].begin();
    auto cEndIt = stiffnessMatrix1[i].end();
    for (; cIt!=cEndIt; ++cIt)
      *cIt = (cIt.index()==i) ? 1.0 : 0.0;
  }
  // FEM1 -- init Dirichlet values
  Vector b1;
  auto dirichletValues1 = [](auto x)
  { return (x[0]< 1e-8 || x[1] < 1e-8) ? 0 : 0.5; };
  Dune::Functions::interpolate(basis, b1, dirichletValues1);
  // init solution vector
  Vector x1(basis.size());
  x1 = b1;
  // init linear operator
  Dune::MatrixAdapter<Matrix,Vector,Vector> linearOperator(stiffnessMatrix1);
  // init ILU
  Dune::SeqILU<Matrix,Vector,Vector> preconditioner(stiffnessMatrix1, 1.0);
  // init conjugate gradient solver
  Dune::CGSolver<Vector> cg1(linearOperator, preconditioner,
                      1e-5, // Desired residual reduction factor
                      50,   // Maximum number of iterations
                      2);   // Verbosity of the solver
  // init solving statistics object
  Dune::InverseOperatorResult statistics1;
  //-----------------------------------------------------------------------
  // end init FEM 1
  // init FEM 2
  //-----------------------------------------------------------------------
  // FEM 2
  // TODO
  //-----------------------------------------------------------------------
  //=======================================================================
  // end Init
  
  // Perform GC
  TIME(Topology topology = Topology(BOTH, basis));
  TIME(Colors clrs1 = gc_greedy(grid));
  // TODO TIME(Colors clrs2 = gc_???(grid))

  // Perform common FEM evaluation
  TIME(cg1.apply(x1, b1, statistics1))
  // TIME(cg2.apply(x2, b2, statistics2))

  // Perform parallel colored FEM evaluation
  TIME(cg1.apply(x1, b1, statistics1, clrs1))
  // TIME(cg2.apply(x2, b2, statistics2, clrs2))

  return 0;
}
