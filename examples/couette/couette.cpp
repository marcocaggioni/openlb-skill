// Couette rheometer — OpenLB case (2D annulus cross-section).
//
// Status: see examples/couette/RUN_REPORT.md for whether this has been
// compiled and run yet. OpenLB version/commit used: see RUN_REPORT.md.
//
// Physics: fluid between two concentric cylinders, inner wall rotating
// at Omega, outer wall fixed. Torque on the inner wall is computed by
// momentum exchange at the inner-wall links, converted to N*m per unit
// length. The post-processing script couette_post.py converts torque ->
// shear stress and Omega -> shear rate and enforces the Newtonian
// validation gate.
//
// Method (all API names verified against the cloned OpenLB 1.9 sources,
// following examples/forBeginners/cylinder2d/cylinder2d.cpp and
// examples/turbulence/aorta3d/aorta3d.cpp):
//  - D2Q9, BGK dynamics (lattice.defineDynamics<BGKdynamics>).
//  - Annulus: material 1 = fluid (Ri < r < Ro), material 2 =
//    outer wall, material 3 = inner wall (the bob).
//  - Outer wall: boundary::set<boundary::BounceBack> (material 2).
//  - Inner wall: Bouzidi interpolated bounce-back with imposed wall
//    velocity — setBouzidiBoundary<BouzidiVelocityPostProcessor>
//    (material 3, analytical inner circle) + setBouzidiVelocity with
//    a RotatingWallVelocity functor (u_wall = Omega x r, returned in
//    lattice units — setBouzidiVelocity expects lattice units,
//    aorta3d pattern). The wall velocity is updated during the
//    quarter-sine ramp in setTemporalValues.
//  - Torque: momentum exchange at the inner-wall links, the same
//    formula as OpenLB's own SuperLatticePhysBoundaryForce2D
//    (src/functors/lattice/latticePhysBoundaryForce2D.hh:87-113):
//    per inner-wall cell and per link to a fluid neighbour,
//    F -= c_i * (f_fluid[iPop] + f_wall[opp(iPop)]), summed over the
//    wall, torque = sum r x F. It is evaluated directly from the
//    block-lattice cells (unambiguous in serial) rather than through
//    the super-functor dispatch. After streaming, the two populations
//    in the sum are the reflected and the incoming one, so the
//    moving-wall correction is already contained in them (the
//    reflected population was built by the velocity boundary).
//    Conversion: 2D torque per unit length, T/L = torque_lb *
//    rho*dx^4/dt^2 (== converter.getPhysForce(torque_lb)), NOT
//    rho*dx^5/dt^2 (that would be a 3D torque). A second, independent
//    torque is computed from the simulated azimuthal velocity profile
//    (shear rate at the inner wall) and both are printed.
//  - Initial condition: rest; the ramp brings the wall up smoothly.
//    Viscous diffusion time gap^2/nu ~ 0.025 s, run time 1 s.
//
// Build (serial, no MPI):
//   cd examples/couette && make            (uses OpenLB's GNU Make build)
// Run:
//   ./couette --gamma 10 --out runs/gamma_10
//
// Defaults below are examples/couette/case.yaml (do not change the case
// by editing numbers here — change case.yaml and keep this file in sync).

#include <olb.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<NavierStokes, Lattice<double, descriptors::D2Q9<>>>;

namespace olb::parameters {

struct INNER_R : public descriptors::FIELD_BASE<1> {};
struct OUTER_R : public descriptors::FIELD_BASE<1> {};

} // namespace olb::parameters

// ---- case parameters (case.yaml) ----
static double Ri = 0.020;    // m, inner radius (rotating)
static double Ro = 0.025;    // m, outer radius (fixed)
static double Lphys = 0.060; // m, wetted length (only used by couette_post.py)
static double rhoPhys = 1000.0;
static double etaPhys = 1.0; // Pa.s, Newtonian validation case
static double gammaTarget = 10.0; // 1/s, from --gamma
static int N_gap = 10;       // cells across the gap
static double omegaTarget = 1.0; // relaxation frequency
static double runTime = 1.0; // s per point
static double rampTime = 0.2; // s, quarter-sine ramp of Omega
static std::string outDir = "runs/gamma_10";

// ---- derived lattice parameters (same formulas as lattice.py) ----
static double dx, dt, nuLb, Omega, torqueConv;
static std::size_t nSteps, rampSteps;
static int nx, ny; // grid size (cells)
static double cxPhys, cyPhys; // annulus centre, physical (m)
static const double PI = std::acos(-1.0);

// Rotating-wall velocity functor (PHYSICAL units — momenta::setVelocity
// scales to lattice units internally, src/case/setters.h). Omega is
// updated during the ramp via setOmega().
template <typename T>
class RotatingWallVelocity : public AnalyticalF2D<T, T> {
private:
  T _cx, _cy, _omega, _scale;
public:
  RotatingWallVelocity(T cx, T cy)
    : AnalyticalF2D<T, T>(2), _cx(cx), _cy(cy), _omega(T(0)), _scale(T(1)) {}
  void setOmega(T omega) { _omega = omega; }
  // setBouzidiVelocity expects the velocity in LATTICE units although
  // it evaluates the functor at physical positions (aorta3d pattern:
  // the functor there is built from converter.getLatticeVelocity).
  // Scale = dt/dx converts the physical wall velocity to lattice units.
  void setScale(T s) { _scale = s; }
  void setCenter(T cx, T cy) { _cx = cx; _cy = cy; }
  bool operator()(T output[], const T input[]) override
  {
    T x = input[0] - _cx, y = input[1] - _cy;
    output[0] = -_scale * _omega * y;
    output[1] = _scale * _omega * x;
    return true;
  }
};

static RotatingWallVelocity<double> wallVel(0.0, 0.0);

// Analytic Couette azimuthal velocity (m/s) at radius r (m).
static double couetteVtheta(double r)
{
  double A = Omega * Ri * Ri / (Ro * Ro - Ri * Ri);
  return A * (Ro * Ro / r - r);
}

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T = MyCase::value_t;
  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  const Vector<T, 2> extent{2 * Ro, 2 * Ro};
  const Vector<T, 2> origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extent, origin);
  Mesh<T, MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();

  const Vector<T, 2> center{Ro, Ro};
  IndicatorCircle2D<T> innerCircle(center, Ri);
  IndicatorCircle2D<T> outerCircle(center, Ro);

  // Material 1 = fluid everywhere, then carve out the walls.
  geometry.rename(0, 1);
  // Inner disc -> material 3 (rotating bob).
  geometry.rename(1, 3, innerCircle);
  // Everything outside the outer circle -> material 2 (fixed cup wall).
  auto domainBox = std::make_shared<IndicatorCuboid2D<T>>(
      Vector<T, 2>{2 * Ro, 2 * Ro}, Vector<T, 2>{0, 0});
  auto outside = domainBox - std::make_shared<IndicatorCircle2D<T>>(center, Ro);
  geometry.rename(1, 2, *outside);

  geometry.clean();
  geometry.checkForErrors();
  geometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes{});

  const T charU = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  lattice.setUnitConverter<UnitConverter<T, DESCRIPTOR>>(
      (T)dx,                 // physDeltaX
      (T)dt,                 // physDeltaT
      (T)(Ro - Ri),          // charPhysLength: the gap
      (T)charU,              // charPhysVelocity: max wall speed this run
      (T)(etaPhys / rhoPhys),// physViscosity: kinematic viscosity
      (T)rhoPhys             // physDensity
  );
  const auto& converter = lattice.getUnitConverter();
  converter.print();
  const T omegaConv = converter.getLatticeRelaxationFrequency();
  clout << "converter omega = " << omegaConv << " (target " << omegaTarget << ")"
        << std::endl;

  // Material 1: bulk BGK dynamics.
  lattice.defineDynamics<BGKdynamics>(geometry, 1);
  // Material 2: fixed outer wall, bounce-back.
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);
  // Material 3: rotating inner wall, Bouzidi interpolated bounce-back
  // with imposed wall velocity (BouzidiVelocityPostProcessor), the
  // moving-wall pattern of examples/turbulence/aorta3d/aorta3d.cpp
  // (setBouzidiBoundary + setBouzidiVelocity). The analytical circle
  // gives the exact wall position, so the wall is NOT a staircase.
  // NOTE: boundary::set<boundary::InterpolatedVelocity> does NOT work
  // here — it iterates over every material-3 cell including the solid
  // disc interior, where the discrete normal is (0,0), and throws
  // "Could not set Boundary" (src/boundary/setBoundary2D.h:127).
  IndicatorCircle2D<T> innerCircle(Vector<T, 2>{Ro, Ro}, Ri);
  setBouzidiBoundary<T, DESCRIPTOR, BouzidiVelocityPostProcessor>(
      lattice, geometry, 3, innerCircle);

  lattice.setParameter<descriptors::OMEGA>(omegaConv);
  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "setInitialValues");
  auto& lattice = myCase.getLattice(NavierStokes{});
  // Rest initial condition; the ramp in setTemporalValues spins the wall up.
  lattice.initialize();
  clout << "Set Initial Values ... OK" << std::endl;
}

// Update the wall velocity during the quarter-sine ramp. Called every
// step before collideAndStream; only re-sets the velocity while the
// ramp is running (afterwards Omega is constant).
void setTemporalValues(MyCase& myCase, std::size_t iT)
{
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice  = myCase.getLattice(NavierStokes{});

  if (iT > rampSteps) return;
  if (iT % 100 != 0 && iT != rampSteps) return;
  const T frac = (iT >= rampSteps) ? T(1)
      : std::sin(0.5 * PI * (T)iT / (T)rampSteps);
  wallVel.setOmega(Omega * frac);
  // Bouzidi wall velocity, in lattice units (functor scale set in
  // main). Evaluated at the physical wall-intersection point of each
  // missing link (src/boundary/setBouzidiBoundary.h:383-432).
  setBouzidiVelocity(lattice, geometry, 3, wallVel);
  // Push the updated boundary velocity to the lattice (aorta3d pattern).
  lattice.setProcessingContext<Array<descriptors::BOUZIDI_VELOCITY>>(
      ProcessingContext::Simulation);
}

// Torque on the inner wall (N*m per unit length) from momentum exchange,
// evaluated directly from the block-lattice cells. Same formula as
// OpenLB's SuperLatticePhysBoundaryForce2D
// (src/functors/lattice/latticePhysBoundaryForce2D.hh:87-113).
// Must be called AFTER collideAndStream (populations just streamed).
static double computeTorque(MyCase& myCase, double rampFrac)
{
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice  = myCase.getLattice(NavierStokes{});
  auto& block    = lattice.getBlock(0);
  auto& blockGeo = geometry.getBlockGeometry(0);

  double torqueLb = 0.0;
  for (int iY = 0; iY < ny; ++iY) {
    for (int iX = 0; iX < nx; ++iX) {
      if (blockGeo.get({iX, iY}) != 3) continue;
      // Position of the cell centre relative to the annulus centre,
      // in lattice units.
      const double xRel = (iX + 0.5) - nx / 2.0;
      const double yRel = (iY + 0.5) - ny / 2.0;
      double Fx = 0.0, Fy = 0.0;
      for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
        const int nX = iX + descriptors::c<DESCRIPTOR>(iPop, 0);
        const int nY = iY + descriptors::c<DESCRIPTOR>(iPop, 1);
        if (nX < 0 || nX >= nx || nY < 0 || nY >= ny) continue;
        if (blockGeo.get({nX, nY}) != 1) continue; // neighbour must be fluid
        const int iOpp = descriptors::opposite<DESCRIPTOR>(iPop);
        const double f = block.get(nX, nY)[iPop] + block.get(iX, iY)[iOpp];
        Fx -= descriptors::c<DESCRIPTOR>(iPop, 0) * f;
        Fy -= descriptors::c<DESCRIPTOR>(iPop, 1) * f;
      }
      torqueLb += xRel * Fy - yRel * Fx;
    }
  }
  (void)rampFrac;
  // Sign convention: the raw sum is the torque the fluid exerts ON the
  // bob, which opposes the rotation (negative for Omega > 0). Return
  // the driving torque instead — the torque that must be applied to
  // the bob to sustain the rotation — so it is positive and directly
  // comparable with the analytic Couette torque and with couette_post.py.
  return -torqueLb * torqueConv;
}

void simulate(MyCase& myCase)
{
  OstreamManager clout(std::cout, "simulate");
  using T = MyCase::value_t;
  auto& lattice   = myCase.getLattice(NavierStokes{});
  auto& geometry  = myCase.getGeometry();
  const auto& converter = lattice.getUnitConverter();

  std::filesystem::create_directories(outDir);
  std::ofstream torqueFile(outDir + "/torque.csv");
  torqueFile << "time_s,torque_N_m_per_m\n";

  // Radial profile bins (1 lattice unit wide), accumulated in steady window.
  const std::size_t steadyStart = (std::size_t)(nSteps * 0.8);
  const int nBins = (int)std::ceil(Ro / dx) + 1;
  std::vector<double> profSum(nBins, 0.0), profCnt(nBins, 0.0);

  const std::size_t sampleEvery = 100, fieldEvery = 1000;
  double torqueSteadySum = 0.0;
  std::size_t torqueSteadyCnt = 0;
  int fieldIdx = 0;

  auto writeField = [&](std::size_t /*iT*/) {
    std::ofstream f(outDir + "/field_" + std::to_string(fieldIdx++) + ".csv");
    f << "iX,iY,ux_m_s,uy_m_s\n";
    auto& block = lattice.getBlock(0);
    for (int iY = 0; iY < ny; iY += 2) {
      for (int iX = 0; iX < nx; iX += 2) {
        double rho = 0.0;
        double u[2] = {0.0, 0.0};
        block.get(iX, iY).computeRhoU(rho, u);
        f << iX << "," << iY << "," << u[0] * dx / dt << "," << u[1] * dx / dt
          << "\n";
      }
    }
  };

  clout << "starting simulation... " << nSteps << " steps" << std::endl;
  for (std::size_t iT = 0; iT < nSteps; ++iT) {
    setTemporalValues(myCase, iT);
    lattice.collideAndStream();

    if (iT % sampleEvery == 0 || iT == nSteps - 1) {
      const double rampFrac = (iT >= rampSteps) ? 1.0
          : std::sin(0.5 * PI * (double)iT / (double)rampSteps);
      const double Tq = computeTorque(myCase, rampFrac);
      torqueFile << iT * dt << "," << Tq << "\n";
      if (iT >= steadyStart) {
        torqueSteadySum += Tq;
        ++torqueSteadyCnt;
      }
    }
    if (iT >= steadyStart && iT % sampleEvery == 0) {
      auto& block = lattice.getBlock(0);
      auto& blockGeo = geometry.getBlockGeometry(0);
      for (int iY = 0; iY < ny; ++iY) {
        for (int iX = 0; iX < nx; ++iX) {
          if (blockGeo.get({iX, iY}) != 1) continue;
          const double x = (iX + 0.5) - nx / 2.0;
          const double y = (iY + 0.5) - ny / 2.0;
          const double r = std::sqrt(x * x + y * y);
          if (r <= 0.0) continue;
          double rho = 0.0;
          double u[2] = {0.0, 0.0};
          block.get(iX, iY).computeRhoU(rho, u);
          const double vThetaLb = (-y * u[0] + x * u[1]) / r;
          const int b = (int)r;
          if (b >= 0 && b < nBins) {
            profSum[b] += vThetaLb * dx / dt; // m/s
            profCnt[b] += 1.0;
          }
        }
      }
    }
    if (iT % fieldEvery == 0) writeField(iT);
    if (iT % 4000 == 0) {
      clout << "step " << iT << "/" << nSteps << " t=" << iT * dt << " s"
            << std::endl;
    }
  }
  writeField(nSteps - 1);

  // Profile from the steady window.
  std::ofstream profFile(outDir + "/profile.csv");
  profFile << "r_m,vtheta_m_s,vtheta_analytic_m_s\n";
  double tauFromProfile = 0.0;
  {
    // Shear rate at the inner wall from the two innermost fluid bins.
    const int b0 = (int)(Ri / dx) + 1, b1 = b0 + 1;
    const double r0 = (b0 + 0.5) * dx, r1 = (b1 + 0.5) * dx;
    const double v0 = profSum[b0] / profCnt[b0];
    const double v1 = profSum[b1] / profCnt[b1];
    const double omega0 = v0 / r0, omega1 = v1 / r1;
    const double domegaDr = (omega1 - omega0) / (r1 - r0);
    // tau_rtheta = eta * r * d(omega)/dr is NEGATIVE (omega falls
    // with r); the torque magnitude uses -tau_rtheta.
    tauFromProfile = -etaPhys * r0 * domegaDr;
    clout << "profile: vtheta(r0=" << r0 << ")=" << v0 << " m/s, r1=" << r1
          << " v1=" << v1 << ", tau_inner=" << tauFromProfile << " Pa"
          << std::endl;
  }
  for (int b = 0; b < nBins; ++b) {
    if (profCnt[b] == 0) continue;
    const double r = (b + 0.5) * dx;
    profFile << r << "," << profSum[b] / profCnt[b] << "," << couetteVtheta(r)
             << "\n";
  }

  const double TqSteady = torqueSteadyCnt ? torqueSteadySum / torqueSteadyCnt : 0.0;
  const double Tanalytic = 4 * PI * Lphys * etaPhys * Omega * Ri * Ri * Ro * Ro
                           / (Ro * Ro - Ri * Ri);
  clout << "steady torque (momentum exchange): " << TqSteady * Lphys
        << " N*m (per-length " << TqSteady << " N*m/m)" << std::endl;
  clout << "analytic torque: " << Tanalytic << " N*m, rel err "
        << std::abs(TqSteady * Lphys - Tanalytic) / Tanalytic * 100 << " %"
        << std::endl;
  clout << "torque from profile: "
        << tauFromProfile * 2 * PI * Ri * Ri * Lphys << " N*m" << std::endl;
  (void)converter;
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--gamma" && i + 1 < argc) gammaTarget = atof(argv[++i]);
    else if (arg == "--out" && i + 1 < argc) outDir = argv[++i];
  }

  // Lattice parameters (identical to lattice.py).
  const double gap = Ro - Ri;
  dx = gap / N_gap;
  const double nuPhys = etaPhys / rhoPhys;
  nuLb = (1.0 / omegaTarget - 0.5) / 3.0;
  dt = nuLb * dx * dx / nuPhys;
  Omega = gammaTarget * (Ro * Ro - Ri * Ri) / (2.0 * Ro * Ro);
  nSteps = (std::size_t)std::ceil(runTime / dt);
  rampSteps = (std::size_t)std::ceil(rampTime / dt);
  torqueConv = rhoPhys * std::pow(dx, 4) / (dt * dt); // lattice torque -> N*m/m
  nx = (int)std::ceil(2.0 * Ro / dx);
  ny = nx;
  cxPhys = Ro;
  cyPhys = Ro;
  wallVel.setCenter(cxPhys, cyPhys);
  wallVel.setScale(dt / dx); // physical -> lattice velocity for Bouzidi

  clout << "Couette (OpenLB): gamma=" << gammaTarget << " 1/s, Omega=" << Omega
        << " rad/s, dx=" << dx << " m, dt=" << dt << " s, omega=" << omegaTarget
        << ", steps=" << nSteps << ", grid=" << nx << "x" << ny << std::endl;

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    using T = MyCase::value_t;
    myCaseParameters.set<RESOLUTION>(N_gap);
    myCaseParameters.set<PHYS_DELTA_X>(dx);
    myCaseParameters.set<PHYS_MAX_T>(runTime);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(gap);
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(Omega * Ri);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(nuPhys);
    myCaseParameters.set<PHYS_DENSITY>(rhoPhys);
    myCaseParameters.set<OVERLAP>(1);
    myCaseParameters.set<INNER_R>(Ri);
    myCaseParameters.set<OUTER_R>(Ro);
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
