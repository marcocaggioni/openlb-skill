// Lid-driven cavity — OpenLB case, THREE-COMPONENT (TC) fluid (2D).
//
// Status: see examples/cavity-tc/RUN_REPORT.md. OpenLB 1.9r0.
//
// Physics: square cavity L x L, top lid moving at U_lid, other walls
// fixed. The TC fluid yields near the lid (tau ~ 93 Pa > sigma_y = 70 Pa
// at the characteristic rate U_lid/L = 10 1/s) while the slow bottom
// corners stay unyielded — the video shows the yielded region vs the
// plug. Same TC rheology and eta-cap regularization as
// examples/couette-tc/ (same tc::ThreeComponent model, same
// visco::OmegaFromCell hook, same gamma_dot_reg = 1e-3 1/s).
//
// Geometry / lid boundary follow OpenLB's own laminar/cavity2d example
// (IndicatorCuboid2D lid, boundary::InterpolatedVelocity on walls and
// lid, lid velocity set with momenta::setVelocity).
//
// Validation: no analytic TC cavity solution exists. The TC dynamics is
// the same code validated in Couette-TC (1.14% max torque error), and
// the cavity setup itself is validated by a Newtonian run of this same
// binary (--newtonian, eta = 1.0 Pa*s, Re = 100) against Ghia, Ghia &
// Shin (1982) centreline velocities (see cavity_post.py).
//
// Build (serial, no MPI):
//   cd examples/cavity-tc && make
// Run:
//   ./cavity --out runs/tc                    (TC, case.yaml defaults)
//   ./cavity --newtonian --eta 1.0 --out runs/newtonian   (Re = 100 check)

#include <olb.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::OMEGA>>>;

// ---- TC rheology (identical to examples/couette-tc/couette.cpp) ----
namespace tc {

struct SIGMA_Y : public descriptors::FIELD_BASE<1> {};
struct GAMMA_C : public descriptors::FIELD_BASE<1> {};
struct ETA_BG  : public descriptors::FIELD_BASE<1> {};
struct ETA_MAX : public descriptors::FIELD_BASE<1> {};

struct ThreeComponent {
  using parameters = meta::list<SIGMA_Y, GAMMA_C, ETA_BG, ETA_MAX>;

  static std::string getName()
  {
    return "tc::ThreeComponent";
  }

  template <typename T, typename PARAMETERS>
  static T computeViscosity(PARAMETERS& parameters, T gamma) any_platform
  {
    const T sigmaY = parameters.template get<SIGMA_Y>();
    const T gammaC = parameters.template get<GAMMA_C>();
    const T etaBg  = parameters.template get<ETA_BG>();
    const T etaMax = parameters.template get<ETA_MAX>();
    if (gamma <= T(1e-12)) return etaMax; // at rest / unyielded: capped
    const T eta = sigmaY / gamma + sigmaY / util::sqrt(gamma * gammaC) + etaBg;
    return util::min(eta, etaMax);
  }
};

} // namespace tc

template <typename T, typename DESCRIPTOR>
using TcBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::BulkTuple,
  equilibria::SecondOrder,
  visco::OmegaFromCell<collision::BGK, tc::ThreeComponent>
>;

// ---- case parameters (case.yaml) ----
static double Lphys = 0.100;  // m, cavity side
static double ULid = 1.0;     // m/s
static double rhoPhys = 1000.0;
static double sigmaYPhys = 70.0;   // Pa
static double gammaCPhys = 100.0;  // 1/s
static double etaBgPhys = 0.05;    // Pa.s
static double gammaRegPhys = 1e-3; // 1/s
static double etaMaxPhys = 70221.4; // Pa.s, eta(gamma_reg)
static double runTime = 5.0;   // s
static double snapInterval = 0.25; // s between field dumps
static int N = 128;            // cells per side
static double omegaTarget = 1.0; // CLI --omega overrides
static bool newtonian = false;
static double etaNewt = 1.0;   // Pa.s, Newtonian validation run only
static std::string outDir = "runs/tc";

// ---- derived lattice parameters (same formulas as cavity.py) ----
static double dx, dt, nuLb, etaCharPhys;
static std::size_t nSteps, snapEvery;
static int nx, ny;
static const double PI = std::acos(-1.0);

static double tcEtaPhys(double gamma)
{
  if (gamma <= 0.0) return etaMaxPhys;
  const double e = sigmaYPhys / gamma
      + sigmaYPhys / std::sqrt(gamma * gammaCPhys) + etaBgPhys;
  return std::min(e, etaMaxPhys);
}

// Local shear rate (1/s) and viscosity (Pa.s) of one fluid cell,
// recovered from its non-equilibrium stress — the same formula
// visco::OmegaFromCell uses (Boyd et al.): gamma_lb = sqrt(2) *
// pre2 * |PiNeq|, pre2 = 0.5 * invCs2 * omega_cell / rho.
template <typename CELL>
static void cellGammaEta(CELL& cell, double& gammaPhys, double& etaPhysOut)
{
  using DESCRIPTOR = typename CELL::descriptor_t;
  double rho = 0.0;
  double u[2] = {0.0, 0.0};
  cell.computeRhoU(rho, u);
  const double om = cell.template getField<descriptors::OMEGA>();
  const double pre2 = 0.5 * descriptors::invCs2<double, DESCRIPTOR>() * om / rho;
  const double gammaLb = std::sqrt(2.0 * pre2 * pre2
      * lbm<DESCRIPTOR>::computePiNeqNormSqr(cell));
  gammaPhys = gammaLb / dt;
  etaPhysOut = newtonian ? etaNewt : tcEtaPhys(gammaPhys);
}

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T = MyCase::value_t;
  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  const Vector<T, 2> extent{Lphys, Lphys};
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
  auto& parameters = myCase.getParameters();

  // cavity2d pattern: everything wall (2), then the interior fluid (1),
  // then the lid row (3).
  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 1});
  const T eps = parameters.get<parameters::PHYS_DELTA_X>();
  Vector<T, 2> extend(Lphys + 2 * eps, 2 * eps);
  Vector<T, 2> origin(-eps, Lphys - eps);
  IndicatorCuboid2D<T> lid(extend, origin);
  // Boundary voxels of material 2 become the lid (3) where fluid (1)
  // is next to them inside the lid indicator.
  geometry.rename(2, 3, 1, lid);
  geometry.clean();
  geometry.innerClean();
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
      (T)Lphys,              // charPhysLength: cavity side
      (T)charU,              // charPhysVelocity: lid velocity
      (T)(etaCharPhys / rhoPhys), // physViscosity: nu at the char rate
      (T)rhoPhys             // physDensity
  );
  const auto& converter = lattice.getUnitConverter();
  converter.print();
  const T omegaConv = converter.getLatticeRelaxationFrequency();
  clout << "converter omega = " << omegaConv << " (target " << omegaTarget << ")"
        << std::endl;

  if (newtonian) {
    // Validation run: constant-viscosity BGK.
    dynamics::set<BGKdynamics>(lattice, geometry.getMaterialIndicator(1));
  } else {
    // TC dynamics — local omega from the TC law (see header comment).
    dynamics::set<TcBGKdynamics<T, DESCRIPTOR>>(
        lattice, geometry.getMaterialIndicator(1));
    lattice.setParameter<tc::SIGMA_Y>(sigmaYPhys * dt * dt / (rhoPhys * dx * dx));
    lattice.setParameter<tc::GAMMA_C>(gammaCPhys * dt);
    lattice.setParameter<tc::ETA_BG>((etaBgPhys / rhoPhys) * dt / (dx * dx));
    lattice.setParameter<tc::ETA_MAX>((etaMaxPhys / rhoPhys) * dt / (dx * dx));
  }
  fields::set<descriptors::OMEGA>(lattice,
                                  geometry.getMaterialIndicator(1),
                                  omegaConv);
  // Walls (2) and lid (3): interpolated velocity boundary (cavity2d).
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 2);
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);

  lattice.setParameter<descriptors::OMEGA>(omegaConv);
  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "setInitialValues");
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  // Lid velocity (physical units — momenta::setVelocity scales to
  // lattice units internally, src/case/setters.h; cavity2d pattern).
  // Walls keep velocity zero by default. Lid starts at full speed
  // at t = 0, as in cavity2d (no ramp).
  AnalyticalConst2D<double, double> uTop(ULid, 0.0);
  momenta::setVelocity(lattice, geometry.getMaterialIndicator(3), uTop);
  lattice.initialize();
  clout << "Set Initial Values ... OK" << std::endl;
}

void simulate(MyCase& myCase)
{
  OstreamManager clout(std::cout, "simulate");
  using T = MyCase::value_t;
  auto& lattice   = myCase.getLattice(NavierStokes{});
  auto& geometry  = myCase.getGeometry();

  std::filesystem::create_directories(outDir);
  std::ofstream energyFile(outDir + "/energy.csv");
  // mean_rho tracks mass conservation (cf. the XLB halfway
  // bounce-back leak at high rate — measured here, not assumed).
  energyFile << "step,time_s,mean_ke_m2_s2,max_speed_m_s,mean_rho\n";

  // meta.json: everything the post/video scripts need.
  {
    std::ofstream meta(outDir + "/meta.json");
    meta << "{\n  \"L_m\": " << Lphys << ",\n  \"U_lid_m_s\": " << ULid
         << ",\n  \"N\": " << N << ",\n  \"dx_m\": " << dx
         << ",\n  \"dt_s\": " << dt << ",\n  \"n_steps\": " << nSteps
         << ",\n  \"run_time_s\": " << runTime
         << ",\n  \"snap_interval_s\": " << snapInterval
         << ",\n  \"snap_every_steps\": " << snapEvery
         << ",\n  \"grid_n\": " << nx
         << ",\n  \"newtonian\": " << (newtonian ? "true" : "false")
         << ",\n  \"eta_newt_Pa_s\": " << etaNewt
         << ",\n  \"sigma_y_Pa\": " << sigmaYPhys
         << ",\n  \"gamma_dot_c_1_s\": " << gammaCPhys
         << ",\n  \"eta_bg_Pa_s\": " << etaBgPhys
         << ",\n  \"eta_max_Pa_s\": " << etaMaxPhys
         << ",\n  \"eta_char_Pa_s\": " << etaCharPhys << "\n}\n";
  }

  auto writeField = [&](int fieldIdx) {
    std::ofstream f(outDir + "/field_" + std::to_string(fieldIdx) + ".csv");
    f << "iX,iY,ux_m_s,uy_m_s,gamma_1_s,eta_Pa_s\n";
    auto& block = lattice.getBlock(0);
    auto& blockGeo = geometry.getBlockGeometry(0);
    for (int iY = 0; iY < ny; ++iY) {
      for (int iX = 0; iX < nx; ++iX) {
        if (blockGeo.get({iX, iY}) != 1) continue;
        double rho = 0.0;
        double u[2] = {0.0, 0.0};
        auto cell = block.get(iX, iY); // Cell is a lightweight view (by value)
        cell.computeRhoU(rho, u);
        double gammaPhys = 0.0, etaPhysOut = 0.0;
        cellGammaEta(cell, gammaPhys, etaPhysOut);
        f << iX << "," << iY << "," << u[0] * dx / dt << ","
          << u[1] * dx / dt << "," << gammaPhys << "," << etaPhysOut << "\n";
      }
    }
  };

  const auto tStart = std::chrono::steady_clock::now();
  int fieldIdx = 0;
  writeField(fieldIdx++);
  double keFirstHalf = 0.0, keSecondHalf = 0.0;
  std::size_t keCntFirst = 0, keCntSecond = 0;

  for (std::size_t iT = 0; iT < nSteps; ++iT) {
    lattice.collideAndStream();
    if ((iT + 1) % 1000 == 0 || iT + 1 == nSteps) {
      // Mean kinetic energy + max speed over fluid cells.
      auto& block = lattice.getBlock(0);
      auto& blockGeo = geometry.getBlockGeometry(0);
      double keSum = 0.0, maxSpeed = 0.0, rhoSum = 0.0;
      std::size_t nFluid = 0;
      for (int iY = 0; iY < ny; ++iY) {
        for (int iX = 0; iX < nx; ++iX) {
          if (blockGeo.get({iX, iY}) != 1) continue;
          double rho = 0.0;
          double u[2] = {0.0, 0.0};
          block.get(iX, iY).computeRhoU(rho, u);
          const double usq = u[0] * u[0] + u[1] * u[1];
          keSum += 0.5 * usq;
          rhoSum += rho;
          maxSpeed = std::max(maxSpeed, std::sqrt(usq));
          ++nFluid;
        }
      }
      const double keMean = keSum / nFluid * (dx / dt) * (dx / dt);
      const double maxSpeedPhys = maxSpeed * dx / dt;
      energyFile << (iT + 1) << "," << (iT + 1) * dt << "," << keMean << ","
                 << maxSpeedPhys << ","
                 << (nFluid ? rhoSum / nFluid : 0.0) << "\n";
      // Convergence bookkeeping: last 10% vs previous 10% of the run.
      if (iT + 1 > nSteps * 0.9) { keSecondHalf += keMean; ++keCntSecond; }
      else if (iT + 1 > nSteps * 0.8) { keFirstHalf += keMean; ++keCntFirst; }
    }
    if ((iT + 1) % snapEvery == 0) writeField(fieldIdx++);
    if ((iT + 1) % 50000 == 0) {
      clout << "step " << (iT + 1) << "/" << nSteps << " t=" << (iT + 1) * dt
            << " s" << std::endl;
    }
  }
  writeField(fieldIdx++);

  // Centreline u_x at x = L/2 (all rows) — Ghia comparison for the
  // Newtonian validation run; also saved for the TC run.
  {
    std::ofstream cl(outDir + "/centerline.csv");
    cl << "y_m,ux_m_s\n";
    auto& block = lattice.getBlock(0);
    const int iXc = nx / 2;
    for (int iY = 0; iY < ny; ++iY) {
      double rho = 0.0;
      double u[2] = {0.0, 0.0};
      block.get(iXc, iY).computeRhoU(rho, u);
      cl << (iY + 0.5) * dx << "," << u[0] * dx / dt << "\n";
    }
  }

  const double ke1 = keCntFirst ? keFirstHalf / keCntFirst : 0.0;
  const double ke2 = keCntSecond ? keSecondHalf / keCntSecond : 0.0;
  const double relChange = ke1 > 0 ? std::abs(ke2 - ke1) / ke1 : 0.0;
  const auto tEnd = std::chrono::steady_clock::now();
  const double wall = std::chrono::duration<double>(tEnd - tStart).count();
  const double mlups = (double)nSteps * nx * ny / wall / 1e6;
  clout << "done: " << nSteps << " steps in " << wall << " s wall ("
        << mlups << " MLUPS, " << nx << "x" << ny << " grid)" << std::endl;
  clout << "mean KE, 80-90% window: " << ke1 << " m2/s2; 90-100%: " << ke2
        << "; rel change " << relChange * 100 << " %" << std::endl;
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--out" && i + 1 < argc) outDir = argv[++i];
    else if (arg == "--newtonian") newtonian = true;
    else if (arg == "--eta" && i + 1 < argc) etaNewt = atof(argv[++i]);
    else if (arg == "--time" && i + 1 < argc) runTime = atof(argv[++i]);
    else if (arg == "--N" && i + 1 < argc) N = atoi(argv[++i]);
    else if (arg == "--ulid" && i + 1 < argc) ULid = atof(argv[++i]);
    else if (arg == "--omega" && i + 1 < argc) omegaTarget = atof(argv[++i]);
  }

  // Lattice parameters (identical to cavity.py).
  etaMaxPhys = sigmaYPhys / gammaRegPhys
      + sigmaYPhys / std::sqrt(gammaRegPhys * gammaCPhys) + etaBgPhys;
  const double gammaChar = ULid / Lphys;
  etaCharPhys = newtonian ? etaNewt : tcEtaPhys(gammaChar);
  dx = Lphys / N;
  const double nuChar = etaCharPhys / rhoPhys;
  nuLb = (1.0 / omegaTarget - 0.5) / 3.0;
  dt = nuLb * dx * dx / nuChar;
  nSteps = (std::size_t)std::ceil(runTime / dt);
  snapEvery = (std::size_t)std::max<std::size_t>(
      1, (std::size_t)std::llround(snapInterval / dt));
  nx = N + 1; // mesh points per side: N intervals + both endpoints
  ny = N + 1;

  clout << "Cavity" << (newtonian ? " (Newtonian validation)" : "-TC")
        << " (OpenLB): L=" << Lphys << " m, U_lid=" << ULid << " m/s, N=" << N
        << ", dx=" << dx << " m, dt=" << dt << " s, omega=" << omegaTarget
        << ", steps=" << nSteps << ", grid=" << nx << "x" << ny << std::endl;
  clout << (newtonian ? "Newtonian eta=" : "TC eta_char=eta(U/L)=")
        << etaCharPhys << " Pa.s, eta_max=" << etaMaxPhys
        << " Pa.s, runTime=" << runTime << " s" << std::endl;

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    using T = MyCase::value_t;
    myCaseParameters.set<RESOLUTION>(N);
    myCaseParameters.set<PHYS_DELTA_X>(dx);
    myCaseParameters.set<PHYS_MAX_T>(runTime);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(Lphys);
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(ULid);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(etaCharPhys / rhoPhys);
    myCaseParameters.set<PHYS_DENSITY>(rhoPhys);
    myCaseParameters.set<OVERLAP>(2); // lid rename needs overlap>=2
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
