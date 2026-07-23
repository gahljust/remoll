// metal_host — GPU transport driver (Apple Metal, Obj-C++). Generates rate
// nodes with the SAME rate_solver_core the CPU pipeline uses, pushes them
// through transport.metal (one thread per trajectory), and reduces per-plane
// accepted rates from the returned bitmasks.
//
//   gpu/build_gpu.sh                       # on the Mac (needs Xcode CLT)
//   gpu/transport_gpu --channel moller ... --field m1 --field m2 \
//       --planes ../transport/planes_full.csv [--cpu-check N] [rate flags]
//
// --cpu-check N: additionally runs the FIRST N nodes through the
// double-precision transport_core on the CPU and prints both totals — the
// on-Mac verification that GPU float32 + the ported algorithm reproduce the
// reference (expect sub-percent differences from float and edge nodes; see
// kernel_ref.cc for the sandbox-side algorithm equivalence test).
//
// v1 limits (see ../DIRECTION.md): accepted-rate totals only (no r/phi/xy
// maps, no eout classes on GPU — trivial to add via extra output fields);
// masks are uploaded whole (atlas ~ a few MB); nodes stream in 2M batches so
// a 16 GB M4 never holds more than ~200 MB of GPU buffers.

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "../rate_solver/rate_solver_core.hh"
#include "../transport/field_map.hh"
#include "../transport/transport_core.hh"
#include "kernel_types.h"

#include <cstdio>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace vacuum_transport;

// Fixed-memory weighted quantiles for compact accepted-kinematics summaries.
// 4096 bins are much finer than the source quadrature and avoid retaining or
// writing millions of accepted-node records merely to find a biasing window.
struct WeightedHistogram {
  double lo = 0.0, hi = 1.0;
  std::vector<double> bins;
  WeightedHistogram() : bins(4096, 0.0) {}
  WeightedHistogram(double a, double b) : lo(a), hi(b), bins(4096, 0.0) {}
  void reset() { std::fill(bins.begin(), bins.end(), 0.0); }
  void add(double x, double weight) {
    if (!(weight > 0.0) || !std::isfinite(x) || !(hi > lo)) return;
    const double u = std::max(0.0, std::min(1.0, (x - lo) / (hi - lo)));
    const size_t i = std::min(bins.size() - 1,
                              size_t(u * double(bins.size())));
    bins[i] += weight;
  }
  double quantile(double q) const {
    double total = 0.0;
    for (double w : bins) total += w;
    if (!(total > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    const double target = std::max(0.0, std::min(1.0, q)) * total;
    double cumulative = 0.0;
    for (size_t i = 0; i < bins.size(); ++i) {
      cumulative += bins[i];
      if (cumulative >= target)
        return lo + (double(i) + 0.5) * (hi - lo) / double(bins.size());
    }
    return hi;
  }
};

static std::vector<PlaneDef> load_planes_csv(const std::string& path) {
  std::ifstream in(path);
  if (!in) { std::fprintf(stderr, "cannot open %s\n", path.c_str()); std::exit(1); }
  std::vector<PlaneDef> planes; std::string line;
  std::getline(in, line);
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> c; std::stringstream ss(line); std::string cell;
    while (std::getline(ss, cell, ',')) c.push_back(cell);
    if (c.size() < 8) continue;
    PlaneDef p; p.name = c[0]; p.z_m = std::stod(c[1]);
    p.r_min_m = std::stod(c[2]); p.r_max_m = std::stod(c[3]);
    p.x_min_m = std::stod(c[4]); p.x_max_m = std::stod(c[5]);
    p.y_min_m = std::stod(c[6]); p.y_max_m = std::stod(c[7]);
    if (c.size() >= 9)  p.is_aperture = (std::stoi(c[8]) != 0);
    if (c.size() >= 10) p.n_septant = std::stoi(c[9]);
    if (c.size() >= 11) p.sept_phi0_deg = std::stod(c[10]);
    if (c.size() >= 12) p.sept_half_deg = std::stod(c[11]);
    if (c.size() >= 13 && !c[12].empty()) p.mask = load_mask(c[12]);
    if (c.size() >= 14 && !c[13].empty()) p.invert = (std::stoi(c[13]) != 0);
    planes.push_back(p);
  }
  std::sort(planes.begin(), planes.end(),
            [](const PlaneDef& a, const PlaneDef& b) { return a.z_m < b.z_m; });
  return planes;
}

int main(int argc, char** argv) { @autoreleasepool {
  solver_model::RunConfig cfg;
  cfg.channel = solver_model::kMoller;
  std::vector<std::string> field_paths;
  std::string planes_path, output = "transport_gpu";
  long long cpu_check = 0;
  double p_min_gev = 0.0;
  std::vector<double> eout_edges;   // --eout-classes, as transport_solver
  // --auto-phi TOL: rerun with phi_nodes doubled until the max relative
  // change of any plane carrying >0.1% of the generated rate is < TOL.
  // The last doubling's per-plane |delta| is reported as the phi error.
  // Cost: sum of levels ~ 2x the final level — cheap at GPU speed, and it is
  // the honest estimator: interleaved-half splits were measured to UNDERstate
  // the true error ~8x (edge correlation), level deltas bracket it from above.
  double auto_phi_tol = 0.0;        // 0 = off (single pass)
  int auto_phi_max = 160;           // per-sector cap
  // --record-impacts NAME[,NAME...]: for each listed APERTURE plane, write
  // <output>.impacts_<NAME>.tsv (one row per node whose trajectory ended on
  // that structure: impact position/direction/momentum + origin kinematics +
  // absolute rate weight) and <output>.impacts_<NAME>.json (stratum manifest:
  // total intercepted rate R_h, the Mode-A kinematic bounding window of the
  // contributing nodes, and the window's box efficiency). This is the input
  // to the secondaries campaign (see ../injection/ and ../TESTPLAN_SECONDARIES.md).
  std::vector<std::string> impact_names;
  std::string accepted_name;       // accepted source-state dump at one monitor
  std::string kinematics_name;     // compact source-kinematics summary plane
  std::string kinematics_output;
  bool write_json = true;
  TransportConfig tref;   // reference CPU config (also source of GPU params)

  for (int i = 1; i < argc; ++i) {
    const std::string o = argv[i];
    auto next = [&]() -> const char* {
      if (++i >= argc) { std::fprintf(stderr, "missing value for %s\n", o.c_str()); std::exit(2); }
      return argv[i];
    };
    if (o == "--channel") cfg.channel = solver_model::channel_from_name(next());
    else if (o == "--beam-mev") cfg.beam_total_mev = atof(next());
    else if (o == "--current-uA") cfg.beam_current_ampere = atof(next()) * 1e-6;
    else if (o == "--target-Z") cfg.target_Z = atof(next());
    else if (o == "--target-A") cfg.target_A = atof(next());
    else if (o == "--thickness-mm") cfg.target_thickness_mm = atof(next());
    else if (o == "--density-g-cm3") cfg.target_density_g_cm3 = atof(next());
    else if (o == "--molar-mass") cfg.target_molar_mass_g_mol = atof(next());
    else if (o == "--radiation-length-mm") cfg.target_radiation_length_mm = atof(next());
    else if (o == "--target-z0-mm") cfg.target_z0_mm = atof(next());
    else if (o == "--target-z-span-mm") cfg.target_z_span_mm = atof(next());
    else if (o == "--theta-min-deg") cfg.theta_min_deg = atof(next());
    else if (o == "--theta-max-deg") cfg.theta_max_deg = atof(next());
    else if (o == "--vertex-nodes") cfg.vertex_nodes = atoi(next());
    else if (o == "--energy-nodes") cfg.energy_nodes = atoi(next());
    else if (o == "--theta-nodes") cfg.theta_nodes = atoi(next());
    else if (o == "--phi-nodes") cfg.phi_nodes = atoi(next());
    else if (o == "--phi-fold") cfg.phi_fold = atoi(next());
    else if (o == "--screening") cfg.screening = true;
    else if (o == "--no-msc") cfg.msc_auto = false;
    else if (o == "--ep-internal-brems") cfg.ep_internal_brems = true;
    else if (o == "--eout-floor-mev") cfg.eout_floor_mev = atof(next());
    else if (o == "--postvertex-material") cfg.postvertex_material = true;
    else if (o == "--exit-window-al-mm") cfg.exit_window_al_mm = atof(next());
    else if (o == "--field") field_paths.push_back(next());
    else if (o == "--planes") planes_path = next();
    else if (o == "--ds-field-m") tref.ds_field_m = atof(next());
    else if (o == "--ds-drift-m") tref.ds_drift_m = atof(next());
    else if (o == "--b-min-tesla") tref.b_min_tesla = atof(next());
    else if (o == "--z-stop-m") tref.z_stop_m = atof(next());
    else if (o == "--r-stop-m") tref.r_stop_m = atof(next());
    else if (o == "--max-path-m") tref.max_path_m = atof(next());
    else if (o == "--max-steps") tref.max_steps = atoi(next());
    else if (o == "--p-min-gev") p_min_gev = atof(next());
    else if (o == "--cpu-check") cpu_check = atoll(next());
    else if (o == "--auto-phi") auto_phi_tol = atof(next());
    else if (o == "--record-impacts") {
      std::stringstream ss(next()); std::string tok;
      while (std::getline(ss, tok, ',')) impact_names.push_back(tok);
    }
    else if (o == "--record-accepted") accepted_name = next();
    else if (o == "--kinematics-summary") kinematics_name = next();
    else if (o == "--kinematics-output") kinematics_output = next();
    else if (o == "--no-json") write_json = false;
    else if (o == "--auto-phi-max") auto_phi_max = atoi(next());
    else if (o == "--eout-classes") {
      std::stringstream ss(next()); std::string tok;
      while (std::getline(ss, tok, ',')) eout_edges.push_back(std::stod(tok));
      std::sort(eout_edges.begin(), eout_edges.end());
    }
    else if (o == "--output") output = next();
    else { std::fprintf(stderr, "unknown option %s\n", o.c_str()); return 2; }
  }
  if (field_paths.empty() || planes_path.empty()) {
    std::fprintf(stderr, "need --field ... --planes csv (+ rate flags)\n"); return 2;
  }

  // ---- node generation (identical path to the CPU pipeline) ----
  const int ncl = int(eout_edges.size()) + 1;
  auto ecl_of = [&](double e_mev) {
    int k = 0;
    for (double edge : eout_edges) { if (e_mev >= edge) ++k; else break; }
    return k;
  };
  struct N { GpuNodeIn in; double rate; int ecl; int par; int phi_index;
             float th_deg, th_lab_deg, ein_mev, zv_mm, phi_deg, eout_mev;
             float x_mm, y_mm, q2_mev2, w_mev; };
  std::vector<N> nodes;
  double total_rate = 0.0, skipped = 0.0;
  auto generate_nodes = [&]() {
    nodes.clear(); total_rate = 0.0; skipped = 0.0;
    solver_model::for_each_node(cfg, [&](const solver_model::NodeOut& m) {
      constexpr double me = 0.510998950;
      total_rate += m.rate;
      if (m.theta_lab_deg > 90.0 || m.eout_mev <= me) { skipped += m.rate; return; }
      const double p = std::sqrt(m.eout_mev * m.eout_mev - me * me) / 1000.0;
      if (p < p_min_gev) { skipped += m.rate; return; }
      N n;
      n.in.rx = float(m.x_mm / 1000.0); n.in.ry = float(m.y_mm / 1000.0);
      n.in.rz = float(m.z_mm / 1000.0);
      n.in.tx = float(m.tx); n.in.ty = float(m.ty);
      n.in.p_gev = float(p);
      n.rate = m.rate;
      n.ecl = ecl_of(m.eout_mev);
      n.par = m.phi_index >= 0 ? (m.phi_index & 1) : -1;
      n.phi_index = m.phi_index;
      n.th_deg = float(m.theta_deg);
      n.th_lab_deg = float(m.theta_lab_deg);
      n.ein_mev = float(m.ein_mev);
      n.zv_mm = float(m.z_mm);
      n.phi_deg = float(std::atan2(m.ty, m.tx) * 180.0 / kPi);
      n.eout_mev = float(m.eout_mev);
      n.x_mm = float(m.x_mm); n.y_mm = float(m.y_mm);
      n.q2_mev2 = float(m.q2_mev2); n.w_mev = float(m.w_mev);
      nodes.push_back(n);
    });
    std::printf("phi_nodes %d: nodes %zu  total rate %.6e /s (skipped %.3e kept)\n",
                cfg.phi_nodes, nodes.size(), total_rate, skipped);
  };
  generate_nodes();

  // ---- geometry/fields -> GPU layouts ----
  auto planes = load_planes_csv(planes_path);
  if ((int)planes.size() > GPU_MAX_PLANES) { std::fprintf(stderr, "too many planes\n"); return 1; }
  int accepted_plane = -1;
  const std::string selected_name = !kinematics_name.empty()
      ? kinematics_name : accepted_name;
  if (!selected_name.empty()) {
    for (size_t j = 0; j < planes.size(); ++j)
      if (planes[j].name == selected_name) accepted_plane = int(j);
    if (accepted_plane < 0) {
      std::fprintf(stderr, "accepted-kinematics plane %s not found\n",
                   selected_name.c_str());
      return 2;
    }
    if (!accepted_name.empty() && !kinematics_name.empty()
        && accepted_name != kinematics_name) {
      std::fprintf(stderr, "--record-accepted and --kinematics-summary must use the same plane\n");
      return 2;
    }
  }
  FieldSet fset;
  for (auto& f : field_paths) fset.add(f);

  std::vector<GpuFieldHead> heads;
  std::vector<float> bflat;   // packed float3
  for (const auto& m : fset.maps) {
    GpuFieldHead h{};
    h.nr = m.nr; h.nphi = m.nphi; h.nz = m.nz;
    h.rmin = m.rmin; h.rmax = m.rmax;
    h.phimin = m.phimin; h.phimax = m.phimax;
    h.zmin = m.zmin + m.z_offset + m.z_map_offset;   // fold offsets into grid
    h.zmax = m.zmax + m.z_offset + m.z_map_offset;
    h.phi_map_offset = m.phi_map_offset;
    h.phi_low = m.phi_low; h.xtant_size = m.xtant_size; h.nxtant = m.nxtant;
    h.scale = m.scale;
    h.data_offset = int(bflat.size() / 3);
    heads.push_back(h);
    for (const auto& b : m.cyl_b) {
      bflat.push_back(float(b.x)); bflat.push_back(float(b.y)); bflat.push_back(float(b.z));
    }
  }

  std::vector<GpuPlane> gpl;
  std::vector<unsigned char> atlas;
  for (const auto& p : planes) {
    GpuPlane g{};
    g.z_m = p.z_m;
    g.r_min = p.r_min_m; g.r_max = float(std::min(p.r_max_m, 1e6));
    g.x_min = float(std::max(p.x_min_m, -1e6)); g.x_max = float(std::min(p.x_max_m, 1e6));
    g.y_min = float(std::max(p.y_min_m, -1e6)); g.y_max = float(std::min(p.y_max_m, 1e6));
    g.dx = p.dx_m; g.dy = p.dy_m;
    g.n_septant = p.n_septant;
    g.sept_phi0_rad = p.sept_phi0_deg * kPi / 180.0;
    g.sept_half_rad = p.sept_half_deg * kPi / 180.0;
    g.flags = (p.is_aperture ? 1 : 0) | (p.invert ? 2 : 0);
    g.mask_offset = -1;
    if (p.mask.loaded) {
      g.mask_offset = int(atlas.size());
      g.mask_nr = p.mask.nr; g.mask_nphi = p.mask.nphi; g.mask_rmax = p.mask.rmax_m;
      atlas.insert(atlas.end(), p.mask.open.begin(), p.mask.open.end());
    }
    gpl.push_back(g);
  }
  if (atlas.empty()) atlas.push_back(0);   // Metal buffers cannot be zero-length

  GpuParams prm{};
  prm.charge = tref.charge;
  prm.ds_field = tref.ds_field_m; prm.ds_drift = tref.ds_drift_m;
  prm.b_min = tref.b_min_tesla;
  prm.z_stop = tref.z_stop_m; prm.r_stop = tref.r_stop_m;
  prm.max_path = tref.max_path_m; prm.max_steps = tref.max_steps;
  prm.n_planes = int(planes.size()); prm.n_fields = int(heads.size());
  prm.record_plane = accepted_plane;

  // ---- Metal setup (compile the kernel from source at runtime) ----
  id<MTLDevice> dev = MTLCreateSystemDefaultDevice();
  if (!dev) { std::fprintf(stderr, "no Metal device\n"); return 1; }
  std::printf("Metal device: %s\n", dev.name.UTF8String);
  NSString* dir = [[NSString stringWithUTF8String:argv[0]]
                   stringByDeletingLastPathComponent];
  NSString* src = [NSString stringWithContentsOfFile:
                   [dir stringByAppendingPathComponent:@"transport.metal"]
                   encoding:NSUTF8StringEncoding error:nil];
  if (!src) { std::fprintf(stderr, "transport.metal not found next to binary\n"); return 1; }
  // kernel_types.h is #included by the shader; inline it for runtime compile
  NSString* hdr = [NSString stringWithContentsOfFile:
                   [dir stringByAppendingPathComponent:@"kernel_types.h"]
                   encoding:NSUTF8StringEncoding error:nil];
  src = [src stringByReplacingOccurrencesOfString:@"#include \"kernel_types.h\""
                                       withString:hdr];
  NSError* err = nil;
  id<MTLLibrary> lib = [dev newLibraryWithSource:src options:nil error:&err];
  if (!lib) { std::fprintf(stderr, "Metal compile failed:\n%s\n",
                           err.localizedDescription.UTF8String); return 1; }
  id<MTLFunction> fn = [lib newFunctionWithName:@"transport_nodes"];
  id<MTLComputePipelineState> pso = [dev newComputePipelineStateWithFunction:fn error:&err];
  if (!pso) { std::fprintf(stderr, "pipeline failed\n"); return 1; }
  id<MTLCommandQueue> q = [dev newCommandQueue];

  auto mkbuf = [&](const void* ptr, size_t bytes) {
    return [dev newBufferWithBytes:ptr length:std::max(bytes, size_t(16))
                           options:MTLResourceStorageModeShared];
  };
  id<MTLBuffer> bPrm   = mkbuf(&prm, sizeof(prm));
  id<MTLBuffer> bHead  = mkbuf(heads.data(), heads.size() * sizeof(GpuFieldHead));
  id<MTLBuffer> bField = mkbuf(bflat.data(), bflat.size() * sizeof(float));
  id<MTLBuffer> bPl    = mkbuf(gpl.data(), gpl.size() * sizeof(GpuPlane));
  id<MTLBuffer> bAtlas = mkbuf(atlas.data(), atlas.size());

  // ---- dispatch in batches, reduce on CPU (one quadrature level) ----
  const size_t BATCH = 2'000'000;
  std::vector<double> accepted(planes.size(), 0.0), reached(planes.size(), 0.0);
  std::vector<double> accepted_cl(planes.size() * size_t(ncl), 0.0);
  std::vector<double> acc_even(planes.size(), 0.0), acc_odd(planes.size(), 0.0);
  std::vector<GpuNodeOut> firstN;
  const std::string accepted_path = accepted_name.empty()
      ? std::string() : output + ".accepted_" + accepted_name + ".tsv";
  const double phi_span_deg = 360.0 / std::max(1, cfg.phi_fold);
  WeightedHistogram kin_theta(cfg.theta_min_deg, cfg.theta_max_deg);
  WeightedHistogram kin_phi(0.0, phi_span_deg);
  WeightedHistogram kin_beamp(0.510998950 / 1000.0,
                              cfg.beam_total_mev / 1000.0);
  WeightedHistogram kin_vertex(0.0, 1.0);
  WeightedHistogram kin_outgoing(0.0, 1.0);
  size_t kin_nodes = 0;
  double kin_sumw = 0.0, kin_sumw2 = 0.0;
  // impact recording: plane index -> rows collected on the CURRENT level
  std::vector<int> impact_idx;               // planes.size(), -1 = not recorded
  struct ImpRow { float x, y, nx, ny, nz; size_t node; };
  std::vector<std::vector<ImpRow>> impacts;  // per recorded structure
  impact_idx.assign(planes.size(), -1);
  for (size_t w = 0; w < impact_names.size(); ++w) {
    bool found = false;
    for (size_t j = 0; j < planes.size(); ++j)
      if (planes[j].name == impact_names[w]) {
        if (!planes[j].is_aperture)
          { std::fprintf(stderr, "--record-impacts %s: not an aperture\n",
                         impact_names[w].c_str()); std::exit(2); }
        impact_idx[j] = int(w); found = true;
      }
    if (!found) { std::fprintf(stderr, "--record-impacts: plane %s not found\n",
                               impact_names[w].c_str()); std::exit(2); }
  }
  impacts.assign(impact_names.size(), {});
  double dt = 0.0;
  auto run_level = [&]() {
  for (auto& v : impacts) v.clear();   // keep only the final level's impacts
  std::fill(accepted.begin(), accepted.end(), 0.0);
  std::fill(reached.begin(), reached.end(), 0.0);
  std::fill(accepted_cl.begin(), accepted_cl.end(), 0.0);
  std::fill(acc_even.begin(), acc_even.end(), 0.0);
  std::fill(acc_odd.begin(), acc_odd.end(), 0.0);
  kin_theta.reset(); kin_phi.reset(); kin_beamp.reset();
  kin_vertex.reset(); kin_outgoing.reset();
  kin_nodes = 0; kin_sumw = 0.0; kin_sumw2 = 0.0;
  firstN.assign(std::min<size_t>(cpu_check, nodes.size()), GpuNodeOut{});
  std::ofstream accepted_file;
  if (!accepted_path.empty()) {
    accepted_file.open(accepted_path);
    accepted_file << "# accepted source states at " << accepted_name << "\n"
                  << "# source variables are weighted quadrature nodes; weight_hz is physical rate\n"
                  << "# x_mm\ty_mm\tz_mm\tphi_deg\ttheta_deg\ttheta_lab_deg\tein_mev\teout_mev\tvertex_z_mm\tq2_mev2\tw_mev\tcol2_x_mm\tcol2_y_mm\tcol2_phi_deg\tweight_hz\n";
  }
  NSDate* t0 = [NSDate date];
  for (size_t off = 0; off < nodes.size(); off += BATCH) {
    const size_t n = std::min(BATCH, nodes.size() - off);
    std::vector<GpuNodeIn> in(n);
    for (size_t i = 0; i < n; ++i) in[i] = nodes[off + i].in;
    id<MTLBuffer> bIn  = mkbuf(in.data(), n * sizeof(GpuNodeIn));
    id<MTLBuffer> bOut = [dev newBufferWithLength:n * sizeof(GpuNodeOut)
                                          options:MTLResourceStorageModeShared];
    id<MTLCommandBuffer> cb = [q commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
    [enc setComputePipelineState:pso];
    [enc setBuffer:bIn    offset:0 atIndex:0];
    [enc setBuffer:bOut   offset:0 atIndex:1];
    [enc setBuffer:bPrm   offset:0 atIndex:2];
    [enc setBuffer:bHead  offset:0 atIndex:3];
    [enc setBuffer:bField offset:0 atIndex:4];
    [enc setBuffer:bPl    offset:0 atIndex:5];
    [enc setBuffer:bAtlas offset:0 atIndex:6];
    NSUInteger tg = pso.maxTotalThreadsPerThreadgroup;
    [enc dispatchThreads:MTLSizeMake(n, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min<NSUInteger>(tg, 256), 1, 1)];
    [enc endEncoding];
    [cb commit];
    [cb waitUntilCompleted];
    const GpuNodeOut* o = (const GpuNodeOut*)bOut.contents;
    for (size_t i = 0; i < n; ++i)
      if (off + i < firstN.size()) firstN[off + i] = o[i];
    for (size_t i = 0; i < n; ++i) {
      const int bp = o[i].blocked_plane;
      if (bp >= 0 && bp < (int)planes.size() && impact_idx[bp] >= 0)
        impacts[impact_idx[bp]].push_back(
            {o[i].imp_x, o[i].imp_y, o[i].imp_nx, o[i].imp_ny, o[i].imp_nz,
             off + i});
    }
    for (size_t i = 0; i < n; ++i) {
      const double r = nodes[off + i].rate;
      for (size_t j = 0; j < planes.size(); ++j) {
        const bool acc = (j < 32) ? (o[i].accepted_lo >> j) & 1u
                                  : (o[i].accepted_hi >> (j - 32)) & 1u;
        const bool rea = (j < 32) ? (o[i].reached_lo >> j) & 1u
                                  : (o[i].reached_hi >> (j - 32)) & 1u;
        if (acc) {
          accepted[j] += r;
          accepted_cl[j * size_t(ncl) + nodes[off + i].ecl] += r;
          if (nodes[off + i].par == 0) acc_even[j] += r;
          else if (nodes[off + i].par == 1) acc_odd[j] += r;
        }
        if (rea) reached[j] += r;
      }
    }
    if (accepted_plane >= 0) {
      for (size_t i = 0; i < n; ++i) {
        const bool acc = accepted_plane < 32
            ? ((o[i].accepted_lo >> accepted_plane) & 1u)
            : ((o[i].accepted_hi >> (accepted_plane - 32)) & 1u);
        if (!acc) continue;
        const N& nd = nodes[off + i];
        if (accepted_file.is_open()) {
          const double cphi = std::atan2(double(o[i].record_y), double(o[i].record_x))
                            * 180.0 / kPi;
          accepted_file << nd.x_mm << "\t" << nd.y_mm << "\t"
                        << nd.in.rz * 1000.0 << "\t" << nd.phi_deg << "\t"
                        << nd.th_deg << "\t" << nd.th_lab_deg << "\t"
                        << nd.ein_mev << "\t" << nd.eout_mev << "\t"
                        << nd.zv_mm << "\t" << nd.q2_mev2 << "\t"
                        << nd.w_mev << "\t" << o[i].record_x * 1000.0f << "\t"
                        << o[i].record_y * 1000.0f << "\t" << cphi << "\t"
                        << nd.rate << "\n";
        }
        if (!kinematics_name.empty()) {
          // phi_index, unlike the transported electron direction, is the
          // generator scattering-plane coordinate. This maps a Moller partner
          // back to the phi value remoll's /bias/phi command samples.
          const double phi_for_bias = cfg.phi_nodes > 0 && nd.phi_index >= 0
              ? phi_span_deg * (double(nd.phi_index) + 0.5) / double(cfg.phi_nodes)
              : 0.0;
          const double vertex_fraction = cfg.target_z_span_mm != 0.0
              ? (double(nd.zv_mm) - cfg.target_z0_mm) / cfg.target_z_span_mm : 0.0;
          const double outgoing_fraction = nd.ein_mev > 0.0
              ? double(nd.eout_mev) / double(nd.ein_mev) : 0.0;
          kin_theta.add(nd.th_deg, nd.rate);
          kin_phi.add(phi_for_bias, nd.rate);
          kin_beamp.add(double(nd.ein_mev) / 1000.0, nd.rate);
          kin_vertex.add(vertex_fraction, nd.rate);
          kin_outgoing.add(outgoing_fraction, nd.rate);
          ++kin_nodes; kin_sumw += nd.rate; kin_sumw2 += nd.rate * nd.rate;
        }
      }
    }
    std::printf("  batch %zu/%zu done\n", off / BATCH + 1,
                (nodes.size() + BATCH - 1) / BATCH);
  }
  dt = -[t0 timeIntervalSinceNow];
  std::printf("GPU transported %zu trajectories in %.2f s  (%.2e traj/s)\n",
              nodes.size(), dt, nodes.size() / std::max(dt, 1e-9));
  };  // run_level

  // ---- single pass, or successive phi-doubling with level-delta errors ----
  std::vector<double> phi_err(planes.size(), 0.0);
  int levels_run = 1;
  run_level();
  if (auto_phi_tol > 0.0) {
    std::vector<double> prev = accepted;
    while (cfg.phi_nodes * 2 <= auto_phi_max) {
      cfg.phi_nodes *= 2;
      generate_nodes();
      run_level();
      ++levels_run;
      double worst = 0.0;
      for (size_t j = 0; j < planes.size(); ++j) {
        phi_err[j] = std::fabs(accepted[j] - prev[j]);
        if (accepted[j] > 1e-3 * total_rate)
          worst = std::max(worst, phi_err[j] / accepted[j]);
      }
      std::printf("auto-phi: phi_nodes %d, worst relative level-delta %.3e "
                  "(tol %.3e)\n", cfg.phi_nodes, worst, auto_phi_tol);
      if (worst < auto_phi_tol) break;
      prev = accepted;
    }
  }

  std::printf("%-18s %14s %14s\n", "plane", "accepted(/s)", "acc/total");
  for (size_t j = 0; j < planes.size(); ++j)
    std::printf("%-18s %14.6e %14.6f\n", planes[j].name.c_str(), accepted[j],
                total_rate > 0 ? accepted[j] / total_rate : 0.0);

  if (!kinematics_name.empty()) {
    if (kinematics_output.empty()) kinematics_output = output + ".kinematics.tsv";
    std::ofstream ks(kinematics_output);
    ks << "plane\taccepted_nodes\taccepted_rate_per_s\tneff"
          "\ttheta_q01_deg\ttheta_q50_deg\ttheta_q99_deg"
          "\tphi_q01_deg\tphi_q50_deg\tphi_q99_deg"
          "\tbeamp_q01_gev\tbeamp_q50_gev\tbeamp_q99_gev"
          "\tvertexz_q01_fraction\tvertexz_q50_fraction\tvertexz_q99_fraction"
          "\toutgoinge_q01_fraction\toutgoinge_q50_fraction\toutgoinge_q99_fraction\n";
    ks.precision(12);
    ks << kinematics_name << "\t" << kin_nodes << "\t" << kin_sumw << "\t"
       << (kin_sumw2 > 0.0 ? kin_sumw * kin_sumw / kin_sumw2 : 0.0)
       << "\t" << kin_theta.quantile(.01) << "\t" << kin_theta.quantile(.50)
       << "\t" << kin_theta.quantile(.99)
       << "\t" << kin_phi.quantile(.01) << "\t" << kin_phi.quantile(.50)
       << "\t" << kin_phi.quantile(.99)
       << "\t" << kin_beamp.quantile(.01) << "\t" << kin_beamp.quantile(.50)
       << "\t" << kin_beamp.quantile(.99)
       << "\t" << kin_vertex.quantile(.01) << "\t" << kin_vertex.quantile(.50)
       << "\t" << kin_vertex.quantile(.99)
       << "\t" << kin_outgoing.quantile(.01) << "\t" << kin_outgoing.quantile(.50)
       << "\t" << kin_outgoing.quantile(.99) << "\n";
    std::printf("wrote compact accepted kinematics: %s\n", kinematics_output.c_str());
  }

  if (write_json) {
    std::ofstream js(output + ".json");
    js << "{\n  \"engine\": \"metal\",\n"
     << "  \"phi_nodes_final\": " << cfg.phi_nodes << ",\n"
     << "  \"auto_phi_levels\": " << levels_run << ",\n"
     << "  \"transported_nodes\": " << nodes.size() << ",\n"
     << "  \"total_rate_per_s\": " << total_rate << ",\n"
     << "  \"transported_rate_per_s\": " << (total_rate - skipped) << ",\n"
     << "  \"skipped_rate_per_s\": " << skipped << ",\n"
     << "  \"eout_class_edges_mev\": [";
  for (size_t k = 0; k < eout_edges.size(); ++k)
    js << eout_edges[k] << (k + 1 < eout_edges.size() ? ", " : "");
  js << "],\n  \"planes\": [\n";
  for (size_t j = 0; j < planes.size(); ++j) {
    js << "    {\"plane\": \"" << planes[j].name << "\", \"z_m\": " << planes[j].z_m
       << ", \"accepted_rate_per_s\": " << accepted[j]
       << ", \"acceptance_fraction\": "
       << (total_rate > 0 ? accepted[j] / total_rate : 0.0);
    if (acc_even[j] + acc_odd[j] > 0.0)
      js << ", \"phi_split_error_per_s\": "
         << 0.5 * std::fabs(acc_even[j] - acc_odd[j]);
    if (auto_phi_tol > 0.0)
      js << ", \"phi_level_delta_per_s\": " << phi_err[j];
    if (ncl > 1) {
      js << ", \"accepted_by_eout\": [";
      for (int k = 0; k < ncl; ++k)
        js << accepted_cl[j * size_t(ncl) + k] << (k + 1 < ncl ? ", " : "");
      js << "]";
    }
    js << "}" << (j + 1 < planes.size() ? "," : "") << "\n";
  }
    js << "  ]\n}\n";
    std::printf("wrote %s.json (schema-compatible with transport_solver)\n",
                output.c_str());
  }

  // ---- impact files: the secondaries-campaign handoff ----
  for (size_t w = 0; w < impact_names.size(); ++w) {
    const std::string base = output + ".impacts_" + impact_names[w];
    int pidx = -1;
    for (size_t j = 0; j < planes.size(); ++j)
      if (impact_idx[j] == (int)w) pidx = int(j);
    const double zplane_m = planes[size_t(pidx)].z_m;
    double Rh = 0.0;
    float th_lo = 1e9f, th_hi = -1e9f, ei_lo = 1e9f, ei_hi = -1e9f;
    float zv_lo = 1e9f, zv_hi = -1e9f, ph_lo = 1e9f, ph_hi = -1e9f;
    {
      std::ofstream f(base + ".tsv");
      f << "# impact states on " << impact_names[w] << " (z=" << zplane_m
        << " m). Units: position mm, momentum MeV, weight Hz (at the run's"
           " current; matrix runs are per uA).\n"
        << "# x_mm\ty_mm\tz_mm\tpx_mev\tpy_mev\tpz_mev\tweight_hz\t"
           "theta_thrown_deg\tein_mev\tvertex_z_mm\tphi_deg\teout_mev\n";
      for (const auto& r : impacts[w]) {
        const N& nd = nodes[r.node];
        const double pm = nd.in.p_gev * 1000.0;
        f << r.x * 1000.0 << "\t" << r.y * 1000.0 << "\t" << zplane_m * 1000.0
          << "\t" << pm * r.nx << "\t" << pm * r.ny << "\t" << pm * r.nz
          << "\t" << nd.rate << "\t" << nd.th_deg << "\t" << nd.ein_mev
          << "\t" << nd.zv_mm << "\t" << nd.phi_deg << "\t" << nd.eout_mev
          << "\n";
        Rh += nd.rate;
        th_lo = std::min(th_lo, nd.th_deg); th_hi = std::max(th_hi, nd.th_deg);
        ei_lo = std::min(ei_lo, nd.ein_mev); ei_hi = std::max(ei_hi, nd.ein_mev);
        zv_lo = std::min(zv_lo, nd.zv_mm);  zv_hi = std::max(zv_hi, nd.zv_mm);
        ph_lo = std::min(ph_lo, nd.phi_deg); ph_hi = std::max(ph_hi, nd.phi_deg);
      }
    }
    // Mode-A box efficiency: what fraction of the bounding window's thrown
    // rate actually ends on this structure (tells you the cost of a stock
    // remoll window run before paying for it).
    double box_rate = 0.0;
    if (!impacts[w].empty())
      for (const auto& nd : nodes)
        if (nd.th_deg >= th_lo && nd.th_deg <= th_hi &&
            nd.ein_mev >= ei_lo && nd.ein_mev <= ei_hi &&
            nd.zv_mm >= zv_lo && nd.zv_mm <= zv_hi)
          box_rate += nd.rate;
    std::ofstream mj(base + ".json");
    mj << "{\n  \"structure\": \"" << impact_names[w] << "\",\n"
       << "  \"z_m\": " << zplane_m << ",\n"
       << "  \"n_impacts\": " << impacts[w].size() << ",\n"
       << "  \"R_h_per_s\": " << Rh << ",\n"
       << "  \"share_of_generated\": " << (total_rate > 0 ? Rh / total_rate : 0.0) << ",\n"
       << "  \"modeA_window\": {\"theta_deg\": [" << th_lo << ", " << th_hi
       << "], \"ein_mev\": [" << ei_lo << ", " << ei_hi
       << "], \"vertex_z_mm\": [" << zv_lo << ", " << zv_hi
       << "], \"phi_deg\": [" << ph_lo << ", " << ph_hi << "]},\n"
       << "  \"modeA_box_rate_per_s\": " << box_rate << ",\n"
       << "  \"modeA_box_efficiency\": " << (box_rate > 0 ? Rh / box_rate : 0.0) << "\n}\n";
    std::printf("impacts on %-14s: %8zu states, R_h %.4e /s "
                "(Mode-A box eff %.3f) -> %s.tsv/.json\n",
                impact_names[w].c_str(), impacts[w].size(), Rh,
                box_rate > 0 ? Rh / box_rate : 0.0, base.c_str());
  }

  // ---- optional CPU double-precision check on the SAME first N nodes ----
  if (cpu_check > 0) {
    const size_t n = firstN.size();
    NSDate* cpu_t0 = [NSDate date];
    std::vector<double> acc_cpu(planes.size(), 0.0), acc_gpu(planes.size(), 0.0);
    long long disagree = 0;
    for (size_t i = 0; i < n; ++i) {
      const auto& nd = nodes[i];
      // GPU side from the saved bitmask
      for (size_t j = 0; j < planes.size(); ++j) {
        const bool acc = (j < 32) ? (firstN[i].accepted_lo >> j) & 1u
                                  : (firstN[i].accepted_hi >> (j - 32)) & 1u;
        if (acc) acc_gpu[j] += nd.rate;
      }
      // CPU double-precision reference (transport_core)
      State s0; s0.r_m = {nd.in.rx, nd.in.ry, nd.in.rz};
      s0.n = unit({nd.in.tx, nd.in.ty, 1.0});
      auto cr = propagate(s0, nd.in.p_gev, fset, planes, tref);
      bool alive = true; unsigned long long cpu_mask = 0;
      for (size_t j = 0; j < planes.size() && alive; ++j) {
        const Crossing* x = nullptr;
        for (const auto& c : cr) if (c.plane_index == (int)j) { x = &c; break; }
        if (!x) break;
        if (in_aperture(planes[j], x->r)) { acc_cpu[j] += nd.rate; cpu_mask |= 1ull << j; }
        else if (planes[j].is_aperture) alive = false;
      }
      const unsigned long long gpu_mask =
          (unsigned long long)firstN[i].accepted_lo |
          ((unsigned long long)firstN[i].accepted_hi << 32);
      if (gpu_mask != cpu_mask) ++disagree;
    }
    const double cpu_dt = -[cpu_t0 timeIntervalSinceNow];
    std::printf("CPU checked %zu trajectories in %.2f s  (%.2e traj/s)\n",
                n, cpu_dt, n / std::max(cpu_dt, 1e-9));
    std::printf("CPU-check/GPU-full time ratio: %.2fx\n",
                cpu_dt / std::max(dt, 1e-9));
    std::printf("\nCPU-vs-GPU check, first %zu nodes: %lld nodes with any "
                "per-plane disagreement (edge/float effects expected at the "
                "sub-percent level)\n", n, disagree);
    std::printf("%-18s %14s %14s %10s\n", "plane", "cpu(/s)", "gpu(/s)", "gpu/cpu");
    for (size_t j = 0; j < planes.size(); ++j)
      if (acc_cpu[j] > 0 || acc_gpu[j] > 0)
        std::printf("%-18s %14.6e %14.6e %10.6f\n", planes[j].name.c_str(),
                    acc_cpu[j], acc_gpu[j],
                    acc_cpu[j] > 0 ? acc_gpu[j] / acc_cpu[j] : 0.0);
  }
  return 0;
} }
