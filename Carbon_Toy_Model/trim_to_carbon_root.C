#include <vector>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "Acts/Definitions/Units.hpp"

int trim_to_carbon_root(const char* infile,
                        const char* outfile,
                        const char* treename)
{
  TFile fin(infile, "READ");
  if (fin.IsZombie()) { fprintf(stderr, "ERROR: cannot open %s\n", infile); return 1; }
  TTree* tin = dynamic_cast<TTree*>(fin.Get(treename));
  if (!tin) { fprintf(stderr, "ERROR: tree '%s' not found\n", treename); return 2; }

  TFile fout(outfile, "RECREATE");

  // define carbon band, 
  struct CylBand { float R; float zmax; };

  //dR_mm = 10 //radial half width if you want to select tracks only near a single layer

  const std::vector<CylBand> bands = {
     { 1200,  3100.0f}   //Change according to which region you want, current set for pixel barrel region
//   { 987.3067993f, 1374.1672f }  // Example of selecting only: Strips::Barrel L8
  };

  auto in_any_band_xyz = [&](float x, float y, float z) -> bool {
    const float r = std::hypot(x, y);
    for (const auto& b : bands) {
      if (std::fabs(r) < b.R && std::fabs(z) < b.zmax) return true;
    }
    return false;
  };
  
  const float Z_C = 6.0f, A_C = 12.011f, rho_C = 2.0f * Acts::UnitConstants::g / Acts::UnitConstants::cm3;;
  const float X0_mass = 42.70f, L0_mass = 85.8f;
  const float X0_len_mm = 213.5f;
  const float L0_len_mm = 429.0f;
  float t_X0 = 0.0f;
  float t_L0 = 0.0f;
  
  std::vector<float>* v_mat_z = nullptr;
  std::vector<float>* v_mat_r = nullptr;
  std::vector<float>* v_mat_x = nullptr;
  std::vector<float>* v_mat_y = nullptr;

  std::vector<float>* v_mat_dx = nullptr;
  std::vector<float>* v_mat_dy = nullptr;
  std::vector<float>* v_mat_dz = nullptr;

  std::vector<float>* v_step = nullptr;
  std::vector<float>* v_Z = nullptr;
  std::vector<float>* v_A = nullptr;
  std::vector<float>* v_rho = nullptr;
  std::vector<float>* v_X0 = nullptr;
  std::vector<float>* v_L0 = nullptr;

  std::vector<std::vector<unsigned int>>* v_Zc = nullptr;
  std::vector<std::vector<float>>* v_Fc = nullptr;

  auto has = [&](const char* name){ return tin->GetBranch(name) != nullptr; };
  auto set = [&](const char* name, void* addr){ if (has(name)) tin->SetBranchAddress(name, addr); };

  set("t_X0", &t_X0);
  set("t_L0", &t_L0);


  set("mat_z", &v_mat_z);
  set("mat_r", &v_mat_r);
  set("mat_x", &v_mat_x);
  set("mat_y", &v_mat_y);

  set("mat_dx", &v_mat_dx);
  set("mat_dy", &v_mat_dy);
  set("mat_dz", &v_mat_dz);

  // required to keep ACTS happy, otherwise it gets messed up
  set("mat_step_length", &v_step);
  set("mat_Z", &v_Z);
  set("mat_A", &v_A);
  set("mat_rho", &v_rho);
  set("mat_X0", &v_X0);
  set("mat_L0", &v_L0);

  if (has("elements")) set("elements", &v_Zc);
  if (has("fraction")) set("fraction", &v_Fc);

  TTree* tout = tin->CloneTree(0);
  tout->SetBranchAddress("elements",  &v_Zc);
  tout->SetBranchAddress("fraction", &v_Fc);

  const Long64_t nentries = tin->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    tin->GetEntry(i);

    if (!v_mat_z || !v_mat_r) { fprintf(stderr, "ERROR: mat_z/r missing.\n"); return 3; }
    const size_t n0 = v_mat_z->size();
    if (v_mat_r->size() != n0 || v_step->size() != n0 || v_Z->size() != n0 ||
        v_A->size() != n0 || v_rho->size() != n0 || v_X0->size() != n0 || v_L0->size() != n0) {
      fprintf(stderr, "ERROR: size mismatch at entry %lld\n", i);
      return 4;
    }
    if (v_Zc && v_Zc->size() != n0) { fprintf(stderr, "ERROR: Zc outer size mismatch.\n"); return 5; }
    if (v_Fc && v_Fc->size() != n0) { fprintf(stderr, "ERROR: Fc outer size mismatch.\n"); return 6; }

  // collect indices to KEEP
  std::vector<size_t> keep;
  keep.reserve(n0);
  for (size_t s = 0; s < n0; ++s) {
    const float x = (*v_mat_x)[s];
    const float y = (*v_mat_y)[s];
    const float z = (*v_mat_z)[s];
    if (in_any_band_xyz(x, y, z)) keep.push_back(s);
  }
  if (keep.empty()) {
  // no steps in any band: skip this event entirely
    continue;
  }

    auto filter1 = [&](std::vector<float>& v){
      std::vector<float> out; out.reserve(keep.size());
      for (auto s : keep) out.push_back(v[s]);
      v.swap(out);
    };
    auto filter2 = [&](std::vector<std::vector<unsigned int>>& v){
      std::vector<std::vector<unsigned int>> out; out.reserve(keep.size());
      for (auto s : keep) out.push_back(v[s]);
      v.swap(out);
    };
    auto filter3 = [&](std::vector<std::vector<float>>& v){
      std::vector<std::vector<float>> out; out.reserve(keep.size());
      for (auto s : keep) out.push_back(v[s]);
      v.swap(out);
    };

    // apply the filter to all step-level
    filter1(*v_mat_z);
    filter1(*v_mat_x);
    filter1(*v_mat_y);
    if (v_mat_dx) filter1(*v_mat_dx);
    if (v_mat_dy) filter1(*v_mat_dy);
    if (v_mat_dz) filter1(*v_mat_dz);

    filter1(*v_mat_r);
    filter1(*v_step);
    filter1(*v_Z);
    filter1(*v_A);
    filter1(*v_rho);
    filter1(*v_X0);
    filter1(*v_L0);
    if (v_Zc) filter2(*v_Zc);
    if (v_Fc) filter3(*v_Fc);

    // set the mat information to be carbon's
    const size_t n = v_mat_z->size();
    const float scale = 0.05;
    for (size_t s = 0; s < n; ++s) {
      float &st = (*v_step)[s];
      if (!std::isfinite(st) || st <= 0.0f) st = 1e-6f;
      st *= scale;
      (*v_Z)[s] = Z_C;
      (*v_A)[s] = A_C;
      (*v_rho)[s] = rho_C;
      (*v_X0)[s] = X0_len_mm;
      (*v_L0)[s] = L0_len_mm;

      if (v_Zc) { (*v_Zc)[s].assign(1, static_cast<unsigned int>(Z_C)); }
      if (v_Fc) { (*v_Fc)[s].assign(1, 1.0f); }
    }

    t_X0 = 0.0f;
    t_L0 = 0.0f;
    for (size_t s = 0; s < v_step->size(); ++s) {
      const float ds = (*v_step)[s];
      const float X0 = (*v_X0)[s];
      const float L0 = (*v_L0)[s];
      if (std::isfinite(ds) && ds > 0.f && std::isfinite(X0) && X0 > 0.f) t_X0 += ds / X0;
      if (std::isfinite(ds) && ds > 0.f && std::isfinite(L0) && L0 > 0.f) t_L0 += ds / L0;
  }

    tout->Fill();
  }

  fout.Write();
  fout.Close();
  fin.Close();
  printf("Wrote trimmed file to %s\n", outfile);
  return 0;
}
