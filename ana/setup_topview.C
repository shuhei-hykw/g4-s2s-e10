// setup_topview.C
//
// Schematic top view (plan view: Z = beam direction horizontal,
// X = horizontal bend-plane coordinate vertical, Y/height projected
// out) of the S-2S spectrometer setup used for the 9Be(K-,pi+)
// 9_Sigma_He acceptance study.
//
// Coordinates are taken directly from
// param/DCGEO/DCGeomParam_e10_9hesigma_kpi (S-2S-coordinate entries
// only; the header notes BC3/BC4/BH1/BFT/BH2/BAC are given from a
// different reference point (VI/VO), so they are NOT included here
// to avoid mixing coordinate frames).  Q1/Q2 have no explicit
// envelope in that file (only the S2SQ1Q2Boundary/S2SQ2D1Boundary
// markers), so their boxes are an approximate schematic, not exact
// engineering dimensions -- labelled accordingly.
//
// This is not a Geant4 3D render (the native Geant4 vis pipeline is
// unavailable in this headless KEKCC session: no working X11 display,
// and the RayTracer driver in this MT-flavoured Geant4 11.2.2 build
// segfaults when used with this app's Serial run manager). It is a
// data-driven schematic drawn directly from the same DCGEO numbers
// the simulation itself uses, so the layout/scale is accurate.
//
// Usage: root -l -b -q ana/setup_topview.C

void setup_topview()
{
  gStyle->SetOptStat(0);

  // Central-trajectory points, in S-2S coordinates (mm), downstream
  // order: {z, x, label}.
  const int NP = 20;
  double pz[NP] = {
    -4900,    -4750,   -4550,    -3885.87, -3246.20, -2776.52, -2303.57,
    -2100.61, -1198.51,-379.893,  279.437,  718.455,  826.841,
     887.2,    1138.3,  1228.3,   1316.20,  1334.09,  1374.59,  1605.93
  };
  double px[NP] = {
    0,        0,       0,        0,        0,        0,        0,
    0,        138.849, 542.544,  1173.72,  1973.94,  2271.73,
    2443.5,   3130.4,  3378.8,   3616.24,  3665.39,  3776.65,  4412.25
  };
  const char* label[NP] = {
    "Target", "SFT", "SDC2", "VP1", "VP2/Q1-Q2", "VP3", "VP4/Q2-D1",
    "VP5/D1-in", "VP6", "VP7", "VP8", "VP9/D1-out", "VP10",
    "SDC3", "SDC4", "SDC5", "TOF", "RKINIT", "AC1", "WC"
  };
  // NOTE: pz/px above are (z,x) in mm; drawn with x-axis=z, y-axis=px
  // (the horizontal bend coordinate) to give a plan/top view.

  auto c = new TCanvas("c", "c", 1400, 500);
  c->SetLeftMargin(0.08); c->SetRightMargin(0.03);
  c->SetTopMargin(0.10);  c->SetBottomMargin(0.14);

  auto frame = new TH2D("frame", "S-2S setup: top view (schematic, from DCGEO);"
                         "Z [mm] (beam direction);X [mm] (bend plane)",
                         10, -5300, 5300, 10, -600, 5000);
  frame->Draw();

  auto traj = new TGraph(NP, pz, px);
  traj->SetLineColor(kBlue+2);
  traj->SetLineWidth(2);
  traj->Draw("L SAME");

  auto pts = new TGraph(NP, pz, px);
  pts->SetMarkerStyle(20);
  pts->SetMarkerSize(0.9);
  pts->SetMarkerColor(kBlue+2);
  pts->Draw("P SAME");

  // Only label the physically meaningful detectors/markers; the
  // VP-only points still show as markers/trajectory but without text,
  // to avoid overlapping labels in the tightly-packed regions.
  bool showLabel[NP] = {
    true, true, true, false, false, false, false,
    false, false, false, false, false, false,
    true, true, true, true, true, true, true
  };
  for(int i = 0; i < NP; ++i){
    if(!showLabel[i]) continue;
    // stagger alternating labels vertically to reduce overlap in the
    // tightly-packed downstream detector cluster.
    double dy = (i % 2 == 0) ? 220 : 420;
    auto* t = new TLatex(pz[i], px[i] + dy, label[i]);
    t->SetTextSize(0.026);
    t->SetTextAlign(21);
    t->Draw();
  }

  // Q1/Q2 schematic envelopes (approximate: no explicit length in
  // DCGEO, only the boundary markers; half-width is a generic
  // placeholder, not a real aperture).
  auto drawBox = [&](double z0, double z1, double halfw,
                      const char* name, Color_t col){
    auto box = new TBox(z0, -halfw, z1, halfw);
    box->SetFillColorAlpha(col, 0.15);
    box->SetLineColor(col);
    box->SetLineWidth(2);
    box->Draw();
    auto t = new TLatex(0.5*(z0+z1), halfw + 120, name);
    t->SetTextSize(0.03);
    t->SetTextAlign(21);
    t->SetTextColor(col);
    t->Draw();
  };
  drawBox(-3885.87, -3246.20, 150, "Q1", kRed+1);
  drawBox(-3246.20, -2303.57, 150, "Q2", kRed+1);

  // D1 dipole: simple polygon through its pole-boundary markers.
  const int nD1 = 5;
  double d1z[nD1] = {-2100.61, -1198.51, 279.437, 718.455, 826.841};
  double d1x[nD1] = {0,         138.849, 1173.72, 1973.94, 2271.73};
  auto d1 = new TGraph(nD1, d1z, d1x);
  d1->SetLineColor(kGreen+2);
  d1->SetLineWidth(4);
  d1->Draw("L SAME");
  auto d1lab = new TLatex(-700, 900, "D1 dipole (~70 deg bend)");
  d1lab->SetTextSize(0.03);
  d1lab->SetTextColor(kGreen+2);
  d1lab->SetTextAlign(21);
  d1lab->Draw();

  auto note = new TLatex(0.01, 0.02,
    "Schematic from DCGeomParam_e10_9hesigma_kpi (S-2S-coordinate "
    "entries only); not a full Geant4 3D render. Q1/Q2 boxes are an "
    "approximate envelope, not the real aperture.");
  note->SetNDC();
  note->SetTextSize(0.023);
  note->Draw();

  c->Print("ana/fig/setup_topview_schematic.pdf");
  c->Print("ana/fig/setup_topview_schematic.png");
  std::cout << "Saved: ana/fig/setup_topview_schematic.{pdf,png}"
            << std::endl;
}
