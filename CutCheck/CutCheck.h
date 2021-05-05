///////////////////////////////////////////////////////////////////////
// Class:       CutCheck
// Plugin Type: ******
// File:        CutCheck.h
////////////////////////////////////////////////////////////////////////
#include "protoduneana/protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "art_root_io/TFileService.h"
//#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FileBlock.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <fstream>
#include <string>
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

#include "StoppingMuonSelection/DataTypes.h"
#include "StoppingMuonSelection/GeometryHelper.h"
#include "StoppingMuonSelection/CutCheck/CutCheckHelper.h"

namespace stoppingcosmicmuonselection {

class CutCheck;

class CutCheck : public art::EDAnalyzer {
public:
  explicit CutCheck(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CutCheck(CutCheck const &) = delete;
  CutCheck(CutCheck &&) = delete;
  CutCheck & operator = (CutCheck const &) = delete;
  CutCheck & operator = (CutCheck &&) = delete;

  // Required functions.
  void analyze(art::Event const &evt) override;

  // Selected optional functions
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const& p);

private:
  size_t evNumber;

  CutCheckHelper           cutCheckHelper; // need configuration
  GeometryHelper           geoHelper;

  // Parameters form FHICL File
  double _trackPitch;
  double _trackPitchTolerance;
  bool _selectAC, _selectCC;
  bool _runCathodeSimple;
  std::string fPFParticleTag, fTrackerTag, fSpacePointTag;

  // Dummy histo to count events.
  TH1D *h_events;

  // Histos
  TH2D *h_dQdxVsRR;
  TH2D *h_dQdxVsRR_TP;

  TH1D *h_startX;
  TH1D *h_startX_signal;
  TH1D *h_startY;
  TH1D *h_startY_signal;
  TH1D *h_startZ;
  TH1D *h_startZ_signal;
  TH1D *h_endX;
  TH1D *h_endX_signal;
  TH1D *h_endY;
  TH1D *h_endY_signal;
  TH1D *h_endZ;
  TH1D *h_endZ_signal;
  TH1D *h_minHitPeakTime;
  TH1D *h_minHitPeakTime_signal;
  TH1D *h_maxHitPeakTime;
  TH1D *h_maxHitPeakTime_signal;
  TH1D *h_endX_Fabio;
  TH1D *h_endX_signal_Fabio;
  TH1D *h_endX_Pandora;
  TH1D *h_endX_signal_Pandora;
  TH1D *h_theta_xz;
  TH1D *h_theta_xz_signal;
  TH1D *h_theta_yz;
  TH1D *h_theta_yz_signal;
  TH1D *h_length;
  TH1D *h_length_signal;

  TH1D *h_startXPriori;
  TH1D *h_startX_signalPriori;
  TH1D *h_startYPriori;
  TH1D *h_startY_signalPriori;
  TH1D *h_startZPriori;
  TH1D *h_startZ_signalPriori;
  TH1D *h_endXPriori;
  TH1D *h_endX_signalPriori;
  TH1D *h_endYPriori;
  TH1D *h_endY_signalPriori;
  TH1D *h_endZPriori;
  TH1D *h_endZ_signalPriori;
  TH1D *h_minHitPeakTimePriori;
  TH1D *h_minHitPeakTime_signalPriori;
  TH1D *h_maxHitPeakTimePriori;
  TH1D *h_maxHitPeakTime_signalPriori;

};

CutCheck::CutCheck(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  reconfigure(p);
}

void CutCheck::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  h_events = tfs->make<TH1D>("h_events","h_events",1,0,1);

  art::TFileDirectory PrioriCutCheckDir = tfs->mkdir("PrioriCutCheck");
  h_startXPriori = PrioriCutCheckDir.make<TH1D>("h_startXPriori","h_startXPriori",400,-400,400);
  h_startX_signalPriori = PrioriCutCheckDir.make<TH1D>("h_startX_signalPriori","h_startX_signalPriori",400,-400,400);
  h_startYPriori = PrioriCutCheckDir.make<TH1D>("h_startYPriori","h_startYPriori",400,-100,700);
  h_startY_signalPriori = PrioriCutCheckDir.make<TH1D>("h_startY_signalPriori","h_startY_signalPriori",400,-100,700);
  h_startZPriori = PrioriCutCheckDir.make<TH1D>("h_startZPriori","h_startZPriori",450,-100,800);
  h_startZ_signalPriori = PrioriCutCheckDir.make<TH1D>("h_startZ_signalPriori","h_startZ_signalPriori",450,-100,800);
  h_endXPriori = PrioriCutCheckDir.make<TH1D>("h_endXPriori","h_endXPriori",400,-400,400);
  h_endX_signalPriori = PrioriCutCheckDir.make<TH1D>("h_endX_signalPriori","h_endX_signalPriori",400,-400,400);
  h_endYPriori = PrioriCutCheckDir.make<TH1D>("h_endYPriori","h_endYPriori",400,-100,700);
  h_endY_signalPriori = PrioriCutCheckDir.make<TH1D>("h_endY_signalPriori","h_endY_signalPriori",400,-100,700);
  h_endZPriori = PrioriCutCheckDir.make<TH1D>("h_endZPriori","h_endZPriori",450,-100,800);
  h_endZ_signalPriori = PrioriCutCheckDir.make<TH1D>("h_endZ_signalPriori","h_endZ_signalPriori",450,-100,800);
  h_minHitPeakTimePriori = PrioriCutCheckDir.make<TH1D>("h_minHitPeakTimePriori","h_minHitPeakTimePriori",3000,0,6000);
  h_minHitPeakTime_signalPriori = PrioriCutCheckDir.make<TH1D>("h_minHitPeakTime_signalPriori","h_minHitPeakTime_signalPriori",3000,0,6000);
  h_maxHitPeakTimePriori = PrioriCutCheckDir.make<TH1D>("h_maxHitPeakTimePriori","h_maxHitPeakTimePriori",3000,0,6000);
  h_maxHitPeakTime_signalPriori = PrioriCutCheckDir.make<TH1D>("h_maxHitPeakTime_signalPriori","h_maxHitPeakTime_signalPriori",3000,0,6000);

  art::TFileDirectory NMinus1Dir = tfs->mkdir("NMinus1");
  // Histograms
  h_dQdxVsRR = NMinus1Dir.make<TH2D>("h_dQdxVsRR","h_dQdxVsRR",200,0,200,800,0,800);
  h_dQdxVsRR_TP = NMinus1Dir.make<TH2D>("h_dQdxVsRR_TP","h_dQdxVsRR_TP",200,0,200,800,0,800);
  h_startX = NMinus1Dir.make<TH1D>("h_startX","h_startX",400,-400,400);
  h_startX_signal = NMinus1Dir.make<TH1D>("h_startX_signal","h_startX_signal",400,-400,400);
  h_startY = NMinus1Dir.make<TH1D>("h_startY","h_startY",400,-100,700);
  h_startY_signal = NMinus1Dir.make<TH1D>("h_startY_signal","h_startY_signal",400,-100,700);
  h_startZ = NMinus1Dir.make<TH1D>("h_startZ","h_startZ",450,-100,800);
  h_startZ_signal = NMinus1Dir.make<TH1D>("h_startZ_signal","h_startZ_signal",450,-100,800);
  h_endX = NMinus1Dir.make<TH1D>("h_endX","h_endX",400,-400,400);
  h_endX_signal = NMinus1Dir.make<TH1D>("h_endX_signal","h_endX_signal",400,-400,400);
  h_endY = NMinus1Dir.make<TH1D>("h_endY","h_endY",400,-100,700);
  h_endY_signal = NMinus1Dir.make<TH1D>("h_endY_signal","h_endY_signal",400,-100,700);
  h_endZ = NMinus1Dir.make<TH1D>("h_endZ","h_endZ",450,-100,800);
  h_endZ_signal = NMinus1Dir.make<TH1D>("h_endZ_signal","h_endZ_signal",450,-100,800);
  h_minHitPeakTime = NMinus1Dir.make<TH1D>("h_minHitPeakTime","h_minHitPeakTime",3000,0,6000);
  h_minHitPeakTime_signal = NMinus1Dir.make<TH1D>("h_minHitPeakTime_signal","h_minHitPeakTime_signal",3000,0,6000);
  h_maxHitPeakTime = NMinus1Dir.make<TH1D>("h_maxHitPeakTime","h_maxHitPeakTime",3000,0,6000);
  h_maxHitPeakTime_signal = NMinus1Dir.make<TH1D>("h_maxHitPeakTime_signal","h_maxHitPeakTime_signal",3000,0,6000);
  h_endX_Fabio = NMinus1Dir.make<TH1D>("h_endX_Fabio","h_endX_Fabio",400,-400,400);
  h_endX_signal_Fabio = NMinus1Dir.make<TH1D>("h_endX_signal_Fabio","h_endX_signal_Fabio",400,-400,400);
  h_endX_Pandora = NMinus1Dir.make<TH1D>("h_endX_Pandora","h_endX_Pandora",400,-400,400);
  h_endX_signal_Pandora = NMinus1Dir.make<TH1D>("h_endX_signal_Pandora","h_endX_signal_Pandora",400,-400,400);
  h_theta_xz = NMinus1Dir.make<TH1D>("h_theta_xz", "h_theta_xz", 180*2,-180,180);
  h_theta_xz_signal = NMinus1Dir.make<TH1D>("h_theta_xz_signal", "h_theta_xz_signal", 180*2,-180,180);
  h_theta_yz = NMinus1Dir.make<TH1D>("h_theta_yz", "h_theta_yz", 180,0,180);
  h_theta_yz_signal = NMinus1Dir.make<TH1D>("h_theta_yz_signal", "h_theta_yz_signal", 180,0,180);
  h_length = NMinus1Dir.make<TH1D>("h_length", "h_length", 800, 0, 800);
  h_length_signal = NMinus1Dir.make<TH1D>("h_length_signal", "h_length_signal", 800, 0, 800);
  
  // Print active volume bounds.
  geoHelper.PrintActiveVolumeBounds();
}

void CutCheck::endJob()
{
  mf::LogVerbatim("CutCheck") << "CutCheck finished job";
}

void CutCheck::reconfigure(fhicl::ParameterSet const& p)
{
  fTrackerTag = p.get<std::string>("TrackerTag");
  fPFParticleTag = p.get<std::string>("PFParticleTag");
  fSpacePointTag = p.get<std::string>("SpacePointTag");
  _selectAC = p.get<bool>("selectAC", true);
  _selectCC = p.get<bool>("selectCC", true);
  _runCathodeSimple = p.get<bool>("runCathodeSimple", false);
  _trackPitch = p.get<double>("trackPitch", 0.75);
  _trackPitchTolerance = p.get<double>("trackPitchTolerance", 0.1);
  cutCheckHelper.reconfigure(p.get<fhicl::ParameterSet>("ConfigSubModules"));
}

} // namespace
