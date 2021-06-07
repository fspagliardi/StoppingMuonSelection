///////////////////////////////////////////////////////////////////////
// Class:       ModBoxModStudyAnode
// Plugin Type: ******
// File:        ModBoxModStudyAnode.h
////////////////////////////////////////////////////////////////////////
#include "protoduneana/protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "art_root_io/TFileService.h"
//#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FileBlock.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include <fstream>
#include <string>
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"

#include "protoduneana/StoppingMuonSelection/DataTypes.h"
#include "protoduneana/StoppingMuonSelection/GeometryHelper.h"
#include "protoduneana/StoppingMuonSelection/SpacePointAlg.h"
#include "protoduneana/StoppingMuonSelection/CalorimetryHelper.h"
#include "protoduneana/StoppingMuonSelection/StoppingMuonSelectionAlg.h"
#include "protoduneana/StoppingMuonSelection/HitHelper.h"
#include "protoduneana/StoppingMuonSelection/HitPlaneAlg.h"
#include "protoduneana/StoppingMuonSelection/CNNHelper.h"
#include "protoduneana/StoppingMuonSelection/SceHelper.h"
#include "protoduneana/StoppingMuonSelection/CalibrationHelper.h"
#include "protoduneana/StoppingMuonSelection/FixCalo.h"

namespace stoppingcosmicmuonselection {

class ModBoxModStudyAnode;

class ModBoxModStudyAnode : public art::EDAnalyzer {
public:
  explicit ModBoxModStudyAnode(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ModBoxModStudyAnode(ModBoxModStudyAnode const &) = delete;
  ModBoxModStudyAnode(ModBoxModStudyAnode &&) = delete;
  ModBoxModStudyAnode & operator = (ModBoxModStudyAnode const &) = delete;
  ModBoxModStudyAnode & operator = (ModBoxModStudyAnode &&) = delete;

  // Required functions.
  void analyze(art::Event const &evt) override;

  // Selected optional functions
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const& p);
  void respondToOpenInputFile(art::FileBlock const &inputFile) override;

  // Function called by the analyze function
  void UpdateTTreeVariableWithTrackProperties(const trackProperties &trackProp);

private:
  // Declare some counters for statistic purposes
  int counter_total_number_tracks = 0;

  GeometryHelper           geoHelper;
  SpacePointAlg            spAlg;        // need configuration
  StoppingMuonSelectionAlg selectorAlg;  // need configuration
  CalorimetryHelper        caloHelper;   // need configuration
  HitHelper                hitHelper;    // need configuration
  CNNHelper             cnnHelper;
  CalibrationHelper        calibHelper;
  SceHelper                *sceHelper;
  FixCalo                 fixCalo;

  // Parameters form FHICL File
  size_t _minNumbMichelLikeHit;
  double _trackPitch;
  double _trackPitchTolerance;
  size_t _numberNeighbors;
  double _michelScoreThreshold;
  double _michelScoreThresholdAvg;
  bool _selectAC, _selectCC;
  std::string fPFParticleTag, fSpacePointTag, fTrackerTag;
  std::string fNNetTag;

  // Utils
  protoana::ProtoDUNEPFParticleUtils   pfpUtil;
  protoana::ProtoDUNETruthUtils        truthUtil;

  // Track Tree stuff
  TTree *fTrackTree;
  // Tree variables
  size_t fEvNumber;
  int    fPdgID = INV_INT;
  double fTrackLength = INV_DBL;
  double fEndX = INV_DBL;
  double fEndY = INV_DBL;
  double fEndZ = INV_DBL;
  double fStartX = INV_DBL;
  double fStartY = INV_DBL;
  double fStartZ = INV_DBL;
  double fRecoTrackID = INV_DBL;
  double fTEndX = INV_DBL;
  double fTEndY = INV_DBL;
  double fTEndZ = INV_DBL;
  double fTStartX = INV_DBL;
  double fTStartY = INV_DBL;
  double fTStartZ = INV_DBL;
  double fTStartT = INV_DBL;
  double fTEndT = INV_DBL;
  double fT0_reco = INV_DBL;
  double fTrackID = INV_DBL;
  double ftheta_xz = INV_DBL;
  double ftheta_yz = INV_DBL;
  double fMinHitPeakTime = INV_DBL;
  double fMaxHitPeakTime = INV_DBL;
  double fEndX_corr = INV_DBL;
  double fEndY_corr = INV_DBL;
  double fEndZ_corr = INV_DBL;
  double fStartX_corr = INV_DBL;
  double fStartY_corr = INV_DBL;
  double fStartZ_corr = INV_DBL;
  double fDistEndPoint = INV_DBL;
  double fDistEndPointNoMichel = INV_DBL;
  bool fIsRecoSelectedCathodeCrosser = false;
  bool fIsRecoSelectedAnodeCrosser = false;
  bool fIsTrueSelectedCathodeCrosser = false;
  bool fIsTrueSelectedAnodeCrosser = false;
  bool fIsAnodePandora = false;
  bool fIsAnodeMine = false;
  std::vector<double> f_michelHitsMichelScore;
  std::vector<double> f_muonHitsMichelScore;
  std::vector<double> fDriftTime;
  std::vector<double> fLifeTimeCorr;
  std::vector<double> fYZcalibFactor;
  std::vector<double> fXcalibFactor;
  std::vector<double> fdQdx;
  std::vector<double> fdEdx;
  std::vector<double> fResRange;
  std::vector<double> fTrackPitch;
  std::vector<double> fHitX;
  std::vector<double> fHitY;
  std::vector<double> fHitZ;
  std::vector<double> fPhis;
  std::vector<double> fHitAmpl;
  std::vector<double> fHitRMS;
  std::vector<double> fEfX;
  std::vector<double> fEfY;
  std::vector<double> fEfZ;
  std::vector<double> fEfield;

  // Objects for TTree
  std::string filename;

  // Histos
  TH2D *h_dQdxVsRR;
  TH2D *h_dQdxVsRR_TP075;

  TH2D *h_hitYZ;
  TH2D *h_hitXZ;
  TH2D *h_hitXY;

};

ModBoxModStudyAnode::ModBoxModStudyAnode(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  reconfigure(p);
}

void ModBoxModStudyAnode::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTrackTree = tfs->make<TTree>("TrackTree", "track by track info");
  fTrackTree->Branch("event", &fEvNumber, "fEvNumber/l");
  fTrackTree->Branch("PdgID", &fPdgID);
  fTrackTree->Branch("trackLength", &fTrackLength, "fTrackLength/d");
  fTrackTree->Branch("endX", &fEndX, "fEndX/d");
  fTrackTree->Branch("endY", &fEndY, "fEndY/d");
  fTrackTree->Branch("endZ", &fEndZ, "fEndZ/d");
  fTrackTree->Branch("startX", &fStartX, "fStartX/d");
  fTrackTree->Branch("startY", &fStartY, "fStartY/d");
  fTrackTree->Branch("startZ", &fStartZ, "fStartZ/d");
  fTrackTree->Branch("recoTrackID", &fRecoTrackID);
  fTrackTree->Branch("TEndX", &fTEndX, "fTEndX/d");
  fTrackTree->Branch("TEndY", &fTEndY, "fTEndY/d");
  fTrackTree->Branch("TEndZ", &fTEndZ, "fTEndZ/d");
  fTrackTree->Branch("TStartX", &fTStartX, "fStartX/d");
  fTrackTree->Branch("TStartY", &fTStartY, "fTStartY/d");
  fTrackTree->Branch("TStartZ", &fTStartZ, "fTStartZ/d");
  fTrackTree->Branch("TStartT", &fTStartT, "fTStartT/d");
  fTrackTree->Branch("TEndT", &fTEndT, "fTEndT/d");
  fTrackTree->Branch("trackID", &fTrackID, "fTrackID/d");
  fTrackTree->Branch("T0_reco", &fT0_reco, "fT0_reco/d");
  fTrackTree->Branch("minHitPeakTime", &fMinHitPeakTime, "fMinHitPeakTime/d");
  fTrackTree->Branch("maxHitPeakTime", &fMaxHitPeakTime, "fMaxHitPeakTime/d");
  fTrackTree->Branch("theta_xz", &ftheta_xz, "ftheta_xz/d");
  fTrackTree->Branch("theta_yz", &ftheta_yz, "ftheta_yz/d");
  fTrackTree->Branch("endX_corr", &fEndX_corr, "fEndX_corr/d");
  fTrackTree->Branch("endY_corr", &fEndY_corr, "fEndY_corr/d");
  fTrackTree->Branch("endZ_corr", &fEndZ_corr, "fEndZ_corr/d");
  fTrackTree->Branch("startX_corr", &fStartX_corr, "fStartX_corr/d");
  fTrackTree->Branch("startY_corr", &fStartY_corr, "fStartY_corr/d");
  fTrackTree->Branch("startZ_corr", &fStartZ_corr, "fStartZ_corr/d");
  fTrackTree->Branch("filename", &filename);
  fTrackTree->Branch("distEndPoint", &fDistEndPoint, "fDistEndPoint/d");
  fTrackTree->Branch("distEndPointNoMichel", &fDistEndPointNoMichel, "fDistEndPointNoMichel/d");
  fTrackTree->Branch("isRecoSelectedCathodeCrosser",&fIsRecoSelectedCathodeCrosser);
  fTrackTree->Branch("isTrueSelectedCathodeCrosser",&fIsTrueSelectedCathodeCrosser);
  fTrackTree->Branch("isRecoSelectedAnodeCrosser",&fIsRecoSelectedAnodeCrosser);
  fTrackTree->Branch("isTrueSelectedAnodeCrosser",&fIsTrueSelectedAnodeCrosser);
  fTrackTree->Branch("isAnodePandora", &fIsAnodePandora);
  fTrackTree->Branch("isAnodeMine", &fIsAnodeMine);
  fTrackTree->Branch("driftTime", &fDriftTime);
  fTrackTree->Branch("lifeTimeCorr", &fLifeTimeCorr);
  fTrackTree->Branch("YZcalibFactor", &fYZcalibFactor);
  fTrackTree->Branch("XcalibFactor", &fXcalibFactor);
  fTrackTree->Branch("dQdx", &fdQdx);
  fTrackTree->Branch("dEdx", &fdEdx);
  fTrackTree->Branch("ResRange", &fResRange);
  fTrackTree->Branch("TrackPitch", &fTrackPitch);
  fTrackTree->Branch("HitX", &fHitX);
  fTrackTree->Branch("HitY", &fHitY);
  fTrackTree->Branch("HitZ", &fHitZ);
  fTrackTree->Branch("Phis", &fPhis);
  fTrackTree->Branch("HitAmpl", &fHitAmpl);
  fTrackTree->Branch("HitRMS", &fHitRMS);
  fTrackTree->Branch("EfX", &fEfX);
  fTrackTree->Branch("EfY", &fEfY);
  fTrackTree->Branch("EfZ", &fEfZ);
  fTrackTree->Branch("Efield", &fEfield);

  // Histograms
  h_dQdxVsRR = tfs->make<TH2D>("h_dQdxVsRR","h_dQdxVsRR",200,0,200,800,0,800);
  h_dQdxVsRR_TP075 = tfs->make<TH2D>("h_dQdxVsRR_TP075","h_dQdxVsRR_TP075",200,0,200,800,0,800);

  h_hitYZ = tfs->make<TH2D>("h_hitYZ","h_hitYZ",800,-50,750,700,-50,650);
  h_hitXZ = tfs->make<TH2D>("h_hitXZ","h_hitXZ",800,-50,750,760,-380,380);
  h_hitXY = tfs->make<TH2D>("h_hitXY","h_hitXY",760,-380,380,700,-50,650);

  // h_dQdx_Cosphi = tfs->make<TH2D>("h_dQdx_Cosphi", "h_dQdx_Cosphi", 48, 0, 1, 800, 0, 800);
  // h_dQdx_Cosphi_corr = tfs->make<TH2D>("h_dQdx_Cosphi_corr", "h_dQdx_Cosphi_corr", 48, 0, 1, 800, 0, 800);
  //
  // h_dQdx_phi = tfs->make<TH2D>("h_dQdx_phi", "h_dQdx_phi", 180, 0, 180, 800, 0, 800);
  // h_dQdx_phi_corr = tfs->make<TH2D>("h_dQdx_phi_corr", "h_dQdx_phi_corr", 180, 0, 180, 800, 0, 800);
  //
  // h_dQdx_thetaxy = tfs->make<TH2D>("h_dQdx_thetaxy", "h_dQdx_thetaxy", 200, -200, 200, 800, 0, 800);
  // h_dQdx_thetaxy_corr = tfs->make<TH2D>("h_dQdx_thetaxy_corr", "h_dQdx_thetaxy_corr", 200, -200, 200, 800, 0, 800);

  // Print active volume bounds.
  geoHelper.PrintActiveVolumeBounds();

}

void ModBoxModStudyAnode::endJob()
{
  mf::LogVerbatim("ModBoxModStudyAnode") << "ModBoxModStudyAnode finished job";
  std::cout << "Total number of tracks: " << counter_total_number_tracks << std::endl;
}

void ModBoxModStudyAnode::respondToOpenInputFile(art::FileBlock const &inputFile) {
  filename = inputFile.fileName();
  std::cout << "Analyzer on file: " << filename << std::endl;
}

void ModBoxModStudyAnode::reconfigure(fhicl::ParameterSet const& p)
{
  fTrackerTag = p.get<std::string>("TrackerTag");
  fPFParticleTag = p.get<std::string>("PFParticleTag");
  fNNetTag = p.get<std::string>("NNetTag");
  fSpacePointTag = p.get<std::string>("SpacePointTag");
  _minNumbMichelLikeHit = p.get<size_t>("minNumbMichelLikeHit", 2);
  _trackPitch = p.get<double>("trackPitch", 0.75);
  _trackPitchTolerance = p.get<double>("trackPitchTolerance", 0.1);
  _numberNeighbors = p.get<size_t>("numberNeighbors", 2);
  _michelScoreThreshold = p.get<double>("michelScoreThreshold", 0.7);
  _michelScoreThresholdAvg = p.get<double>("michelScoreThresholdAvg", 0.5);
  _selectAC = p.get<bool>("selectAC", true);
  _selectCC = p.get<bool>("selectCC", true);
  spAlg.reconfigure(p.get<fhicl::ParameterSet>("SpacePointAlg"));
  selectorAlg.reconfigure(p.get<fhicl::ParameterSet>("StoppingMuonSelectionAlg"));
  caloHelper.reconfigure(p.get<fhicl::ParameterSet>("CalorimetryHelper"));
  hitHelper.reconfigure(p.get<fhicl::ParameterSet>("HitHelper"));
}

void ModBoxModStudyAnode::UpdateTTreeVariableWithTrackProperties(const trackProperties &trackInfo) {
  fEvNumber       = trackInfo.evNumber;
  fT0_reco        = trackInfo.trackT0;
  fStartX         = trackInfo.recoStartPoint.X();
  fStartY         = trackInfo.recoStartPoint.Y();
  fStartZ         = trackInfo.recoStartPoint.Z();
  fEndX           = trackInfo.recoEndPoint.X();
  fEndY           = trackInfo.recoEndPoint.Y();
  fEndZ           = trackInfo.recoEndPoint.Z();
  ftheta_xz       = trackInfo.theta_xz;
  ftheta_yz       = trackInfo.theta_yz;
  fMinHitPeakTime = trackInfo.minHitPeakTime;
  fMaxHitPeakTime = trackInfo.maxHitPeakTime;
  fTrackLength    = trackInfo.trackLength;
  fRecoTrackID    = trackInfo.trackID;
  fPdgID          = trackInfo.pdg;
  fTStartX        = trackInfo.trueStartPoint.X();
  fTStartY        = trackInfo.trueStartPoint.Y();
  fTStartZ        = trackInfo.trueStartPoint.Z();
  fTEndX          = trackInfo.trueEndPoint.X();
  fTEndY          = trackInfo.trueEndPoint.Y();
  fTEndZ          = trackInfo.trueEndPoint.Z();
  fTStartT        = trackInfo.trueStartT;;
  fTEndT          = trackInfo.trueEndT;
  fTrackID        = trackInfo.trueTrackID;
  fIsAnodePandora = trackInfo.isAnodeCrosserPandora;
  fIsAnodeMine    = trackInfo.isAnodeCrosserMine;
}

} // namespace
