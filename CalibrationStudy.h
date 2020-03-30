///////////////////////////////////////////////////////////////////////
// Class:       CalibrationStudy
// Plugin Type: ******
// File:        CalibrationStudy.h
////////////////////////////////////////////////////////////////////////
#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
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
#include "TProfile2D.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"

#include "DataTypes.h"
#include "GeometryHelper.h"
#include "SpacePointAlg.h"
#include "CalorimetryHelper.h"
#include "StoppingMuonSelectionAlg.h"
#include "HitHelper.h"
#include "HitPlaneAlg.h"
#include "CNNHelper.h"

namespace stoppingcosmicmuonselection {

class CalibrationStudy;

class CalibrationStudy : public art::EDAnalyzer {
public:
  explicit CalibrationStudy(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CalibrationStudy(CalibrationStudy const &) = delete;
  CalibrationStudy(CalibrationStudy &&) = delete;
  CalibrationStudy & operator = (CalibrationStudy const &) = delete;
  CalibrationStudy & operator = (CalibrationStudy &&) = delete;

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
  int counter_T0_tagged_tracks = 0;
  int counter_total_number_events = 0;
  int counter_total_number_tracks = 0;

  GeometryHelper           geoHelper;
  SpacePointAlg            spAlg;        // need configuration
  StoppingMuonSelectionAlg selectorAlg;  // need configuration
  CalorimetryHelper        caloHelper;   // need configuration
  HitHelper                hitHelper;    // need configuration
  CNNHelper             cnnHelper;

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
  // Objects for TTree
  std::string filename;
  // TH1s
  TH1D *fh_progressiveDistance = nullptr;
  // TProfiles
  TGraph2D *fg_imageCollection = nullptr;
  TGraph2D *fg_imageScore = nullptr;
  TGraph2D *fg_imageCollectionNoMichel = nullptr;

  // Graphs
  TGraphErrors *fg_wireID = nullptr;
  TGraphErrors *fg_Q = nullptr;
  TGraphErrors *fg_Dqds = nullptr;
  TGraphErrors *fg_QSmooth = nullptr;
  TGraphErrors *fg_DqdsSmooth = nullptr;
  TGraphErrors *fg_LocalLin = nullptr;
  TGraphErrors *fg_CnnScore = nullptr;

  // Histos
  TH2D *h_dQdxVsRR;
  TH2D *h_dQdxVsRR_TP075;
  TH2D *h_dQdxVsRR_LTCorr;
  TH2D *h_dQdxVsRR_TP075_LTCorr;
  TH2D *h_dQdEVsRR_TP075_LTCorr_MC;
  TH2D *h_dQdEVsRR_TP075_LTCorr_LV;
  TH2D *h_dQdxVsRR_NoMichel;
  TH2D *h_dQdxVsRR_NoMichelTP;
};

CalibrationStudy::CalibrationStudy(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  reconfigure(p);
}

void CalibrationStudy::beginJob()
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
  fTrackTree->Branch("distEndPoint", &fDistEndPoint, "fDistEndPoint/d");
  fTrackTree->Branch("distEndPointNoMichel", &fDistEndPointNoMichel, "fDistEndPointNoMichel/d");
  fTrackTree->Branch("isRecoSelectedCathodeCrosser",&fIsRecoSelectedCathodeCrosser);
  fTrackTree->Branch("isTrueSelectedCathodeCrosser",&fIsTrueSelectedCathodeCrosser);
  fTrackTree->Branch("isRecoSelectedAnodeCrosser",&fIsRecoSelectedAnodeCrosser);
  fTrackTree->Branch("isTrueSelectedAnodeCrosser",&fIsTrueSelectedAnodeCrosser);
  fTrackTree->Branch("isAnodePandora", &fIsAnodePandora);
  fTrackTree->Branch("isAnodeMine", &fIsAnodeMine);
  fTrackTree->Branch("filename", &filename);
  fTrackTree->Branch("g_imageCollection",&fg_imageCollection);
  fTrackTree->Branch("g_imageCollectionNoMichel",&fg_imageCollectionNoMichel);
  fTrackTree->Branch("g_imageScore",&fg_imageScore);
  fTrackTree->Branch("h_progressiveDistance","TH1D",&fh_progressiveDistance);
  fTrackTree->Branch("g_wireID", &fg_wireID);
  fTrackTree->Branch("g_Q", &fg_Q);
  fTrackTree->Branch("g_Dqds", &fg_Dqds);
  fTrackTree->Branch("g_QSmooth", &fg_QSmooth);
  fTrackTree->Branch("g_DqdsSmooth", &fg_DqdsSmooth);
  fTrackTree->Branch("g_LocalLin", &fg_LocalLin);
  fTrackTree->Branch("g_CnnScore", &fg_CnnScore);
  fTrackTree->Branch("michelHitsMichelScore", &f_michelHitsMichelScore);
  fTrackTree->Branch("muonHitsMichelScore", &f_muonHitsMichelScore);

  // Init the graph for the hits
  fg_imageCollection = new TGraph2D();
  fg_imageScore = new TGraph2D();
  fg_imageCollectionNoMichel = new TGraph2D();

  fh_progressiveDistance = new TH1D("h_progressiveDistance","h_progressiveDistance",200,0,200);

  // Histograms
  h_dQdxVsRR = tfs->make<TH2D>("h_dQdxVsRR","h_dQdxVsRR",200,0,200,800,0,800);
  h_dQdxVsRR_TP075 = tfs->make<TH2D>("h_dQdxVsRR_TP075","h_dQdxVsRR_TP075",200,0,200,800,0,800);
  h_dQdxVsRR_LTCorr = tfs->make<TH2D>("h_dQdxVsRR_LTCorr","h_dQdxVsRR_LTCorr",200,0,200,800,0,800);
  h_dQdxVsRR_TP075_LTCorr = tfs->make<TH2D>("h_dQdxVsRR_TP075_LTCorr","h_dQdxVsRR_TP075_LTCorr",200,0,200,800,0,800);
  h_dQdEVsRR_TP075_LTCorr_MC = tfs->make<TH2D>("h_dQdEVsRR_TP075_LTCorr_MC","h_dQdEVsRR_TP075_LTCorr_MC",200,0,200,800,0,800);
  h_dQdEVsRR_TP075_LTCorr_LV = tfs->make<TH2D>("h_dQdEVsRR_TP075_LTCorr_LV","h_dQdEVsRR_TP075_LTCorr_LV",200,0,200,800,0,800);
  h_dQdxVsRR_NoMichel = tfs->make<TH2D>("h_dQdxVsRR_NoMichel","h_dQdxVsRR_NoMichel",200,0,200,800,0,800);
  h_dQdxVsRR_NoMichelTP = tfs->make<TH2D>("h_dQdxVsRR_NoMichelTP","h_dQdxVsRR_NoMichelTP",200,0,200,800,0,800);

  // Graphs
  fg_wireID = new TGraphErrors();
  fg_Q = new TGraphErrors();
  fg_Dqds = new TGraphErrors();
  fg_QSmooth = new TGraphErrors();
  fg_DqdsSmooth = new TGraphErrors();
  fg_LocalLin = new TGraphErrors();
  fg_CnnScore = new TGraphErrors();

  // Print active volume bounds.
  geoHelper.PrintActiveVolumeBounds();

}

void CalibrationStudy::endJob()
{
  mf::LogVerbatim("CalibrationStudy") << "CalibrationStudy finished job";
  std::cout << "Total number of events: " << counter_total_number_events << std::endl;
  std::cout << "Total number of tracks: " << counter_total_number_tracks << std::endl;
  std::cout << "Number of T0-tagged tracks: " << counter_T0_tagged_tracks << std::endl;
}

void CalibrationStudy::respondToOpenInputFile(art::FileBlock const &inputFile) {
  filename = inputFile.fileName();
  std::cout << "Analyzer on file: " << filename << std::endl;
}

void CalibrationStudy::reconfigure(fhicl::ParameterSet const& p)
{
  fTrackerTag = p.get<std::string>("TrackerTag");
  fPFParticleTag = p.get<std::string>("PFParticleTag");
  fSpacePointTag = p.get<std::string>("SpacePointTag");
  fNNetTag = p.get<std::string>("NNetTag");
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

void CalibrationStudy::UpdateTTreeVariableWithTrackProperties(const trackProperties &trackInfo) {
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
