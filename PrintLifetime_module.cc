///////////////////////////////////////////////////////////////////////
// Class:       PrintLifetime
// Plugin Type: ******
// File:        PrintLifetime.h
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
#include "TTimeStamp.h"

#include "protoduneana/StoppingMuonSelection/DataTypes.h"
#include "protoduneana/StoppingMuonSelection/GeometryHelper.h"
#include "protoduneana/StoppingMuonSelection/CalibrationHelper.h"

namespace stoppingcosmicmuonselection {

class PrintLifetime;

class PrintLifetime : public art::EDAnalyzer {
public:
  explicit PrintLifetime(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PrintLifetime(PrintLifetime const &) = delete;
  PrintLifetime(PrintLifetime &&) = delete;
  PrintLifetime & operator = (PrintLifetime const &) = delete;
  PrintLifetime & operator = (PrintLifetime &&) = delete;

  // Required functions.
  void analyze(art::Event const &evt) override;

  // Selected optional functions
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const& p);
  void respondToOpenInputFile(art::FileBlock const &inputFile) override;


private:
  std::string filename;

};

PrintLifetime::PrintLifetime(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  reconfigure(p);
}

void PrintLifetime::beginJob()
{
  std::cout << "BEGIN JOB" << std::endl;
  art::ServiceHandle<calib::LifetimeCalibService> lifetimecalibHandler;
  calib::LifetimeCalibService & lifetimecalibService = *lifetimecalibHandler;
  calib::LifetimeCalib *lifetimecalib = lifetimecalibService.provider();
  double fLifetime = lifetimecalib->GetLifetime()*1000.0; // [ms]*1000.0 -> [us]
  std::cout << "LIFETIME BEGIN JOB: " << fLifetime << std::endl;

}

void PrintLifetime::endJob()
{
  mf::LogVerbatim("PrintLifetime") << "PrintLifetime finished job";
}

void PrintLifetime::reconfigure(fhicl::ParameterSet const& p)
{}

void PrintLifetime::respondToOpenInputFile(art::FileBlock const &inputFile) {
  filename = inputFile.fileName();
  std::cout << "Analyzer on file: " << filename << std::endl;
}

void PrintLifetime::analyze(art::Event const &evt)
{
  // increase counter and store event number
  double fEvNumber = evt.id().event();
  std::cout << "PrintLifetime_module is on event: " << fEvNumber << std::endl;
  mf::LogVerbatim("PrintLifetime") << "PrintLifetime module on event " << fEvNumber;
  art::ServiceHandle<calib::LifetimeCalibService> lifetimecalibHandler;
  calib::LifetimeCalibService & lifetimecalibService = *lifetimecalibHandler;
  calib::LifetimeCalib *lifetimecalib = lifetimecalibService.provider();
  double fLifetime = lifetimecalib->GetLifetime()*1000.0; // [ms]*1000.0 -> [us]
  std::cout << "LIFETIME: " << fLifetime << std::endl;

  // Timing stuff
  const char * timestamp;
  art::Timestamp ts = evt.time();
  if (ts.timeHigh()==0) {
    TTimeStamp ts2(ts.timeLow());
    timestamp = ts2.AsString();
  }
  else {
    TTimeStamp ts2(ts.timeHigh(), ts.timeLow());
    timestamp = ts2.AsString();
  }
  std::cout << "TIMESTAMP: "  << timestamp << std::endl;

}

} // namespace
DEFINE_ART_MODULE(stoppingcosmicmuonselection::PrintLifetime)
