#include "RecoMuon/StandAloneTrackFinder/interface/ExhaustiveMuonTrajectoryBuilder.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHitBuilder.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"

ExhaustiveMuonTrajectoryBuilder::ExhaustiveMuonTrajectoryBuilder(const edm::ParameterSet& pset,
                                                                 const MuonServiceProxy* proxy,
                                                                 edm::ConsumesCollector& iC)
    : theTrajBuilder(pset, proxy, iC), theSeeder(), theService(proxy) {}

ExhaustiveMuonTrajectoryBuilder::~ExhaustiveMuonTrajectoryBuilder() {}

MuonTrajectoryBuilder::TrajectoryContainer ExhaustiveMuonTrajectoryBuilder::trajectories(const TrajectorySeed& seed) {
  std::cout << "In ExhaustiveMuonTrajectoryBuilder::trajectories" << std::endl;
  LocalTrajectoryParameters localTrajectoryParameters(seed.startingState().parameters());
  LocalVector p(localTrajectoryParameters.momentum());
  int rawId = seed.startingState().detId();
  DetId detId(rawId);
  bool isBarrel = (detId.subdetId() == 1);
  // homemade local-to-global
  double pt = (isBarrel) ? -p.z() : p.perp();
  pt *= localTrajectoryParameters.charge();
  float err00 = seed.startingState().error(0);
  //   float p_err = sqr(sptmean/(ptmean*ptmean));
  //  mat[0][0]= p_err;
  float sigmapt = sqrt(err00) * pt * pt;
  TrajectoryContainer result;
  std::cout << "Preliminary parameters set, looping over recHits" << std::endl;
  // Make a new seed based on each segment, using the original pt and sigmapt
  for (auto const& recHit : seed.recHits()) {
    auto geoId = recHit.geographicalId();

    if (geoId.subdetId() != MuonSubdetId::GEM) {
      std::cout << "RecHit in seed not in GEM: " << geoId.rawId() << std::endl;
    } else if (geoId.subdetId() == MuonSubdetId::GEM) {
      GEMDetId chamberId(geoId.rawId());
      std::cout << "RecHit in seed in GEM: \n"
                << "Global Id: " << geoId.rawId() << ", chamber Id: " << chamberId << " Is ME0? "
                << (chamberId.isME0() ? "Yes" : "No") << std::endl;
    }

    const GeomDet* geomDet = theService->trackingGeometry()->idToDet(recHit.geographicalId());
    std::cout << geomDet << std::endl;
    auto muonRecHit = MuonTransientTrackingRecHit::specificBuild(geomDet, &recHit);
    TrajectorySeed tmpSeed(theSeeder.createSeed(pt, sigmapt, muonRecHit));
    std::cout << "Number of hits in the created seed: " << tmpSeed.nHits() << std::endl;
    TrajectoryContainer trajectories(theTrajBuilder.trajectories(tmpSeed));
    std::cout << "Trajectories found: " << trajectories.size() << std::endl;
    result.insert(
        result.end(), std::make_move_iterator(trajectories.begin()), std::make_move_iterator(trajectories.end()));
  }
  return result;
}

MuonTrajectoryBuilder::CandidateContainer ExhaustiveMuonTrajectoryBuilder::trajectories(const TrackCand&) {
  return CandidateContainer();
}

void ExhaustiveMuonTrajectoryBuilder::setEvent(const edm::Event& event) { theTrajBuilder.setEvent(event); }

void ExhaustiveMuonTrajectoryBuilder::clean(TrajectoryContainer& trajectories) const {
  // choose the one with the most hits, and the smallest chi-squared
  if (trajectories.empty()) {
    return;
  }
  int best_nhits = 0;
  unsigned best = 0;
  unsigned ntraj = trajectories.size();
  for (unsigned i = 0; i < ntraj; ++i) {
    int nhits = trajectories[i]->foundHits();
    if (nhits > best_nhits) {
      best_nhits = nhits;
      best = i;
    } else if (nhits == best_nhits && trajectories[i]->chiSquared() < trajectories[best]->chiSquared()) {
      best = i;
    }
  }
  TrajectoryContainer result;
  result.reserve(1);
  result.emplace_back(std::move(trajectories[best]));
  trajectories.swap(result);
}
