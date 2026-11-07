#include "DataFormats/VertexReco/interface/Vertex.h"

#include "PhysicsTools/RDsNano/plugins/ParametrizedEngine/OAEParametrizedMagneticField.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexPrimitives/interface/CachingVertex.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include <memory>


class PVRefitter {

  public:
    PVRefitter() : paramField(std::make_unique<OAEParametrizedMagneticField>("3_8T")) {};
    virtual ~PVRefitter() {};

    // constructed from reco::TrackRef
    reco::TransientTrack getTransientTrack(const reco::TrackRef& trackRef) const {    
      reco::TransientTrack transientTrack(trackRef, paramField.get());
      return transientTrack;
    }

    // constructed from reco::Track
    reco::TransientTrack getTransientTrack(const reco::Track& track) const {    
      reco::TransientTrack transientTrack(track, paramField.get());
      return transientTrack;
    }


    reco::Vertex fit(const std::vector<reco::Track>& tracks, const reco::BeamSpot& bs) const {
    //TransientVertex fit(const std::vector<reco::Track>& tracks, const reco::BeamSpot& bs){
    //CachingVertex<5> fit(const std::vector<reco::Track>& tracks, const reco::BeamSpot& bs){
    
      std::vector<reco::TransientTrack> transient_tracks;
      transient_tracks.reserve(tracks.size());
      
      for (reco::TrackCollection::const_iterator itrk=tracks.begin(); itrk!=tracks.end(); itrk++){
        reco::TransientTrack tmpTransientTrack = getTransientTrack(*itrk); 
        tmpTransientTrack.setBeamSpot(bs);
        transient_tracks.push_back(tmpTransientTrack);
      }
      
      AdaptiveVertexFitter fitter;
      TransientVertex vertex = fitter.vertex(transient_tracks, bs);
      
//       return vertex; 
      return reco::Vertex(vertex); 
    }

  private:
    // owned here: the raw "new" without delete leaked one field object per PVRefitter instance
    std::unique_ptr<OAEParametrizedMagneticField> paramField;

};
