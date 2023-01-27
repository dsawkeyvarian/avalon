
#include "TB02_CellFlux.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleDefinition.hh"
#include <stdlib.h>

TB02_CellFlux::TB02_CellFlux(G4String name, G4String particleName, G4double Elow, G4double Ehigh,G4int steps, G4int depth)
:G4VPrimitiveScorer(name,depth),HCID(-1),weighted(true)
{
  energy_low_limit = Elow;
  energyscale = Ehigh - Elow;
  stepsize = energyscale / steps;
  step = steps;
  G4cout << "Flux particles " << particleName << " and stepsize (MeV) " << stepsize << G4endl;
  G4cout << "Energy low limit " << energy_low_limit << " MeV"<< G4endl;
  if (particleName == "gamma") {
    pinfo = 22;
  } else if (particleName == "electron"){
    pinfo = 11;
  } else if (particleName == "O16"){
    pinfo = 1000080160;
  } else if (particleName == "proton"){
    pinfo = 2212;
  } else ;
  DefineUnitAndCategory();
  SetUnit("percm2");
  //verboseLevel = 10;
}

//TB02_CellFlux::TB02_CellFlux(G4String name, const G4String& unit, G4int depth)
//    :G4VPrimitiveScorer(name,depth),HCID(-1),weighted(true)
//{
  
  //G4ParticleDefinition* pd = G4ParticleTable::GetParticleTable()->FindParticle(particleName);
//  DefineUnitAndCategory();
//  SetUnit(unit);
//}

TB02_CellFlux::~TB02_CellFlux()
{;}

G4bool TB02_CellFlux::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double stepLength = aStep->GetStepLength();
  if ( stepLength == 0. ) return FALSE;
  G4int pd = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  if  (abs(pd) == pinfo) { 
	//G4String prosname = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
	//if (prosname != "compt" ) return FALSE;
    G4double totalEnergy = aStep->GetTrack()->GetTotalEnergy();
    for (G4int i=1; i < (step+1) ; i++){
      G4double e_low = energy_low_limit + (i-1)*stepsize;
      G4double e_high = energy_low_limit + i*stepsize;
      if ((e_low <= totalEnergy) && (totalEnergy < e_high)){ 
        G4int idx = ((G4TouchableHistory*)
        (aStep->GetPreStepPoint()->GetTouchable()))->GetReplicaNumber(indexDepth);
        G4double cubicVolume = ComputeVolume(aStep, idx);
                    
        G4double CellFlux = stepLength / cubicVolume;
        if (weighted) CellFlux *= aStep->GetPreStepPoint()->GetWeight(); 
        G4int index = GetIndex(aStep);
        EvtMap->add(index*step +i-1,CellFlux);
        return TRUE;
      }
    } 
    return FALSE;
  }
  else return FALSE;
}
    
void TB02_CellFlux::Initialize(G4HCofThisEvent* HCE)
{
EvtMap = new G4THitsMap<G4double>(detector->GetName(),
                                GetName());
  if ( HCID < 0 ) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID,EvtMap);
}

void TB02_CellFlux::EndOfEvent(G4HCofThisEvent*)
{;}

void TB02_CellFlux::clear(){
  EvtMap->clear();
}
 
void TB02_CellFlux::DrawAll()
{;}
 
void TB02_CellFlux::PrintAll()
{
 G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
 G4cout << " PrimitiveScorer " << GetName() <<G4endl; 
 G4cout << " Number of entries " << EvtMap->entries() << G4endl;
 std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
 for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
           << "  cell flux : " << *(itr->second)/GetUnitValue() 
           << " [" << GetUnit() << "]"
           << G4endl;
  }
}
 
void TB02_CellFlux::SetUnit(const G4String& unit)
{
  CheckAndSetUnit(unit,"Per Unit Surface");
}
 
void TB02_CellFlux::DefineUnitAndCategory(){
   // Per Unit Surface
  new G4UnitDefinition("percentimeter2","percm2","Per Unit Surface",(1./cm2));
  new G4UnitDefinition("permillimeter2","permm2","Per Unit Surface",(1./mm2));
  new G4UnitDefinition("permeter2","perm2","Per Unit Surface",(1./m2));
}
 
G4double TB02_CellFlux::ComputeVolume(G4Step* aStep, G4int idx){
 
  G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid* solid = 0;
  if(physParam)
  { // for parameterized volume
    if(idx<0)
    {
      G4ExceptionDescription ED;
      ED << "Incorrect replica number --- GetReplicaNumber : " << idx << G4endl;
       G4Exception("G4PSCellFlux::ComputeVolume","DetPS0001",JustWarning,ED);
    }
    solid = physParam->ComputeSolid(idx, physVol);
    solid->ComputeDimensions(physParam,idx,physVol);
  }
  else
  { // for ordinary volume
      solid = physVol->GetLogicalVolume()->GetSolid();
  }
  return solid->GetCubicVolume();
}
