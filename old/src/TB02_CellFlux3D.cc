//* http://www.apc.univ-paris7.fr/~franco/g4doxy/html/G4PSCellFlux_8cc-source.html
// Generally this is the average of the distance that each particle travels through a cell 
//    cm / cm^3 = particles/cm^2
//    Tally by photon and electron energy ranges and you have the exact result from XB. 
//    (XB when uses a flux-to-dose conversion factor to turn flux into Gy. 
//    If the fluxes agree, then there is a difference in how Geant tallies dose.

#include "TB02_CellFlux3D.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TB02_CellFlux3D::TB02_CellFlux3D(G4String name, G4String particleName,G4double Elow, G4double Ehigh,G4int steps, G4int ni, G4int nj, G4int nk, G4int depi, G4int depj, G4int depk)
:TB02_CellFlux(name, particleName,Elow,Ehigh,steps), fDepthi(depi),fDepthj(depj),fDepthk(depk)
                                //G4int depi, G4int depj, G4int depk)
{
    fNi=ni;
    fNj=nj;
    fNk=nk;
}

//TB02_CellFlux3D::TB02_CellFlux3D(G4String name, G4String& unit, G4int ni, G4int nj, G4int nk, G4int depi, G4int depj, G4int depk)
//    :TB02_CellFlux(name), fDepthi(depi),fDepthj(depj),fDepthk(depk)
                                //G4int depi, G4int depj, G4int depk)
//{
//    fNi=ni;
//    fNj=nj;
//    fNk=nk;
//    SetUnit(unit);
//}


TB02_CellFlux3D::~TB02_CellFlux3D()
{;} 

G4int TB02_CellFlux3D::GetIndex(G4Step* aStep)
{
    const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int i = touchable->GetReplicaNumber(fDepthi);
    G4int j = touchable->GetReplicaNumber(fDepthj);
    G4int k = touchable->GetReplicaNumber(fDepthk);

    return i*fNj*fNk+j*fNk+k;
}