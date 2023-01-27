#ifndef TB02_CellFlux3D_h
#define TB02_CellFlux3D_h 1


#include "TB02_CellFlux.hh"

class TB02_CellFlux3D : public TB02_CellFlux
{
    public: // with description
        TB02_CellFlux3D(G4String name, G4String particleName, G4double Elow, G4double Ehigh,G4int steps, G4int ni=1, G4int nj=1, G4int nk=1, G4int depi=2, G4int depj=1, G4int depk=0);
        //TB02_CellFlux3D(G4String name, const G4String& unit, G4int ni, G4int nj, G4int nk, G4int depi, G4int depj, G4int depk);
        virtual ~TB02_CellFlux3D();

    protected: // with description
        virtual G4int GetIndex(G4Step*);

    private:
        G4int fDepthi, fDepthj, fDepthk;
};
#endif