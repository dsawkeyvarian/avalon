#ifndef TB02_CellFlux_h
#define TB02_CellFlux_h 1
#include "G4VPrimitiveScorer.hh"
#include "G4THitsMap.hh"

class G4VSolid;
 
class TB02_CellFlux : public G4VPrimitiveScorer
{
    public: // with description
        TB02_CellFlux(G4String name, G4String particleName ="gamma",G4double Elow =0.0, G4double Ehigh=1.0,G4int steps=1, G4int depth=0);
        //TB02_CellFlux(G4String name, const G4String& unit, G4int depth=0);
        virtual ~TB02_CellFlux();
    
        inline void Weighted(G4bool flg=true) { weighted = flg; }
      // Multiply track weight

    protected: // with description
        virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
        virtual G4double ComputeVolume(G4Step*, G4int idx);

   public: 
       virtual void Initialize(G4HCofThisEvent*);
       virtual void EndOfEvent(G4HCofThisEvent*);
       virtual void clear();
       virtual void DrawAll();
       virtual void PrintAll();
 
       virtual void SetUnit(const G4String& unit);    
 
   protected:
       virtual void DefineUnitAndCategory();

   private:
       G4int HCID;
       G4THitsMap<G4double>* EvtMap;
       G4bool  weighted;
       G4int pinfo;
       G4double energyscale;
       G4int step;
       G4double stepsize;
       G4double energy_low_limit;
};
#endif