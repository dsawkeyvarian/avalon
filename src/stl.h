#include <string>
#include <filesystem>
namespace fs = std::filesystem;

class G4TessellatedSolid;
namespace CLHEP {
  class Hep3Vector;
}

G4TessellatedSolid* getTessalatedSolidFromSTL(
  const std::string& name, 
  const fs::path& stl_path,
  const CLHEP::Hep3Vector& offset
);