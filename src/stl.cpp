#include "stl.h"

#include "G4Exception.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4ThreeVector.hh"
#include <fstream>
#include <fmt/core.h>

struct STLTriangle {
  float normal[3];
  float vertex1[3];
  float vertex2[3];
  float vertex3[3];
};

G4TessellatedSolid* getTessalatedSolidFromSTL(
  const std::string& name, 
  const fs::path& stl_path,
  const G4ThreeVector& offset)
{
  std::ifstream is(stl_path.generic_u8string(), std::ios::binary );
  if(!is.good())
      throw std::runtime_error(fmt::format("Unable to read file {}", stl_path.generic_u8string()));

  auto volume_solid = new G4TessellatedSolid(name);

  uint8_t header[80];
  uint32_t triangle_count;
  is.read(reinterpret_cast<char*>(header), sizeof(header));
  is.read(reinterpret_cast<char*>(&triangle_count), sizeof(triangle_count));
  const float x_offset = static_cast<float>(offset.getX());
  const float y_offset = static_cast<float>(offset.getY());
  const float z_offset = static_cast<float>(offset.getZ());
  for(uint32_t i = 0; i < triangle_count; ++i) {
    STLTriangle triangle;
    uint16_t attribute_byte_count;
    is.read(reinterpret_cast<char*>(&triangle), sizeof(STLTriangle));
    is.read(reinterpret_cast<char*>(&attribute_byte_count), sizeof(attribute_byte_count));

    auto t = new G4TriangularFacet(
      G4ThreeVector(triangle.vertex1[0] + x_offset, triangle.vertex1[1] + y_offset, triangle.vertex1[2] + z_offset), 
      G4ThreeVector(triangle.vertex2[0] + x_offset, triangle.vertex2[1] + y_offset, triangle.vertex2[2] + z_offset), 
      G4ThreeVector(triangle.vertex3[0] + x_offset, triangle.vertex3[1] + y_offset, triangle.vertex3[2] + z_offset), 
      ABSOLUTE);
    volume_solid->AddFacet((G4VFacet *)t);
  }
  volume_solid->SetSolidClosed(true);
  if (volume_solid->GetNumberOfFacets() == 0) {
    throw std::runtime_error("The loaded mesh has 0 faces.");
  }
  return volume_solid;
}