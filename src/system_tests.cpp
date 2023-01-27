#include "external/doctest.h"

#define JSON_DIAGNOSTICS 1
#include "external/json.hpp"
using json = nlohmann::json;

#include <filesystem>
namespace fs = std::filesystem;
#include "linac.h"
#include "logger.h"
#include "simulation.h"
#include "test_utils.h"
#include "test.h"

double calculateMeanEnergy(const PhaseSpace& phsp)
{
  size_t num_particles = phsp.getNumParticles();
  size_t summed_particles = 0;
  double sum = 0.0;
  for(size_t i = 0; i < num_particles; ++i) {
    auto p = phsp.getParticle(i);
    bool isGoingPlusZ = p->weight > 0.0;
    if(isGoingPlusZ) //We want particles travelling towards phantom, not scattering backwards
      continue;
    summed_particles += 1;
    sum += p->energy;
  }
  return sum/static_cast<double>(summed_particles);
}

TEST_SUITE_BEGIN("System");

TEST_CASE("UniformPointSource"
  * doctest::description("Test uniform point source phasespace and dose"))
{
  const std::string json_input = R"(
{
    "physics_model": "QGSP_BIC", //Inaccurate but fast to initialize for testing
    "treatment_head": {
        "type": "no_head",
        "source": "uniform_point"
    },
    "uniform_point_source": {
        "particle": "PHOTON",
        "z": [1000.0, "mm"],
        "spectrum": {
            "unit": "MeV",
            "energies": [0.00000, 0.07500, 0.22500, 0.37500, 0.52500, 0.67500, 0.82500, 0.97500, 1.12500, 1.27500, 1.42500, 1.57500, 1.72500, 1.87500, 2.02500, 2.17500, 2.32499, 2.47499, 2.62499, 2.77499, 2.92499, 3.07499, 3.22499, 3.37499, 3.52499, 3.67499, 3.82499, 3.97499, 4.12499, 4.27499, 4.42499, 4.57499, 4.72499, 4.87499, 5.02499, 5.17499, 5.32499, 5.47499, 5.62499, 5.77499, 5.92499, 5.99999, 6.07499],
            "intensities": [0.00000, 0.22270, 1.38521, 1.49894, 1.20472, 0.94601, 0.74626, 0.61317, 0.51813, 0.43576, 0.37460, 0.33335, 0.29606, 0.26020, 0.22824, 0.21108, 0.19138, 0.16633, 0.15929, 0.13970, 0.12627, 0.11495, 0.10895, 0.09709, 0.08718, 0.08227, 0.07822, 0.06882, 0.06269, 0.06091, 0.05453, 0.05016, 0.04290, 0.04609, 0.04031, 0.02844, 0.03023, 0.02151, 0.01822, 0.01866, 0.01341, 0.01063, 0.00000]
        },
        "field_size": [25, "mm"]
    },
    "dose_grids": [
        {
            "filename": "out.dose",
            "center_position": [0.0, 0.0, -15.0, "cm"],
            "size": [1.0, 1.0, 30.0, "cm"],
            "num_voxels": [1, 1, 60]
        }
    ],
    "phase_space_recorders": [
        {
            "type": "plane",
            "phsp_filename": "out.phsp",
            "center_position": [0.0, 0.0, 0.1, "cm"],
            "size": [200.0, 200.0, 0.01, "cm"],
            "kill_particles": false
        }
    ],
    "phantom": {
        "type": "water_box",
        "center_position": [0.0, 0.0, -15.0, "cm"],
        "size": [30.0, 30.0, 30.0, "cm"]
    },
    "simulation_points": [
        {
            "weight": 1.0,
            "angles": {"unit": "deg", "gantry": 0.0, "collimator": 0.0}
        }
    ]
}
)";
  auto root = getVirtualLinacRoot();
  auto out_dir = getOutputFolder();
  CommandLineAndEnvironmentArguments args;
  args.virtuallinac_root = root;
  args.num_particles = 2000000;
  args.input_json = json_input;
  args.output_directory = out_dir;
  args.dry_run = false;
  COutAndFileLogger log;
  CHECK(doLinacSimulation(args, &runG4Simulation, &log) == 0);

  auto dose_path = out_dir / "out.dose";
  auto phsp_path = out_dir / "out.phsp";
  const double expected_mean_energy = 1.265324050558101;
  std::vector<double> expected_pdd = { 1.226900e-13, 1.286862e-13, 1.318843e-13, 1.362851e-13, 1.409652e-13, 1.457558e-13, 1.499649e-13, 1.544982e-13, 1.613125e-13, 1.657265e-13, 1.706243e-13, 1.760191e-13, 1.832583e-13, 1.893167e-13, 1.947981e-13, 1.998765e-13, 2.072935e-13, 2.144586e-13, 2.194103e-13, 2.281435e-13, 2.350508e-13, 2.428037e-13, 2.486735e-13, 2.574436e-13, 2.677929e-13, 2.764909e-13, 2.875417e-13, 2.955126e-13, 3.048334e-13, 3.154172e-13, 3.256514e-13, 3.369810e-13, 3.492808e-13, 3.620325e-13, 3.750503e-13, 3.881585e-13, 4.003133e-13, 4.145702e-13, 4.297368e-13, 4.437233e-13, 4.610352e-13, 4.752642e-13, 4.918447e-13, 5.094487e-13, 5.294235e-13, 5.503553e-13, 5.685581e-13, 5.900446e-13, 6.081120e-13, 6.318013e-13, 6.518756e-13, 6.748198e-13, 6.987591e-13, 7.213659e-13, 7.476569e-13, 7.703234e-13, 7.920803e-13, 7.943376e-13, 7.469919e-13, 5.054220e-13 };

  Dose dose(dose_path);
  CHECK(expected_pdd.size() == dose.getTotalVoxelCount());
  for(size_t i = 0; i < expected_pdd.size(); ++i) {
    double rel_diff = std::abs(expected_pdd[i] - dose.m_doses[i])/expected_pdd[i];
    CHECK(rel_diff < 0.09); //Relatively loose criterion due to anount of noise, only tests that dose grid is working and roughly in correct place
  }
  PhaseSpace phsp(phsp_path);
  double mean = calculateMeanEnergy(phsp);
  double relative_mean_error = std::abs((expected_mean_energy - mean)/expected_mean_energy);
  CHECK_MESSAGE(relative_mean_error < 0.001, "Relative error of spectrum's mean energy");
}

TEST_CASE("WaterCylinderAndGantry"
  * doctest::description("Test WaterCylinder phantom and gantry rotation with dose"))
{
  const std::string json_input = R"(
{
	//Rotates the source around Y-axis (gantry angle) from 0 to 90 degrees
    //results visible on XZ-dose plane
	"num_particles": 2000000,
    "physics_model": "QGSP_BIC",
    "treatment_head": {
        "type": "no_head",
        "source": "uniform_point"
    },
    "uniform_point_source": {
        "particle": "PHOTON",
        "z": [1000.0, "mm"],
        "spectrum": {
            "unit": "MeV",
            "energies": [0.00000, 0.07500, 0.22500, 0.37500, 0.52500, 0.67500, 0.82500, 0.97500, 1.12500, 1.27500, 1.42500, 1.57500, 1.72500, 1.87500, 2.02500, 2.17500, 2.32499, 2.47499, 2.62499, 2.77499, 2.92499, 3.07499, 3.22499, 3.37499, 3.52499, 3.67499, 3.82499, 3.97499, 4.12499, 4.27499, 4.42499, 4.57499, 4.72499, 4.87499, 5.02499, 5.17499, 5.32499, 5.47499, 5.62499, 5.77499, 5.92499, 5.99999, 6.07499],
            "intensities": [0.00000, 0.22270, 1.38521, 1.49894, 1.20472, 0.94601, 0.74626, 0.61317, 0.51813, 0.43576, 0.37460, 0.33335, 0.29606, 0.26020, 0.22824, 0.21108, 0.19138, 0.16633, 0.15929, 0.13970, 0.12627, 0.11495, 0.10895, 0.09709, 0.08718, 0.08227, 0.07822, 0.06882, 0.06269, 0.06091, 0.05453, 0.05016, 0.04290, 0.04609, 0.04031, 0.02844, 0.03023, 0.02151, 0.01822, 0.01866, 0.01341, 0.01063, 0.00000]
        },
        "field_size": [25, "mm"]
    },
    "phantom": {
        "type": "water_cylinder",
        "center_position": [0.0, 0.0, 0.0, "cm"],
        "height_axis": "y",
        "height": [35.0, "cm"],
        "radius": [10.0, "cm"]
    },
    "dose_grids": [
        {
            "filename": "out.dose",
            "center_position": [0.0, 0.0, 0.0, "cm"],
            "size": [20.0, 1.0, 20.0, "cm"],
            "num_voxels": [20, 1, 20]
        }
    ],
    "simulation_points": [
        {
            "weight": 0.5,
            "angles": {"unit": "deg", "gantry": 0.0, "collimator": 0.0}
        },
        {
            "weight": 0.5,
            "angles": {"unit": "deg", "gantry": 90.0, "collimator": 0.0}
        }
    ]
}
)";
  auto root = getVirtualLinacRoot();
  auto out_dir = getOutputFolder();
  CommandLineAndEnvironmentArguments args;
  args.virtuallinac_root = root;
  args.num_particles = 2000000;
  args.input_json = json_input;
  args.output_directory = out_dir;
  args.dry_run = false;
  COutAndFileLogger log;
  CHECK(doLinacSimulation(args, &runG4Simulation, &log) == 0);

  auto dose_path = out_dir / "out.dose";
  Dose dose(dose_path);
  double expected_center_dose = 5.50e-13;
  std::vector<double> expected_slice_dose = { 1.36889e-13, 1.47207e-13, 1.59090e-13, 1.68104e-13, 1.80593e-13, 1.94805e-13, 2.09743e-13, 2.27784e-13, 3.17399e-13, 5.29472e-13, 5.50337e-13, 3.73393e-13, 3.22880e-13, 3.44406e-13, 3.67796e-13, 3.93242e-13, 4.23010e-13, 4.49739e-13, 4.73547e-13, 3.72682e-13 };
  size_t center_x_idx = 10;
  size_t center_z_idx = 10;
  double center_dose = dose.getPointDose(center_x_idx, 0, center_z_idx);
  double relative_err_center_dose = std::abs(center_dose - expected_center_dose)/expected_center_dose;
  CHECK(relative_err_center_dose < 0.05);
  auto x_slice = dose.getSlice1D(Axis::X, 0, center_z_idx);
  auto z_slice = dose.getSlice1D(Axis::Z, center_x_idx, 0);
  const size_t num_pts = expected_slice_dose.size();
  REQUIRE(x_slice.size() == expected_slice_dose.size());
  REQUIRE(z_slice.size() == expected_slice_dose.size());
  for(size_t i = 0; i < num_pts; ++i) {
    double relative_err_x = std::abs(x_slice[i] - expected_slice_dose[i])/expected_slice_dose[i];
    double relative_err_z = std::abs(z_slice[i] - expected_slice_dose[i])/expected_slice_dose[i];
    CHECK(relative_err_x < 0.10); //Arbitary/Loose criterion due to speed
    CHECK(relative_err_z < 0.10);
  }
}

TEST_CASE("CollimatorAndGantry"
  * doctest::description("Test dose into water box with gantry and collimator angles"))
{
  const std::string json_input = R"(
{
    "physics_model": "QGSP_BIC",
    "treatment_head": {
        "type": "no_head",
        "source": "uniform_point"
    },
    "uniform_point_source": {
        "particle": "PHOTON",
        "z": [1000.0, "mm"],
        "spectrum": { 
            "unit": "MeV",
            "energies": [0.00000, 0.07500, 0.22500, 0.37500, 0.52500, 0.67500, 0.82500, 0.97500, 1.12500, 1.27500, 1.42500, 1.57500, 1.72500, 1.87500, 2.02500, 2.17500, 2.32499, 2.47499, 2.62499, 2.77499, 2.92499, 3.07499, 3.22499, 3.37499, 3.52499, 3.67499, 3.82499, 3.97499, 4.12499, 4.27499, 4.42499, 4.57499, 4.72499, 4.87499, 5.02499, 5.17499, 5.32499, 5.47499, 5.62499, 5.77499, 5.92499, 5.99999, 6.07499],
            "intensities": [0.00000, 0.22270, 1.38521, 1.49894, 1.20472, 0.94601, 0.74626, 0.61317, 0.51813, 0.43576, 0.37460, 0.33335, 0.29606, 0.26020, 0.22824, 0.21108, 0.19138, 0.16633, 0.15929, 0.13970, 0.12627, 0.11495, 0.10895, 0.09709, 0.08718, 0.08227, 0.07822, 0.06882, 0.06269, 0.06091, 0.05453, 0.05016, 0.04290, 0.04609, 0.04031, 0.02844, 0.03023, 0.02151, 0.01822, 0.01866, 0.01341, 0.01063, 0.00000]
        },
        "field_size": {
            "unit": "mm",
            "x1": -5.0, "x2": 5.0,
            "y1": -20.0, "y2": 20.0
        }
    },
    "phantom": {
        "type": "water_box",
        "center_position": [0.0, 0.0, 0.0, "cm"],
        "size": [20.0, 10.0, 10.0, "cm"]
    },
    "dose_grids": [
        {
            "filename": "out.dose",
            "center_position": [0.0, 0.0, 0.0, "cm"],
            "size": [1.0, 10.0, 10.0, "cm"],
            "num_voxels": [1, 20, 20]
        }
    ],
    "simulation_points": [
        {
            "weight": 1.0,
            "angles": {"unit": "deg", "gantry": 90.0, "collimator": 45.0}
        }
    ]
}
)";
  auto root = getVirtualLinacRoot();
  auto out_dir = getOutputFolder();
  CommandLineAndEnvironmentArguments args;
  args.virtuallinac_root = root;
  args.num_particles = 2000000;
  args.input_json = json_input;
  args.output_directory = out_dir;
  args.dry_run = false;
  COutAndFileLogger log;
  CHECK(doLinacSimulation(args, &runG4Simulation, &log) == 0);

  auto dose_path = out_dir / "out.dose";
  Dose dose(dose_path);
  //YZ-plane with Z-slice on each row 
  std::vector<double> expected_dose = {
    1.549880e-15, 1.868747e-15, 2.177388e-15, 2.391537e-15, 2.696297e-15, 2.963480e-15, 3.156509e-15, 3.268299e-15, 3.476847e-15, 3.431342e-15, 3.417806e-15, 3.224779e-15, 3.023701e-15, 2.809650e-15, 2.577416e-15, 2.292631e-15, 2.062813e-15, 1.860586e-15, 1.683909e-15, 1.360874e-15,
    1.900753e-15, 2.192767e-15, 2.489973e-15, 2.918105e-15, 3.306556e-15, 3.734347e-15, 3.940136e-15, 4.181001e-15, 4.330415e-15, 4.353174e-15, 4.397743e-15, 3.979295e-15, 3.825426e-15, 3.394854e-15, 3.117764e-15, 2.726172e-15, 2.403017e-15, 2.192394e-15, 1.893400e-15, 1.606729e-15,
    2.128407e-15, 2.573169e-15, 2.882555e-15, 3.485149e-15, 4.078819e-15, 4.688265e-15, 5.118821e-15, 5.505049e-15, 5.631576e-15, 5.618009e-15, 5.465691e-15, 5.016367e-15, 4.687604e-15, 4.155971e-15, 3.739375e-15, 3.326581e-15, 2.905197e-15, 2.513607e-15, 2.138861e-15, 1.832679e-15,
    2.407183e-15, 2.905359e-15, 3.421708e-15, 4.180076e-15, 4.902377e-15, 5.775702e-15, 6.570691e-15, 7.204457e-15, 7.335122e-15, 7.302898e-15, 7.097466e-15, 6.420927e-15, 5.744177e-15, 5.161761e-15, 4.576566e-15, 3.900799e-15, 3.449465e-15, 2.889302e-15, 2.378894e-15, 2.061921e-15,
    2.738079e-15, 3.372267e-15, 4.110399e-15, 4.943707e-15, 6.000146e-15, 7.279659e-15, 8.660282e-15, 9.921540e-15, 1.038190e-14, 9.690017e-15, 9.208477e-15, 8.301098e-15, 7.305138e-15, 6.353240e-15, 5.429554e-15, 4.587623e-15, 3.880629e-15, 3.380102e-15, 2.774715e-15, 2.261460e-15,
    2.959538e-15, 3.676659e-15, 4.466075e-15, 5.740829e-15, 7.442982e-15, 9.711098e-15, 1.331557e-14, 1.890846e-14, 1.886889e-14, 1.489243e-14, 1.239844e-14, 1.062761e-14, 9.278152e-15, 7.686506e-15, 6.434657e-15, 5.411412e-15, 4.508686e-15, 3.772080e-15, 3.141664e-15, 2.528184e-15,
    3.283019e-15, 3.925484e-15, 5.047829e-15, 6.423734e-15, 8.651931e-15, 1.348763e-14, 3.079549e-14, 1.861190e-13, 1.234679e-13, 2.985065e-14, 1.844988e-14, 1.422624e-14, 1.156711e-14, 9.397033e-15, 7.627523e-15, 6.263188e-15, 5.176571e-15, 4.281586e-15, 3.374407e-15, 2.739183e-15,
    3.320009e-15, 4.090034e-15, 5.489881e-15, 7.213968e-15, 9.810158e-15, 1.862455e-14, 1.859742e-13, 6.826069e-13, 6.188993e-13, 1.442254e-13, 3.379655e-14, 1.970827e-14, 1.457099e-14, 1.164220e-14, 9.188610e-15, 7.290939e-15, 5.765504e-15, 4.690714e-15, 3.708035e-15, 2.872013e-15,
    3.394096e-15, 4.227173e-15, 5.563248e-15, 7.460991e-15, 1.057266e-14, 1.880270e-14, 1.228378e-13, 6.180164e-13, 7.821434e-13, 6.297949e-13, 1.474646e-13, 3.504620e-14, 2.000015e-14, 1.398445e-14, 1.064625e-14, 8.221503e-15, 6.481591e-15, 5.091230e-15, 4.036540e-15, 3.205754e-15,
    3.422572e-15, 4.396289e-15, 5.616747e-15, 7.305909e-15, 1.012207e-14, 1.489999e-14, 2.986548e-14, 1.447000e-13, 6.297452e-13, 7.828227e-13, 6.321390e-13, 1.471005e-13, 3.395555e-14, 1.858067e-14, 1.274415e-14, 9.232946e-15, 7.042466e-15, 5.448908e-15, 4.393496e-15, 3.299836e-15,
    3.362929e-15, 4.242827e-15, 5.433437e-15, 6.860658e-15, 9.149371e-15, 1.247787e-14, 1.832242e-14, 3.383882e-14, 1.472008e-13, 6.293109e-13, 7.821849e-13, 6.283391e-13, 1.440031e-13, 2.977893e-14, 1.488010e-14, 9.876764e-15, 7.391632e-15, 5.718590e-15, 4.460364e-15, 3.398774e-15,
    3.186422e-15, 4.000894e-15, 5.011028e-15, 6.399873e-15, 8.305839e-15, 1.064006e-14, 1.399997e-14, 1.991349e-14, 3.482984e-14, 1.466170e-13, 6.283786e-13, 7.770253e-13, 6.203307e-13, 1.226561e-13, 1.852027e-14, 1.048221e-14, 7.429070e-15, 5.577039e-15, 4.352980e-15, 3.390351e-15,
    3.025938e-15, 3.653261e-15, 4.669301e-15, 5.747186e-15, 7.185038e-15, 9.052240e-15, 1.155166e-14, 1.453477e-14, 1.982545e-14, 3.353470e-14, 1.446167e-13, 6.182953e-13, 6.833852e-13, 1.853605e-13, 1.888342e-14, 9.822484e-15, 6.917426e-15, 5.300270e-15, 4.293381e-15, 3.313113e-15,
    2.830220e-15, 3.436585e-15, 4.286734e-15, 5.116959e-15, 6.394434e-15, 7.683826e-15, 9.437472e-15, 1.152375e-14, 1.440433e-14, 1.820003e-14, 2.990396e-14, 1.231563e-13, 1.860400e-13, 3.076658e-14, 1.328980e-14, 8.733683e-15, 6.526033e-15, 4.978404e-15, 3.913684e-15, 3.162386e-15,
    2.529521e-15, 3.102223e-15, 3.663489e-15, 4.518928e-15, 5.426137e-15, 6.417170e-15, 7.661090e-15, 9.011296e-15, 1.052560e-14, 1.247162e-14, 1.485668e-14, 1.864840e-14, 1.867550e-14, 1.349048e-14, 9.814851e-15, 7.425544e-15, 5.666654e-15, 4.590396e-15, 3.624742e-15, 2.939183e-15,
    2.345514e-15, 2.834199e-15, 3.226772e-15, 3.948647e-15, 4.556548e-15, 5.488983e-15, 6.293976e-15, 7.279549e-15, 8.275092e-15, 9.357016e-15, 1.002872e-14, 1.058603e-14, 1.021052e-14, 8.821714e-15, 7.464140e-15, 6.138332e-15, 5.020229e-15, 4.164014e-15, 3.296193e-15, 2.550786e-15,
    1.996131e-15, 2.413364e-15, 2.906872e-15, 3.402019e-15, 3.921773e-15, 4.518753e-15, 5.016893e-15, 5.750184e-15, 6.338501e-15, 6.988912e-15, 7.519259e-15, 7.442787e-15, 6.961621e-15, 6.288288e-15, 5.655677e-15, 4.965761e-15, 4.248054e-15, 3.519352e-15, 2.944735e-15, 2.477601e-15,
    1.743661e-15, 2.167203e-15, 2.572638e-15, 2.874798e-15, 3.247574e-15, 3.675433e-15, 4.219741e-15, 4.789285e-15, 5.090961e-15, 5.439455e-15, 5.690678e-15, 5.586485e-15, 5.400467e-15, 5.140181e-15, 4.732021e-15, 3.999850e-15, 3.478098e-15, 3.008694e-15, 2.568491e-15, 2.136709e-15,
    1.627963e-15, 1.941321e-15, 2.219560e-15, 2.450487e-15, 2.761036e-15, 3.142898e-15, 3.525155e-15, 3.817899e-15, 3.999179e-15, 4.189602e-15, 4.312994e-15, 4.484723e-15, 4.371357e-15, 4.043887e-15, 3.717357e-15, 3.190780e-15, 2.891582e-15, 2.495830e-15, 2.156409e-15, 1.821577e-15,
    1.386365e-15, 1.657862e-15, 1.838565e-15, 2.098135e-15, 2.304468e-15, 2.562286e-15, 2.766384e-15, 2.954199e-15, 3.140546e-15, 3.352550e-15, 3.450517e-15, 3.449328e-15, 3.217795e-15, 3.109277e-15, 2.906213e-15, 2.665865e-15, 2.390481e-15, 2.104569e-15, 1.902619e-15, 1.574330e-15 };
  const size_t num_pts = expected_dose.size();
  REQUIRE(dose.getTotalVoxelCount() == num_pts);
  double max_expected_dose = 0.0;
  for(size_t i = 0; i < num_pts; ++i)
    max_expected_dose = std::max(max_expected_dose, expected_dose[i]);
  std::vector<double> rel_errors;
  for(size_t i = 0; i < num_pts; ++i) {
    if(expected_dose[i]/max_expected_dose < 0.2)
      continue;
    double rel_error = std::abs(dose.m_doses[i] - expected_dose[i])/expected_dose[i];
    CHECK(rel_error < 0.1);
  }
}

TEST_SUITE_END();
