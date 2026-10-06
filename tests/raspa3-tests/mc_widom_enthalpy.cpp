#include <gtest/gtest.h>

import std;

import archive;
import double3;
import units;
import atom;
import pseudo_atom;
import vdwparameters;
import forcefield;
import component;
import system;
import simulationbox;
import monte_carlo;
import property_widom;
import mc_moves_probabilities;
import mc_moves_move_types;

namespace
{
constexpr double epsilonK = 158.5;  // Lennard-Jones epsilon of the test particle [K]
constexpr double sigma = 3.72;      // [A]
constexpr double cutOff = 12.0;     // [A], truncated (not shifted), no tail corrections
constexpr double boxLength = 24.0;
constexpr double temperature = 150.0;

// Exact Boltzmann average of the test-particle energy over the box, by numerical quadrature on a
// regular grid: <U> = int U exp(-U/T) dr / int exp(-U/T) dr (energies in Kelvin).
double gridBoltzmannAverageEnergy(const std::vector<double3>& fixedPositions, double spacing)
{
  std::size_t n = static_cast<std::size_t>(boxLength / spacing);
  double h = boxLength / static_cast<double>(n);
  double sumWeight = 0.0;
  double sumWeightedEnergy = 0.0;
  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      for (std::size_t k = 0; k < n; ++k)
      {
        double3 p((static_cast<double>(i) + 0.5) * h, (static_cast<double>(j) + 0.5) * h,
                  (static_cast<double>(k) + 0.5) * h);
        double energy = 0.0;
        for (const double3& q : fixedPositions)
        {
          double3 d = p - q;
          d.x -= boxLength * std::round(d.x / boxLength);
          d.y -= boxLength * std::round(d.y / boxLength);
          d.z -= boxLength * std::round(d.z / boxLength);
          double rr = double3::dot(d, d);
          if (rr < cutOff * cutOff)
          {
            double s6 = std::pow(sigma * sigma / rr, 3);
            energy += 4.0 * epsilonK * (s6 * s6 - s6);
          }
        }
        if (energy / temperature > 500.0) continue;  // overlap: exp(-U/T) is zero to double precision
        double w = std::exp(-energy / temperature);
        sumWeight += w;
        sumWeightedEnergy += w * energy;
      }
    }
  }
  return sumWeightedEnergy / sumWeight;
}
}  // namespace

// Widom insertions into a fixed configuration of 60 Lennard-Jones particles: the energy-weighted Widom
// average <W dU>/<W> must reproduce the exact Boltzmann-averaged energy of a single test particle in
// that configuration, obtained by grid quadrature, and the enthalpy must equal it minus k_B T.
TEST(MC_WIDOM_ENTHALPY, energy_weighted_average_matches_grid_quadrature)
{
  ForceField forceField =
      ForceField({{"CH4", false, 16.04246, 0.0, 0.0, 6, false}}, {{epsilonK, sigma}},
                 ForceField::MixingRule::Lorentz_Berthelot, cutOff, cutOff, cutOff, false, false, false);

  MCMoveProbabilities probabilities = MCMoveProbabilities();
  probabilities.setProbability(Move::Types::Widom, 1.0);  // Widom only: the configuration never changes

  Component methane =
      Component(forceField, "methane", 190.564, 45599200, 0.01142, {Atom({0, 0, 0}, 0.0, 1.0, 0, 0, 0, false, false)},
                {}, {}, 5, 21, probabilities, std::nullopt, false);

  System system = System(forceField, SimulationBox(boxLength, boxLength, boxLength), false, temperature, 1e4, 1.0, {},
                         {methane}, {}, {60}, 5);

  std::size_t numberOfProductionCycles{5000};
  std::size_t numberOfInitializationCycles{0};
  std::size_t numberOfEquilibrationCycles{0};
  std::size_t printEvery{100000};
  std::size_t writeBinaryRestartEvery{1000000};
  std::size_t rescaleWangLandauEvery{100000};
  std::size_t optimizeMCMovesEvery{100000};
  MonteCarlo mc = MonteCarlo({numberOfProductionCycles, 0, numberOfInitializationCycles, numberOfEquilibrationCycles,
                              printEvery, writeBinaryRestartEvery, rescaleWangLandauEvery, optimizeMCMovesEvery},
                             {system}, 42uz, 5, false);
  mc.run();

  const System& s = mc.systems.front();
  std::vector<double3> positions;
  for (const Atom& atom : s.spanOfMoleculeAtoms()) positions.push_back(atom.position);
  ASSERT_EQ(positions.size(), 60uz);

  double reference = gridBoltzmannAverageEnergy(positions, 0.1);

  const PropertyWidom& widom = s.components.front().averageRosenbluthWeights;
  auto [insertionEnergy, insertionEnergyError] = widom.insertionEnergyResult();
  insertionEnergy *= Units::EnergyToKelvin;
  insertionEnergyError *= Units::EnergyToKelvin;

  // The configuration must give a non-trivial reference (attractive by more than k_B T on average).
  EXPECT_LT(reference, -temperature);
  EXPECT_TRUE(std::isfinite(insertionEnergy));
  EXPECT_GT(insertionEnergyError, 0.0);
  EXPECT_NEAR(insertionEnergy, reference, std::max(3.0 * insertionEnergyError, 0.02 * std::abs(reference)));

  auto [enthalpy, enthalpyError] = widom.enthalpyResult(s.beta);
  // k_B T in internal units converted back to Kelvin agrees with T to ~1e-6 relative (unit constants)
  EXPECT_NEAR(Units::EnergyToKelvin * enthalpy, insertionEnergy - temperature, 1e-3);
  EXPECT_NEAR(Units::EnergyToKelvin * enthalpyError, insertionEnergyError, 1e-8 * std::abs(insertionEnergyError));
}

// The energy-weighted channel survives a binary-restart round trip.
TEST(MC_WIDOM_ENTHALPY, restart_round_trip)
{
  PropertyWidom widom(5);
  for (std::size_t block = 0; block < 5; ++block)
  {
    for (std::size_t i = 0; i < 10; ++i)
    {
      double w = 0.1 * static_cast<double>(i + 1 + block);
      widom.addWidomSample(block, w, -3.0 + 0.2 * static_cast<double>(i), 0, 1000.0, 1.0);
    }
  }

  std::filesystem::path path = std::filesystem::temp_directory_path() / "raspa3_widom_enthalpy_restart.bin";
  {
    std::ofstream stream(path, std::ios::binary);
    Archive<std::ofstream> archive(stream);
    archive << widom;
  }
  PropertyWidom restored;
  {
    std::ifstream stream(path, std::ios::binary);
    Archive<std::ifstream> archive(stream);
    archive >> restored;
  }
  std::filesystem::remove(path);

  auto [mean, error] = widom.insertionEnergyResult();
  auto [restoredMean, restoredError] = restored.insertionEnergyResult();
  EXPECT_DOUBLE_EQ(mean, restoredMean);
  EXPECT_DOUBLE_EQ(error, restoredError);
}
