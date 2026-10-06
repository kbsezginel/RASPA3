module;

export module mc_moves_widom;

import std;

import double3;
import randomnumbers;
import running_energy;
import atom;
import system;

export namespace MC_Moves
{
/**
 * \brief Result of a single Widom test insertion.
 *
 * 'insertionEnergy' is the intermolecular energy of the configuration selected by the CBMC growth
 * (framework, other molecules, external field, Ewald Fourier, tail and polarization) and 'intraEnergy'
 * its intramolecular energy. Weighted by 'rosenbluthWeight' they give Boltzmann averages over the
 * inserted molecule, <W A>/<W>. For flexible components an isolated (ideal-gas) molecule is grown as
 * well; 'idealGasWeight' and 'idealGasIntraEnergy' give the ideal-gas reference <U_intra>_IG the same
 * way. Rigid components report intraEnergy = 0, idealGasWeight = 1, idealGasIntraEnergy = 0.
 */
struct WidomInsertion
{
  double rosenbluthWeight{0.0};     ///< Rosenbluth weight normalized by the ideal-gas Rosenbluth weight.
  double insertionEnergy{0.0};      ///< Intermolecular energy of the inserted configuration.
  double intraEnergy{0.0};          ///< Intramolecular energy of the inserted configuration.
  double idealGasWeight{1.0};       ///< Rosenbluth weight of the isolated (ideal-gas) growth.
  double idealGasIntraEnergy{0.0};  ///< Intramolecular energy of the isolated (ideal-gas) growth.
};

/**
 * \brief Performs a Widom insertion move for the specified component.
 *
 * Attempts to insert a molecule of the selected component into the system using
 * Configurational Bias Monte Carlo (CBMC) method. Calculates the energy differences
 * and computes the insertion weight used for chemical potential estimation.
 *
 * \param random Reference to the random number generator.
 * \param system Reference to the simulation system.
 * \param selectedComponent Index of the component to perform the Widom move on.
 * \return The Widom insertion weight and the intermolecular energy of the inserted configuration
 *         (both zero if the growth failed).
 */
WidomInsertion WidomMove(RandomNumber& random, System& system, std::size_t selectedComponent);
}  // namespace MC_Moves
