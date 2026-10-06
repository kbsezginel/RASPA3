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
 * (framework, other molecules, external field, Ewald Fourier, tail and polarization; intramolecular
 * terms excluded). Weighted by 'rosenbluthWeight' it gives the Boltzmann average of the insertion
 * energy, <W dU>/<W>.
 */
struct WidomInsertion
{
  double rosenbluthWeight{0.0};  ///< Rosenbluth weight normalized by the ideal-gas Rosenbluth weight.
  double insertionEnergy{0.0};   ///< Intermolecular energy of the inserted configuration.
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
