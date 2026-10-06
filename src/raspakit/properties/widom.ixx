module;

export module property_widom;

import std;

import archive;
import int3;
import averages;
export import property_block_average;

export struct WidomData
{
  WidomData():
    total(0.0),
    excess(0.0),
    idealGas(0.0)
  {
  };

  WidomData(double total, double excess, double idealGas):
    total(total),
    excess(excess),
    idealGas(idealGas)
  {
  };

  inline WidomData& operator+=(const WidomData& b)
  {
    total += b.total;
    excess += b.excess;
    idealGas += b.idealGas;

    return *this;
  }

  std::uint64_t versionNumber{1};

  double total{};
  double excess{};
  double idealGas{};

  friend Archive<std::ofstream>& operator<<(Archive<std::ofstream>& archive, const WidomData& l);
  friend Archive<std::ifstream>& operator>>(Archive<std::ifstream>& archive, WidomData& l);
};

export inline WidomData operator+(const WidomData& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a.total + b.total;
  m.excess = a.excess + b.excess;
  m.idealGas = a.idealGas + b.idealGas;

  return m;
}

export inline WidomData operator-(const WidomData& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a.total - b.total;
  m.excess = a.excess - b.excess;
  m.idealGas = a.idealGas - b.idealGas;

  return m;
}

export inline WidomData operator*(const WidomData& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a.total * b.total;
  m.excess = a.excess * b.excess;
  m.idealGas = a.idealGas * b.idealGas;

  return m;
}

export inline WidomData operator*(const double& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a * b.total;
  m.excess = a * b.excess;
  m.idealGas = a * b.idealGas;

  return m;
}

export inline WidomData operator/(const WidomData& a, const double& b)
{
  WidomData m{}; 

  double inv_b = 1.0 / b;
  m.total = inv_b * a.total;
  m.excess = inv_b * a.excess;
  m.idealGas = inv_b * a.idealGas;

  return m;
}

export inline WidomData operator/(const WidomData& a, const std::array<double,3>& b)
{
  WidomData m{}; 

  m.total = a.total / b[0];
  m.excess = a.excess / b[1];
  m.idealGas = a.idealGas / b[2];

  return m;
}

export inline WidomData operator/(const double& a, const WidomData& b)
{
  WidomData m{}; 

  m.total = a / b.total;
  m.excess = a / b.excess;
  m.idealGas = a / b.idealGas;

  return m;
}

export inline WidomData sqrt(const WidomData& a)
{
  WidomData m{};

  m.total = std::sqrt(a.total);
  m.excess = std::sqrt(a.excess);
  m.idealGas = std::sqrt(a.idealGas);

  return m;
}

export inline WidomData log(const WidomData& a)
{
  WidomData m{};

  m.total = std::log(a.total);
  m.excess = std::log(a.excess);
  m.idealGas = std::log(a.idealGas);

  return m;
}

/**
 * \brief Raw Widom terms for the energy-weighted insertion averages.
 *
 * Each Widom insertion contributes its (normalized) Rosenbluth weight W, W * dU with dU the
 * intermolecular energy of the inserted configuration, and W * U_intra with U_intra its
 * intramolecular energy. The CBMC selection probability times W is proportional to the Boltzmann
 * factor of the full energy of the grown molecule, so <W A> / <W> is the Boltzmann average of A over
 * the configurations of a single inserted molecule (Frenkel & Smit, "Understanding Molecular
 * Simulation", Ch. 13). For flexible molecules the ideal-gas reference <U_intra>_IG is sampled the
 * same way from a growth of an isolated molecule (weight W_IG); rigid molecules contribute W_IG = 1
 * and U_intra = 0. All sums are accumulated per block; the ratios are formed afterwards.
 */
export struct WidomEnergyTerms
{
  WidomEnergyTerms() = default;

  WidomEnergyTerms(double weight, double weightedEnergy, double weightedIntraEnergy, double idealGasWeight,
                   double idealGasWeightedIntraEnergy)
      : weight(weight),
        weightedEnergy(weightedEnergy),
        weightedIntraEnergy(weightedIntraEnergy),
        idealGasWeight(idealGasWeight),
        idealGasWeightedIntraEnergy(idealGasWeightedIntraEnergy)
  {
  }

  inline WidomEnergyTerms &operator+=(const WidomEnergyTerms &b)
  {
    weight += b.weight;
    weightedEnergy += b.weightedEnergy;
    weightedIntraEnergy += b.weightedIntraEnergy;
    idealGasWeight += b.idealGasWeight;
    idealGasWeightedIntraEnergy += b.idealGasWeightedIntraEnergy;
    return *this;
  }

  bool operator==(WidomEnergyTerms const &) const = default;

  std::uint64_t versionNumber{2};

  double weight{};                       ///< Rosenbluth weight W of the insertion.
  double weightedEnergy{};               ///< W times the intermolecular insertion energy.
  double weightedIntraEnergy{};          ///< W times the intramolecular energy of the inserted molecule.
  double idealGasWeight{};               ///< Rosenbluth weight W_IG of the isolated (ideal-gas) growth.
  double idealGasWeightedIntraEnergy{};  ///< W_IG times the intramolecular energy of the isolated molecule.

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const WidomEnergyTerms &l);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, WidomEnergyTerms &l);
};

export inline WidomEnergyTerms operator+(const WidomEnergyTerms &a, const WidomEnergyTerms &b)
{
  return WidomEnergyTerms(a.weight + b.weight, a.weightedEnergy + b.weightedEnergy,
                          a.weightedIntraEnergy + b.weightedIntraEnergy, a.idealGasWeight + b.idealGasWeight,
                          a.idealGasWeightedIntraEnergy + b.idealGasWeightedIntraEnergy);
}

export inline WidomEnergyTerms operator-(const WidomEnergyTerms &a, const WidomEnergyTerms &b)
{
  return WidomEnergyTerms(a.weight - b.weight, a.weightedEnergy - b.weightedEnergy,
                          a.weightedIntraEnergy - b.weightedIntraEnergy, a.idealGasWeight - b.idealGasWeight,
                          a.idealGasWeightedIntraEnergy - b.idealGasWeightedIntraEnergy);
}

export inline WidomEnergyTerms operator*(const WidomEnergyTerms &a, const WidomEnergyTerms &b)
{
  return WidomEnergyTerms(a.weight * b.weight, a.weightedEnergy * b.weightedEnergy,
                          a.weightedIntraEnergy * b.weightedIntraEnergy, a.idealGasWeight * b.idealGasWeight,
                          a.idealGasWeightedIntraEnergy * b.idealGasWeightedIntraEnergy);
}

export inline WidomEnergyTerms operator*(const double &a, const WidomEnergyTerms &b)
{
  return WidomEnergyTerms(a * b.weight, a * b.weightedEnergy, a * b.weightedIntraEnergy, a * b.idealGasWeight,
                          a * b.idealGasWeightedIntraEnergy);
}

export inline WidomEnergyTerms operator/(const WidomEnergyTerms &a, const double &b)
{
  return WidomEnergyTerms(a.weight / b, a.weightedEnergy / b, a.weightedIntraEnergy / b, a.idealGasWeight / b,
                          a.idealGasWeightedIntraEnergy / b);
}

export inline WidomEnergyTerms sqrt(const WidomEnergyTerms &a)
{
  return WidomEnergyTerms(std::sqrt(a.weight), std::sqrt(a.weightedEnergy), std::sqrt(a.weightedIntraEnergy),
                          std::sqrt(a.idealGasWeight), std::sqrt(a.idealGasWeightedIntraEnergy));
}

/**
 * \brief Widom-insertion statistics: Rosenbluth weight, chemical potential and fugacity.
 *
 * Three block-averaged channels are sampled: the bare Rosenbluth weight, the raw chemical
 * potential terms (excess Rosenbluth weight and ideal-gas density), and the energy-weighted
 * insertion terms (W and W * dU). The chemical potential, fugacity, insertion energy and enthalpy of
 * adsorption at infinite dilution are non-linear functions of those averages, propagated through
 * BlockAverage::statistics().
 */
export struct PropertyWidom
{
  PropertyWidom() = default;

  PropertyWidom(std::size_t numberOfBlocks)
      : numberOfBlocks(numberOfBlocks),
        rosenbluthWeight(numberOfBlocks),
        chemicalPotentialTerms(numberOfBlocks),
        insertionEnergyTerms(numberOfBlocks)
  {
  }

  std::uint64_t versionNumber{2};

  std::size_t numberOfBlocks;
  BlockAverage<double> rosenbluthWeight;
  BlockAverage<WidomData> chemicalPotentialTerms;
  BlockAverage<WidomEnergyTerms> insertionEnergyTerms;

  std::string writeAveragesRosenbluthWeightStatistics(double temperature, double volume,
                                                      std::optional<double> frameworkMass,
                                                      std::optional<int3> number_of_unit_cells) const;
  std::string writeAveragesChemicalPotentialStatistics(double beta, std::optional<double> imposedChemicalPotential,
                                                       std::optional<double> imposedFugacity) const;
  std::string writeAveragesEnthalpyStatistics(double beta, bool hasFramework, bool rigidComponent) const;

  /// \param insertionEnergy intermolecular energy of the inserted (test) configuration.
  /// \param intraEnergy intramolecular energy of the inserted configuration (zero for rigid molecules).
  /// \param idealGasWeight, idealGasIntraEnergy Rosenbluth weight and intramolecular energy of an isolated
  ///        (ideal-gas) growth of the same molecule (1 and 0 for rigid molecules).
  /// The energies enter only through the energy-weighted channel.
  inline void addWidomSample(std::size_t blockIndex, double RosenbluthValue, double insertionEnergy,
                             double intraEnergy, double idealGasWeight, double idealGasIntraEnergy, std::size_t N,
                             double V, double weight)
  {
    rosenbluthWeight.addSample(blockIndex, RosenbluthValue, weight);
    chemicalPotentialTerms.addSample(blockIndex, WidomData(0.0, RosenbluthValue, static_cast<double>(N) / V), weight);
    insertionEnergyTerms.addSample(
        blockIndex,
        WidomEnergyTerms(RosenbluthValue, RosenbluthValue * insertionEnergy, RosenbluthValue * intraEnergy,
                         idealGasWeight, idealGasWeight * idealGasIntraEnergy),
        weight);
  }

  //====================================================================================================================

  double averagedRosenbluthWeight() const { return rosenbluthWeight.averaged(); }
  double averagedRosenbluthWeight(std::size_t blockIndex) const { return rosenbluthWeight.averaged(blockIndex); }

  std::pair<double, double> result() const { return rosenbluthWeight.average(); }

  //====================================================================================================================

  /// Chemical potential as a non-linear function of the averaged raw terms.
  static WidomData chemicalPotentialTransform(const WidomData &terms, double beta)
  {
    return WidomData(-(1.0 / beta) * std::log(terms.excess) + (1.0 / beta) * std::log(terms.idealGas),
                     -(1.0 / beta) * std::log(terms.excess),
                      (1.0 / beta) * std::log(terms.idealGas));
  }

  WidomData averagedChemicalPotential(double beta) const
  {
    return chemicalPotentialTransform(chemicalPotentialTerms.averaged(), beta);
  }

  WidomData averagedChemicalPotential(std::size_t blockIndex, double beta) const
  {
    return chemicalPotentialTransform(chemicalPotentialTerms.averaged(blockIndex), beta);
  }

  std::pair<WidomData, WidomData> chemicalPotentialResult(double beta) const
  {
    return chemicalPotentialTerms.statistics(
        [beta](const WidomData &terms) { return chemicalPotentialTransform(terms, beta); });
  }

  //====================================================================================================================

  static double fugacityTransform(const WidomData &terms, double beta)
  {
    return std::exp(beta * chemicalPotentialTransform(terms, beta).total) / beta;
  }

  double averagedFugacity(double beta) const { return fugacityTransform(chemicalPotentialTerms.averaged(), beta); }

  double averagedFugacity(std::size_t blockIndex, double beta) const
  {
    return fugacityTransform(chemicalPotentialTerms.averaged(blockIndex), beta);
  }

  std::pair<double, double> fugacityResult(double beta) const
  {
    return chemicalPotentialTerms.statistics([beta](const WidomData &terms) { return fugacityTransform(terms, beta); });
  }

  //====================================================================================================================

  /// Boltzmann-weighted intermolecular energy of a single inserted molecule, <W dU> / <W>.
  static double insertionEnergyTransform(const WidomEnergyTerms &terms)
  {
    return terms.weightedEnergy / terms.weight;
  }

  /// Change of the intramolecular energy on insertion, <W U_intra> / <W> - <W_IG U_intra> / <W_IG>
  /// (zero for rigid molecules).
  static double intraEnergyChangeTransform(const WidomEnergyTerms &terms)
  {
    return terms.weightedIntraEnergy / terms.weight - terms.idealGasWeightedIntraEnergy / terms.idealGasWeight;
  }

  /// Enthalpy of adsorption at infinite dilution, <W dU> / <W> + [<U_intra> - <U_intra>_IG] - k_B T
  /// (adsorption from the ideal gas into an otherwise empty, rigid host; the k_B T accounts for the
  /// p V of the removed gas molecule).
  static double enthalpyTransform(const WidomEnergyTerms &terms, double beta)
  {
    return insertionEnergyTransform(terms) + intraEnergyChangeTransform(terms) - 1.0 / beta;
  }

  std::pair<double, double> insertionEnergyResult() const
  {
    return insertionEnergyTerms.statistics([](const WidomEnergyTerms &terms) { return insertionEnergyTransform(terms); });
  }

  std::pair<double, double> intraEnergyChangeResult() const
  {
    return insertionEnergyTerms.statistics(
        [](const WidomEnergyTerms &terms) { return intraEnergyChangeTransform(terms); });
  }

  std::pair<double, double> enthalpyResult(double beta) const
  {
    return insertionEnergyTerms.statistics(
        [beta](const WidomEnergyTerms &terms) { return enthalpyTransform(terms, beta); });
  }

  double averagedEnthalpy(std::size_t blockIndex, double beta) const
  {
    return enthalpyTransform(insertionEnergyTerms.averaged(blockIndex), beta);
  }

  //====================================================================================================================

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyWidom &w);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyWidom &w);
};
