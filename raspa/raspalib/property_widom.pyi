from raspalib.widom_data import *

class PropertyWidom:
    """
    Configure and access Widom-insertion sampling results in RASPA.
    """

    def __init__(self) -> None:
        ...
        """
        Initialize a :class:`PropertyWidom`.
        """

    def chemical_potential_result(self, temperature: float) -> tuple[WidomData, WidomData]:
        """Return Widom chemical-potential aggregates at a given temperature.

        Args:
            temperature: Temperature used to convert Widom statistics.

        Returns:
            Pair of :class:`WidomData` aggregates for chemical potential.
        """
        ...

    def fugacity_result(self, temperature: float) -> float:
        """Return fugacity derived from Widom statistics at a temperature.

        Args:
            temperature: Temperature used to convert Widom statistics.

        Returns:
            Fugacity in Pascal.
        """
        ...



    def insertion_energy_result(self) -> tuple[float, float]:
        """Return the Boltzmann-weighted Widom insertion energy, <W dU>/<W>.

        dU is the intermolecular energy of the inserted test molecule (intramolecular terms
        excluded).

        Returns:
            Pair of mean and error in Kelvin.
        """
        ...

    def enthalpy_of_adsorption_result(self, temperature: float) -> tuple[float, float]:
        """Return the enthalpy of adsorption at infinite dilution, <W dU>/<W> - k_B T.

        Valid for Widom insertions into an otherwise empty host.

        Args:
            temperature: Temperature in Kelvin.

        Returns:
            Pair of mean and error in Kelvin.
        """
        ...
