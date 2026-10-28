"""Results retain physical fields and their provenance."""


class SimulationResult:
    def __init__(
        self,
        *,
        fluence,
        flux=None,
        energy_balance=None,
        outgoing_power_by_face=None,
        diagnostics=None,
    ):
        self.fluence = fluence

        self.flux = flux

        self.energy_balance = dict(energy_balance or {})

        self.outgoing_power_by_face = dict(outgoing_power_by_face or {})

        self.diagnostics = dict(diagnostics or {})
