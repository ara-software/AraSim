# AraSim/data/antennas/

This folder contains data on the various antennas used in AraSim.  There are differences between ideal and realized gain that must be considered.

# Realized Gain

This is the default gain used in AraSim and is the most useful.  This has an implicit consideration for the impedance mismatch in the circuit.  When using any of these gain values, whether for Tx or Rx, you should keep the impedance of your antenna to the default value of 50 Ohms.  This is dictated by the settings IMPEDANCE_TX, IMPEDANCE_RX_VPOL, IMPEDANCE_RX_VPOL_TOP, and IMPEDANCE_RX_HPOL; which should be set to the default value of 0.

# Ideal Gain

The ideal gain represents the gain of the antenna when there is no impedance mismatch, which means that your formula for effective height MUST include the impedance of the antenna as a function of frequency.  AraSim does not natively use ideal gain, but we felt the need to include this note for any future implementations.

# Impedance

These are the frequency-dependent impedances of various antennas in AraSim that are used in conjunction with the ideal gain.  At present, these impedances are not used with the realized gain implementation in AraSim.  The data is provided for the sake of completeness.