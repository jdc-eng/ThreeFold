"""Inventory of numerical constants, their units, and what they are used for."""


G = 6.6743*10**(-20)        # Grav Constant

# Earth Constants
earthMass = 5.9722 * 10 ** 24 # mass of Earth
earthD = 12742 * 10 ** 3 
AU = 1.495978707 * 10 ** 11

# Moon Constants
moonMass = .07346 * 10 ** 24
moonD = 3474.8 * 10 ** 3

# Earth-Moon System Values
mustar = moonMass/(moonMass + earthMass)
lstar = 3.844 * 10 ** 8
tstar = 382981
EMmu = G*(earthMass+moonMass)


# Sun Constants
ws = 0.925195985520347
sunMass = (1988500 * 10 ** 24)/ (earthMass + moonMass) # Non-dimensionalized mass of of sun in EM system: ms/(me+mm)

# Sail Constants
Beta = 0.01
# pr = 
a_0 = 0.0114  # Nondimensional characteristic acceleration
