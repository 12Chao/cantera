import cantera as ct

gas = ct.Solution('gri30.xml')
gas()
print("reaction expression:", gas.reaction(2), "k_f:",gas.forward_rate_constants[2], "Delta H:", gas.delta_enthalpy[2])
#output reaction expression: H2 + O <=> H + OH k_f: 5195.447373103844 Delta H: 8170274.380677491
#According to Blowers Masel approximation
# Ea = 0 when delta Hrxn < -4Ea0
# Ea = deltaHrxn when delta Hrxn > 4Ea0
# Ea =  (w0+ (deltaH/2))(Vp-2w0+deltaHrxn)^2 / Vp^2 -4w0^2+deltaH^2 ohterwise
# so we input different delta H to see how the forward rate cosntants chagned
delta_h_1 =  gas.delta_enthalpy[2] * 5
delta_h_2 = gas.delta_enthalpy[2] * (-5)
delta_h_3 = gas.delta_enthalpy[2] * 2
gas.add_delta_enthalpy(2, delta_h_1)
print(gas.forward_rate_constants[2])
gas.add_delta_enthalpy(2, delta_h_2)
print(gas.forward_rate_constants[2])
gas.add_delta_enthalpy(2, delta_h_3)
print(gas.forward_rate_constants[2])
#output 52493.11115133174
#       680485879718.8324
#       7143603.71754471
