from QuantyLF.QuantyLF import QuantyLF


quantyLF = QuantyLF()
quantyLF.set_quanty_command('/Users/phr542/Documents/Quanty/Quanty_macOS', 'Darwin')

quantyLF.load_exp_xas('XAS_Exp.dat')
quantyLF.load_exp_rixs('RIXS_Exp.dat', [638, 639.35, 640.16, 640.75])

quantyLF.config_edge_jump([[637.7, 0.14, 4], [648.2, 0.006, 8]], [600, 700,], display=False)

# print(qunatyLF.available_cases())
quantyLF.load_case('Oh_3d', manual=False)
# qunatyLF.load_custom_case('./src/QuantyLF/cases/Td_3d.lua')

# Set up ion and oxidation state
quantyLF.add_par('ion', 22, from_file=False)
quantyLF.add_par('oxy', 4, from_file=False)
quantyLF.add_par('Gamma1', 0.4120250470353196, [0.4, 1])

# # Crystal field contribution in D4h symmetry
quantyLF.add_par('tenDq', 0.19, [-0.2, 0.2])
quantyLF.add_par('tenDqF', 0.6815489483432551, [0.01, 1.0])


quantyLF.add_par('Ds', 0.999, [0.1, 1])
quantyLF.add_par('Dt', 0.99, [0.1, 1])
quantyLF.add_par('DsF', 0.999, [0.1, 1])
quantyLF.add_par('DtF', 0.999, [0.1, 1])

# # Multiplet contribution
# # spin orbit coupling
quantyLF.add_par('zeta_2p', 1.0196625781428472, [0.8, 1.02])
quantyLF.add_par('zeta_3d', 0.8403012992370478, [0.8, 1.02])
quantyLF.add_par('Xzeta_3d', 1.0, [0.8, 1.02])

# # Slater integrals (Coulomb repulsion/exchange correlation)
quantyLF.add_par('Fdd', 0.9397729329705585, [0.8, 1.0])
quantyLF.add_par('XFdd', 0.8137253445941214, [0.8, 1.0])
quantyLF.add_par('Fpd', 0.8098173584848158, [0.8, 1.0])
quantyLF.add_par('Gpd', 0.8053014352519605, [0.8, 1.0])


# # Ligand field contribution
# # on-site energies (usually drops out of the equation in crystal field theory)
quantyLF.add_par('Udd', 6.543685631427877, [2.0, 7.0])
quantyLF.add_par('Upd_Udd', 4.001467895225598, [0.5, 5.0])

# # Crystal field contribution of ligand site
quantyLF.add_par('tenDqL', 0.022975132006965073, [0.01, 1.0])

# # Charge transfer contribution
quantyLF.add_par('Delta', 4.040660314729548, [1.0, 5.0])

# # Hybridization
quantyLF.add_par('VfScale', 0.9775495515653867, [0.8, 1.0])


# # XAS and RIXS broadening
quantyLF.add_broadening('XAS', [[-3.7, 0.5], [3, 0.7], [9, 0.7]], gamma=0.5)
quantyLF.add_broadening('RIXS', [[0, 0.3], [2.7, 0.5], [8, 1]])

quantyLF.fit('RIXS')