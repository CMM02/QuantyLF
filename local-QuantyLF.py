from QuantyLF.QuantyLF import QuantyLF


qunatyLF = QuantyLF()
qunatyLF.set_quanty_command('/Users/phr542/Documents/Quanty/Quanty_macOS')

qunatyLF.load_exp_xas('XAS_Exp.dat')
qunatyLF.load_exp_rixs('RIXS_Exp.dat', [638, 639.35, 640.16, 640.75])

qunatyLF.config_edge_jump([[637.7, 0.14, 4], [648.2, 0.006, 8]], [600, 700,], display=False)

# print(qunatyLF.available_cases())
qunatyLF.load_case('D3h_3d', manual=True)
# qunatyLF.load_custom_case('./src/QuantyLF/cases/Td_3d.lua')

# Set up ion and oxidation state
qunatyLF.add_par('ion', 25)
qunatyLF.add_par('oxy', 2)
qunatyLF.add_par('Gamma1', 0.4120250470353196, [0.4, 1])

# # Crystal field contribution in D4h symmetry
qunatyLF.add_par('tenDq', -1.4494920445039643, [-2, 0.1])
qunatyLF.add_par('tenDqF', 0.6815489483432551, [0.01, 1.0])

# # Multiplet contribution
# # spin orbit coupling
qunatyLF.add_par('zeta_2p', 1.0196625781428472, [0.8, 1.02])
qunatyLF.add_par('zeta_3d', 0.8403012992370478, [0.8, 1.02])
qunatyLF.add_par('Xzeta_3d', 1.0, [0.8, 1.02])

# # Slater integrals (Coulomb repulsion/exchange correlation)
qunatyLF.add_par('Fdd', 0.9397729329705585, [0.8, 1.0])
qunatyLF.add_par('XFdd', 0.8137253445941214, [0.8, 1.0])
qunatyLF.add_par('Fpd', 0.8098173584848158, [0.8, 1.0])
qunatyLF.add_par('Gpd', 0.8053014352519605, [0.8, 1.0])


# # Ligand field contribution
# # on-site energies (usually drops out of the equation in crystal field theory)
qunatyLF.add_par('Udd', 6.543685631427877, [2.0, 7.0])
qunatyLF.add_par('Upd_Udd', 4.001467895225598, [0.5, 5.0])

# # Crystal field contribution of ligand site
qunatyLF.add_par('tenDqL', 0.022975132006965073, [0.01, 1.0])

# # Charge transfer contribution
qunatyLF.add_par('Delta', 4.040660314729548, [1.0, 5.0])

# # Hybridization
qunatyLF.add_par('VfScale', 0.9775495515653867, [0.8, 1.0])


# # XAS and RIXS broadening
qunatyLF.add_broadening('XAS', [[-3.7, 0.5], [3, 0.7], [9, 0.7]], gamma=0.5)
qunatyLF.add_broadening('RIXS', [[0, 0.3], [2.7, 0.5], [8, 1]])

qunatyLF.fit('RIXS')