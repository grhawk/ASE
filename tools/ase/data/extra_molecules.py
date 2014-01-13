"""
Database of molecules outside the G2_1 set

"""

molecule_names = ['Be2','C7NH5','BDA','biphenyl','C60']

data = {
'Be2': {
    'description': "Diatomic Beryllium",
    'name': "Be_2",
    'enthalpy': 155.1,
    'ZPE': 1.0000,
    'thermal correction': 5.0600,
    'symbols': 'BeBe',
    'magmoms': None,
    'positions': [[ 0.  ,  0.  ,  1.0106],
                  [ 0.  ,  0.  , -1.0106]]},
'C7NH5': {
    'description': "Benzonitride",
    'name': "C_7NH_5",
    'symbols': 'C7NH5',
    'magmoms': None,
    'positions': [[ -1.593581, -1.142601, 0.],
                  [ -2.235542,  0.095555, 0.],
                  [ -0.204885, -1.210726, 0.],
                  [  0.549645, -0.025355, 0.],
                  [  1.976332, -0.085321, 0.],
                  [ -0.099258,  1.220706, 0.],
                  [ -1.488628,  1.273345, 0.],
                  [  3.136871, -0.128138, 0.],
                  [ -2.177996, -2.060896, 0.],
                  [ -3.323594,  0.141242, 0.],
                  [  0.301694, -2.173705, 0.],
                  [  0.488716,  2.136782, 0.],
                  [ -1.987765,  2.240495, 0.]]},
'BDA': {
    'description': "1,4-Benzodiamine",
    # aka p-Aminoaniline; p-Benzenediamine; p-Diaminobenzene;
    #     p-Phenylenediamine; Paraphenylen-diamine
    'name': "BDA",
    # PBE-gpaw relaxed
    'symbols': 'C6H4N2H4',
    'magmoms': None,
    'positions': [[ 0.004212,  1.406347,  0.061073],
                  [ 1.193490,  0.687096,  0.029481],
                  [ 1.190824, -0.690400, -0.028344],
                  [ 0.000295, -1.406191, -0.059503],
                  [-1.186974, -0.685668, -0.045413],
                  [-1.185376,  0.690203,  0.009452],
                  [ 2.147124,  1.219997,  0.064477],
                  [ 2.141593, -1.227477, -0.054266],
                  [-2.138408, -1.222814, -0.095050],
                  [-2.137740,  1.226930,  0.023036],
                  [-0.006314,  2.776024,  0.186278],
                  [-0.007340, -2.777839, -0.159936],
                  [ 0.844710, -3.256543,  0.110098],
                  [-0.854965, -3.253324,  0.130125],
                  [ 0.845826,  3.267270, -0.055549],
                  [-0.854666,  3.254654, -0.092676]]},
'biphenyl': {
    'description': "Biphenyl",
    'name': "biphenyl",
    # PBE-gpaw relaxed
    'ionization energy': 8.16,
    'symbols': 'C6H5C6H5',
    'magmoms': None,
    'positions': [[-0.74081, -0.00000, -0.00003],
                  [-1.46261, -1.20370, -0.00993],
                  [-2.85531, -1.20350, -0.00663],
                  [-3.55761, -0.00000, -0.00003],
                  [-2.85531,  1.20350,  0.00667],
                  [-1.46261,  1.20370,  0.00997],
                  [-0.92071, -2.14850,  0.00967],
                  [-3.38981, -2.15110, -0.00083],
                  [-4.64571, -0.00000, -0.00003],
                  [-3.38981,  2.15110,  0.00077],
                  [-0.92071,  2.14850, -0.00963],
                  [ 3.55849, -0.00000, -0.00003],
                  [ 2.85509, -0.86640, -0.83553],
                  [ 1.46289, -0.87000, -0.83153],
                  [ 0.73969, -0.00000, -0.00003],
                  [ 1.46289,  0.87000,  0.83157],
                  [ 2.85509,  0.86640,  0.83547],
                  [ 4.64659, -0.00000, -0.00003],
                  [ 3.39189, -1.53770, -1.50253],
                  [ 0.91869, -1.53310, -1.50263],
                  [ 0.91869,  1.53310,  1.50267],
                  [ 3.39189,  1.53770,  1.50257]]},
'C60': {
    'description': "Buckminsterfullerene, I*h symm.",
    'name': "C_{60}",
    # The Buckyball has two degrees of freedom, the C-C bond, and the C=C bond.
    # This is an LDA-gpaw relaxed structure with bond lengths 1.437 and 1.385.
    # Experimentally, the two bond lengths are 1.45 and 1.40 Angstrom.
    'symbols': 'C60',
    'magmoms': None,
    'positions': [[ 2.2101953,  0.5866631,  2.6669504],
                  [ 3.1076393,  0.1577008,  1.6300286],
                  [ 1.3284430, -0.3158939,  3.2363232],
                  [ 3.0908709, -1.1585005,  1.2014240],
                  [ 3.1879245, -1.4574599, -0.1997005],
                  [ 3.2214623,  1.2230966,  0.6739440],
                  [ 3.3161210,  0.9351586, -0.6765151],
                  [ 3.2984981, -0.4301142, -1.1204138],
                  [-0.4480842,  1.3591484,  3.2081020],
                  [ 0.4672056,  2.2949830,  2.6175264],
                  [-0.0256575,  0.0764219,  3.5086259],
                  [ 1.7727917,  1.9176584,  2.3529691],
                  [ 2.3954623,  2.3095689,  1.1189539],
                  [-0.2610195,  3.0820935,  1.6623117],
                  [ 0.3407726,  3.4592388,  0.4745968],
                  [ 1.6951171,  3.0692446,  0.1976623],
                  [-2.1258394, -0.8458853,  2.6700963],
                  [-2.5620990,  0.4855202,  2.3531715],
                  [-0.8781521, -1.0461985,  3.2367302],
                  [-1.7415096,  1.5679963,  2.6197333],
                  [-1.6262468,  2.6357030,  1.6641811],
                  [-3.2984810,  0.4301871,  1.1204208],
                  [-3.1879469,  1.4573895,  0.1996030],
                  [-2.3360261,  2.5813627,  0.4760912],
                  [-0.5005210, -2.9797771,  1.7940308],
                  [-1.7944338, -2.7729087,  1.2047891],
                  [-0.0514245, -2.1328841,  2.7938830],
                  [-2.5891471, -1.7225828,  1.6329715],
                  [-3.3160705, -0.9350636,  0.6765268],
                  [-1.6951919, -3.0692581, -0.1976564],
                  [-2.3954901, -2.3096853, -1.1189862],
                  [-3.2214182, -1.2231835, -0.6739581],
                  [ 2.1758234, -2.0946263,  1.7922529],
                  [ 1.7118619, -2.9749681,  0.7557198],
                  [ 1.3130656, -1.6829416,  2.7943892],
                  [ 0.3959024, -3.4051395,  0.7557638],
                  [-0.3408219, -3.4591883, -0.4745610],
                  [ 2.3360057, -2.5814499, -0.4761050],
                  [ 1.6263757, -2.6357349, -1.6642309],
                  [ 0.2611352, -3.0821271, -1.6622618],
                  [-2.2100844, -0.5868636, -2.6670300],
                  [-1.7726970, -1.9178969, -2.3530466],
                  [-0.4670723, -2.2950509, -2.6175105],
                  [-1.3283500,  0.3157683, -3.2362375],
                  [-2.1759882,  2.0945383, -1.7923294],
                  [-3.0909663,  1.1583472, -1.2015749],
                  [-3.1076090, -0.1578453, -1.6301627],
                  [-1.3131365,  1.6828292, -2.7943639],
                  [ 0.5003224,  2.9799637, -1.7940203],
                  [-0.3961148,  3.4052817, -0.7557272],
                  [-1.7120629,  2.9749122, -0.7557988],
                  [ 0.0512824,  2.1329478, -2.7937450],
                  [ 2.1258630,  0.8460809, -2.6700534],
                  [ 2.5891853,  1.7227742, -1.6329562],
                  [ 1.7943010,  2.7730684, -1.2048262],
                  [ 0.8781323,  1.0463514, -3.2365313],
                  [ 0.4482452, -1.3591061, -3.2080510],
                  [ 1.7416948, -1.5679557, -2.6197714],
                  [ 2.5621724, -0.4853529, -2.3532026],
                  [ 0.0257904, -0.0763567, -3.5084446]]},
}
