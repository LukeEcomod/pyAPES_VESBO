""" Created in 6.2.2020 by AJ Kieloaho

Gathering parameters needed to establish relationships between moss height [m], bulk density [g m-3] and dry weigth (or dry density) [g dw m-2]

check: 
Rixen and Mulder (2005) Oecologia

"""
nan = None
reference = 'reference'

height = 'height'  # [m]
dry_mass = 'dry_mass'  # [kg m-2]
max_water_content = 'max_water_content'  # [g H2O g-1 dw]

measured_head = 'measured_pressure_head'  # units separately
target_head = 'target_pressure_head'  # units separately

pressure_head = 'pressure_head'  # units separately
water_content = 'water_content'  # units separately
bulk_density = 'bulk_density'  # units separately
hydraulic_conductivity = 'hydraulic_conductivity'  # units separately

bd_unit = 'bulk_density_unit'
ph_unit = 'pressure_head_unit'
wc_unit = 'water_content_unit'
hc_unit = 'hydraulic_conductivity_unit'

water_retention = {
    'S.rubellum and S.fuscum': {
        measured_head: [
            [-0.036968577, -20, -39.96303142, -59.92606285, -80.07393715, -100.0369686, -120],
        ],
        target_head: [
            [0.0, -20.0, -40.0, -60.0, -80.0, -100.0, -120.0],
        ],
        water_content: [
            [96.33401222, 43.54378819, 37.31160896, 35.47861507, 33.09572301, 32.54582485, 30.52953157],
        ],
        bulk_density: [
            0.035,
        ],
        ph_unit: [
            'mb',
        ],
        wc_unit: [
            'percent',
        ],
        bd_unit: [
            'g cm-3'
        ],
        reference: [
            'Cagampan et al. 2007',
        ]
    },
    'Sphagnum fuscum': {
        measured_head: [
            [0.0, -4.006677796, -8.013355593, -12.02003339, -16.02671119, -20.03338898, -24.04006678, -27.9933222, -32.0],
            [0.025020851, -3.999165972, -8.002502085, -16.00917431, -23.99499583],
            [0.025020851, -3.999165972, -7.981651376, -15.9883236, -24.01584654],
        ],
        target_head: [
            [0.0, -4.0, -8.0, -12.0, -16.0, -20.0, -24.0, -28.0, -32.0],
            [0.0, -4.0, -8.0, -16.0, -24.0],
            [0.0, -4.0, -8.0, -16.0, -24.0],
        ],
        water_content: [
            [0.965901639, 0.291803279, 0.241967213, 0.220983607, 0.210491803, 0.2, 0.194754098, 0.189508197, 0.184262295],
            [0.939056025, 0.627724434, 0.475037743, 0.375710311, 0.339259944],
            [0.925513099, 0.711884047, 0.53888216, 0.430848561, 0.399236567],
        ],
        hydraulic_conductivity: [
            [4.073527706, 0.00429963, 0.000861954, 0.000290632, 0.00012412, 7.37964E-05, 4.38762E-05, 1.62609E-05, 9.22167E-06],
            [nan, 2.001124239,  0.560974733, 0.295284527, 0.157836863],
            [nan, 8.018792128, 3.054891137, 0.328752033, 0.174379886]
        ],
        
        bulk_density: [
            0.018, 0.026803419, 0.032752137, 
        ],
        ph_unit: [
            'cm', 'cm', 'cm',
        ],
        wc_unit: [
            'fraction', 'fraction', 'fraction',
        ],
        hc_unit: [
            'cm s-1', 'cm hr-1', 'cm hr-1'
        ],
        bd_unit: [
            'g cm-3', 'g cm-3', 'g cm-3'
        ],
        reference: [
            'Golubev et al. 2018', 'McCarter et al. 2014', 'McCarter et al. 2014',
        ]
    },
    'Sphagnum magellanicum': {
        measured_head: [
            [0.0, -4.976588629, -10.00668896, -17.01672241, -24.02675585, -32.0],
            [0.004170142, -3.978315263, -8.002502085, -15.9883236, -23.99499583],
            [0.004170142, -3.978315263, -7.981651376, -15.9883236, -23.99499583],
            [0.004170142, -4.020016681, -7.981651376, -15.9883236, -23.99499583],
        ],
        target_head: [
            [0.0, -5.0, -10.0, -17.0, -24.0, -32.0],
            [0.0, -4.0, -8.0, -16.0, -24.0],
            [0.0, -4.0, -8.0, -16.0, -24.0],
            [0.0, -4.0, -8.0, -16.0, -24.0],
        ],
        water_content: [
            [0.947540984, 0.45704918, 0.375737705, 0.291803279, 0.26557377, 0.247213115],
            [0.980652963, 0.347191583, 0.298979702, 0.190945296, 0.163201902],
            [0.981620314, 0.287215767, 0.239003079, 0.158055332, 0.137083402],
            [0.980652963, 0.362670827, 0.305750358, 0.20352087, 0.180614236],
        ],
        hydraulic_conductivity: [
            [2.203476912, 0.005194485, 0.00109175, 0.00015721, 1.11409E-05, 4.98824E-06],
            [nan, 7.656957644, 1.47407335, 0.177975149, 0.015811323],
            [nan, 20.91045326, 4.551087243, 0.439910916, 0.175696125],
            [nan, 10.24748514, 1.496888508, 0.349500755, 0.159026539]
        ],
        bulk_density: [
            0.0255, 0.014564103, 0.012923077, 0.012923077,
        ],
        ph_unit: [
            'cm', 'cm', 'cm', 'cm',
        ],
        wc_unit: [
            'fraction', 'fraction', 'fraction', 'fraction',
        ],
        hc_unit:[
            'cm s-1', 'cm hr-1', 'cm hr-1', 'cm hr-1'
        ],
        bd_unit: [
            'g cm-3', 'g cm-3', 'g cm-3', 'g cm-3'
        ],
        reference: [
            'Golubev et al. 2018', 'McCarter et al. 2014', 'McCarter et al. 2014', 'McCarter et al. 2014',
        ]
    },
    'Sphagnum rubellum': {
         measured_head: [
            [0.004170142, -4.020016681, -8.002502085, -16.00917431, -24.01584654],
            [0.004170142, -3.999165972, -8.002502085, -16.00917431, -24.01584654],
            [0.004170142, -4.020016681, -8.002502085, -16.00917431, -24.01584654],
            [0.023419204, 6.908665105, 12.529274, 20.04683841, 26.08899297, 32.76346604, 39.08665105],
            [0.222222222, 5.111111111, 10.22222222, 15.33333333, 20.22222222, 25.11111111, 30.0, 35.11111111, 40.44444444,  60.22222222, 80.22222222, 100.4444444, 120.4444444, 200.2222222],
        ],
        target_head: [
            [0.0, -4.0, -8.0, -16.0, -24.0],
            [0.0, -4.0, -8.0, -16.0, -24.0],
            [0.0, -4.0, -8.0, -16.0, -24.0],
            [0.0, -6.0, -10.0, -14.0, -25.0, -32.0, -40.0],
            [0.0, -5.0, -10.0, -15.0, -20.0, -25.0, -30.0, -35.0, -40.0, -60.0, -80.0, -100.0, -120.0, -200.0],
        ],
        water_content: [
            [0.911003628, 0.491328626, 0.353151407, 0.260595438, 0.217374414],
            [0.940024184, 0.676092028, 0.525340041, 0.462771979, 0.43212653],
            [0.940024184, 0.616117018, 0.503090947, 0.408600275, 0.353771028],
            [0.90730897, 0.318272425, 0.279401993, 0.252491694, 0.240531561, 0.222591362, 0.213621262],
            [86.04166667, 73.95833333, 54.16666667, 38.95833333, 33.125, 29.375, 26.04166667, 25.20833333, 23.95833333, 22.08333333, 19.79166667, 18.125, 17.08333333, 15.0],
        ],
        hydraulic_conductivity: [
            [nan, 1.188134936, 0.449181639, 0.178028235, 0.156652749],
            [nan, 2.211264233, 0.733799569, 0.204380218, 0.147330814],
            [nan, 3.449949789, 1.463279585, 0.39222979, 0.189758776],
            [],
            []
        ],
        bulk_density: [
            0.022, 0.031111111, 0.023316239, 0.033162393, 0.04, 24.0
        ],
        ph_unit: [
            'cm', 'cm', 'cm', 'cm', 'mb',
        ],
        wc_unit: [
            'fraction', 'fraction', 'fraction', 'fraction', 'percent'
        ],
        hc_unit: [
            'cm hr-1', 'cm hr-1', 'cm hr-1', nan, nan
        ],
        bd_unit: [
            'g cm-3','g cm-3','g cm-3', 'g cm-3', 'g cm-3', 'kg m-3'
        ],
        reference: [
            'McCarter et al. 2014', 'McCarter et al. 2014', 'McCarter et al. 2014', 'Price et al. 2008', 'Waddington et al. 2011'
        ]
    },
    'Sphagnum capillifolium': {
        measured_head: [
            [9.652650823, 19.72577697, 29.79890311, 39.87202925, 49.94515539, 74.95429616, 99.96343693, 149.9817185, 200.0],
        ],
        target_head: [
            [-10.0, -20.0, -30.0, -40.0, -50.0, -75.0, -100.0, -150.0, -200.0]
        ],
        water_content: [
            [0.337435897, 0.275726496, 0.24974359, 0.233504274, 0.223760684, 0.207521368, 0.197777778, 0.191282051, 0.178290598],
        ],
        bulk_density: [
            14.0,
        ],
        ph_unit: [
            'mb'
        ],
        bd_unit: [
            'kg m-3'
        ],
        wc_unit: [
            'fraction'
        ],
        reference: [
            'Moore et al. 2018',
        ]
    },
    'Sphagnum papillosum': {
        measured_head: [
            [9.652650823, 19.72577697, 29.79890311, 39.52468007, 49.597806215722, 74.95429616, 99.61608775, 149.9817185, 200.0],
        ],
        target_head: [
            [-10.0, -20.0, -30.0, -40.0, -50.0, -75.0, -100.0, -150.0, -200.0]
        ],
        water_content: [
            [0.238376068, 0.150683761, 0.131196581, 0.123076923, 0.116581197, 0.11008547, 0.105213675, 0.098717949, 0.088974359],
        ],
        bulk_density: [
            14.0,
        ],
        ph_unit: [
            'mb'
        ],
        wc_unit: [
            'fraction'
        ],
        bd_unit: [
            'kg m-3'
        ],
        reference: [
            'Moore et al. 2018',
        ]
    },
    'Sphagnum angustifolium': {
        measured_head: [
            [9.652650823, 19.72577697, 29.79890311, 39.87202925, 49.59780622, 74.95429616, 99.96343693, 149.9817185, 200.0],
        ],
        target_head: [
            [-10.0, -20.0, -30.0, -40.0, -50.0, -75.0, -100.0, -150.0, -200.0]
        ],
        water_content: [
            [0.121452991, 0.085726496, 0.069487179, 0.066239316, 0.066239316, 0.064615385, 0.064615385, 0.064615385, 0.061367521],
        ],
        bulk_density: [
            14.0,
        ],
        ph_unit: [
            'mb'
        ],
        wc_unit: [
            'fraction'
        ],
        bd_unit: [
            'kg m-3'
        ],
        reference: [
            'Moore et al. 2018',
        ]
    },
    'Polytrichum piliferum': {
        measured_head: [
            [9.722222222, 19.79166667, 29.86111111, 39.93055556, 49.65277778, 74.65277778, 99.65277778, 149.6527778, 199.6527778],
            [-100.0344944, -50.11560621, -15.00231212, -5.120417932, -1.109069776]
        ],
        target_head: [
            [-10.0, -20.0, -30.0, -40.0, -50.0, -75.0, -100.0, -150.0, -200.0],
            [-98.0, -51.0, -15.0, -5.0, -1.0]
        ],
        water_content: [
            [0.435294118, 0.341176471, 0.302231237, 0.282758621, 0.263286004, 0.25030426, 0.240567951, 0.227586207, 0.212981744],
            [0.070018622, 0.111731844, 0.274115456, 0.426070764, 0.472253259],
        ],
        bulk_density: [
            26.0, 0.020372168
        ],
        ph_unit: [
            'mb', 'cm'
        ],
        wc_unit: [
            'fraction', 'fraction'
        ],
        bd_unit: [
            'kg m-3', 'g cm-3'
        ],
        reference: [
            'Moore et al. 2018', 'Voortman et al. 2014',
        ]
    },
    'Polytrichum strictum': {
        measured_head: [
            [9.722222222, 19.79166667, 29.51388889, 39.93055556, 49.65277778, 74.65277778, 99.65277778, 149.6527778, 199.6527778]
        ],
        target_head: [
            [-10.0, -20.0, -30.0, -40.0, -50.0, -75.0, -100.0, -150.0, -200.0]
        ],
        water_content: [
            [0.393103448, 0.346044625, 0.323326572, 0.30872211, 0.295740365, 0.277890467, 0.263286004, 0.245436105, 0.225963489],
        ],
        bulk_density: [
            26.0,
        ],
        ph_unit: [
            'mb'
        ],
        wc_unit: [
            'fraction'
        ],
        bd_unit: [
            'kg m-3'
        ],
        reference: [
            'Moore et al. 2018',
        ]
    },
    'Polytrichum commune': {
        measured_head: [
            [9.722222222, 19.79166667, 29.86111111, 39.93055556, 50.0, 74.65277778, 99.65277778, 150.0, 199.6527778],
        ],
        target_head: [
            [-10.0, -20.0, -30.0, -40.0, -50.0, -75.0, -100.0, -150.0, -200.0]
        ],
        water_content: [
            [0.352535497, 0.305476673, 0.279513185, 0.261663286, 0.247058824, 0.232454361, 0.219472617, 0.20811359, 0.193509128],
        ],
        bulk_density: [
            26.0,
        ],
        ph_unit: [
            'mb'
        ],
        wc_unit: [
            'fraction'
        ],
        bd_unit: [
            'kg m-3'
        ],
        reference: [
            'Moore et al. 2018',
        ]
    },
    'Hypnum cupressiforme': {
        measured_head: [
            [-99.73244147, -50.23411371, -15.18394649, -5.284280936, -1.27090301]
        ],
        target_head: [
            [-98.0, -51.0, -15.0, -5.0, -1.0]
        ],
        water_content: [
            [0.026815642, 0.041713222, 0.09981378, 0.17132216, 0.201117318]
        ],
        bulk_density: [
            0.011308882,
        ],
        ph_unit: [
            'cm'
        ],
        wc_unit: [
            'fraction'
        ],
        bd_unit: [
            'g cm-3'
        ],
        reference: [
            'Voortman et al. 2014'
        ]
    },
    'Campolypus introflexus': {
        measured_head: [
            [-100.3199279, -50.23944526, -15.02634774, -5.119597462, -1.091710162]
        ],
        target_head: [
            [-98.0, -51.0, -15.0, -5.0, -1.0]
        ],
        water_content: [
            [0.118959108, 0.179925651, 0.397026022, 0.573977695, 0.611152416]
        ],
        bulk_density: [
            0.026673139
        ],
        ph_unit: [
            'cm'
        ],
        wc_unit: [
            'fraction'
        ],
        bd_unit: [
            'g cm-3'
        ],
        reference: [
            'Voortman et al. 2014'
        ]
    },
    'Syntrichia ruralis': {
        measured_head: [
            [-99.73154362, -50.06711409, -14.89932886, -4.966442953, -0.939597315]
        ],
        target_head: [
            [-98.0, -51.0, -15.0, -5.0, -1.0]
        ],
        water_content: [
            [0.028252788, 0.068401487, 0.27063197, 0.447583643, 0.495167286],
        ],
        bulk_density: [
            0.0098274
        ],
        ph_unit: [
            'cm'
        ],
        wc_unit: [
            'fraction'
        ],
        bd_unit: [
            'g cm-3'
        ],
        reference: [
            'Voortman et al. 2014'
        ]
    }
}

library = {
    'Aulacomniun palustre': {
        height: [0.044, 0.045, 0.034],
        dry_mass: [1.03, 0.8, 0.510],
        bulk_density:[23.5, 17.5, 14.86],
        max_water_content: [nan, 20.25, nan],
        reference: ['Soudziloskaia et al. 2013', 'Elumeeva et al. 2011', 'Michel et al. 2012']
    },
    'A.palustre and C.stygium': {
        height: [0.027],
        dry_mass: [0.377],
        bulk_density: [14.18],
        reference: ['Michel et al. 2012']
    },
    'A.palustre and D.scoparium': {
        height: [0.028],
        dry_mass: [0.678],
        bulk_density: [23.95],
        reference: ['Michel et al. 2012']
    },
    'A.palustre and H.splendens': {
        height: [0.034],
        dry_mass: [0.484],
        bulk_density: [14.41],
        reference: ['Michel et al. 2012']
    },
    'Aulacomnium turgidum': {
        height: [0.042, 0.034,],
        dry_mass: [1.74, 1.1],
        bulk_density: [41.0, 37.6],
        max_water_content: [nan, 14.52],
        reference: ['Soudziloskaia et al. 2013', 'Elumeeva et al., 2011']
    },
    'Campylus clavatus': {
        height: [0.0131],
        dry_mass: [0.054],
        bulk_density: [4.12],
        max_water_content: [8.2],
        reference: ['Michel et al. 2013']
    },
    'Cinclidium stygium': {
        height: [0.054],
        dry_mass: [0.341],
        bulk_density: [6.35],
        reference: ['Michel et al. 2012']
    },
    'Dicranum dummondii': {
        height: [0.028],
        dry_mass: [2.58],
        bulk_density: [85.70],
        max_water_content: [6.03],
        reference: ['Lett et al. 2017']
    },
    'Dicranum scoparium': {
        height: [0.036, 0.034, 0.022, 0.025],
        dry_mass: [1.59, 0.8, 0.642, 0.61],
        bulk_density: [44.0, 22.4, 29.72, 37.1],
        max_water_content: [nan, 12.42, nan, nan],
        reference: ['Soudziloskaia et al. 2013', 'Elumeeva et al. 2011', 'Michel et al. 2012', 'Soudziloskaia et al. 2011']
    },
    'Drepanocladus cossonii': {
        height: [0.044],
        dry_mass: [0.98],
        bulk_density: [22.0],
        max_water_content: [16.73],
        reference: ['Soudziloskaia et al. 2013']
    },
    'Hylocomium splendens': {
        height: [0.060, 0.066, 0.046, 0.078, 0.073, 0.045],
        dry_mass: [0.86, 0.9, 0.606, 0.651, 1.45, 1.00],
        bulk_density: [14.3, 13.2, 13.1, 7.30, 19.7, 23.41],
        max_water_content: [nan, 9.44, nan, 9.2, 9.93, nan],
        reference: ['Soudziloskaia et al. 2013', 'Elumeeva et al., 2011', 'Michel et al. 2012', 'Stuiver et al. 2014', 'Soudziloskaia et al. 2011', 'Lett et al. 2013']
    },
    'H.splendens and L.lycopodioides': {
        height: [0.021],
        dry_mass: [0.550],
        bulk_density: [25.84],
        reference: ['Michel et al. 2012']
    },
    'H.splendens and P.schreberi': {
        height: [0.047],
        dry_mass: [0.522],
        bulk_density: [11.17],
        reference: ['Michel et al. 2012']
    },
    'H.splendens and P.commune': {
        height: [0.057],
        dry_mass: [0.668],
        bulk_density: [11.81],
        reference: ['Michel et al. 2012']
    },
    'Hypnum cupressiforme': {
        height: [0.061],
        dry_mass: [1.012],
        bulk_density: [16.70],
        max_water_content: [14.7],
        reference: ['Michel et al. 2013']
    },
    'Leptotheca gaudihaudii': {
        height: [0.0414],
        dry_mass: [0.152],
        bulk_density: [3.67],
        max_water_content: [11.4],
        reference: ['Michel et al. 2013']
    },
    'Limprichtia cossonii': {
        height: [0.051],
        dry_mass: [1.3],
        bulk_density: [28.5],
        max_water_content: [17.33],
        reference: ['Elumeeva et al. 2011']
    },
    'Lohozia floearkii': {
        height: [0.032],
        dry_mass: [1.69],
        bulk_density: [50.79],
        max_water_content: [8.67],
        reference: ['Lett et al. 2017']
    },
    'Oncophorus wahlenbergii': {
        height: [0.035, 0.021,],
        dry_mass: [2.00, 1.3,],
        bulk_density: [57.0, 61.1,],
        max_water_content: [nan, 14.2],
        reference: ['Soudziloskaia et al. 2013', 'Elumeeva et al., 2011']
    },
    'Paludella squarrosa': {
        height: [0.06, 0.03],
        dry_mass: [1.18, 0.9],
        bulk_density: [19.7, 29.4],
        max_water_content: [nan, 22.45],
        reference: ['Soudziloskaia et al. 2013', 'Elumeeva et al., 2011']
    },
    'Philonotis coespitosa': {
        height: [0.036],
        dry_mass: [1.7],
        bulk_density: [47.7],
        max_water_content: [12.94],
        reference: ['Elumeeva et al., 2011']
    },
    'Pleurozium schreberi': {
        height: [0.095, 0.074, 0.069, 0.042, 0.075, 0.054],
        dry_mass: [1.62, 1.1, 0.651, 0.581, 1.40, 1.49],
        bulk_density: [17.1, 15.2, 9.47, 14.12, 18.8, 28.74],
        max_water_content: [nan, 10.43, 11.8, nan, 7.91, nan],
        reference: ['Soudziloskaia et al. 2013', 'Elumeeva et al., 2011', 'Stuiver et al. 2014', 'Michel et al. 2012', 'Soudziloskaia et al. 2011', 'Lett et al. 2017']
    },
    'Plogiomnium ellipticum': {
        height: [0.048],
        dry_mass: [0.6],
        bulk_density: [12.4],
        max_water_content: [21.38],
        reference: ['Elumeeva et al. 2011']
    },
    'Polytrichastrum sexangulare': {
        height: [0.037, 0.038,],
        dry_mass: [2.04, 2.4,],
        bulk_density: [54.8, 64.8,],
        max_water_content: [nan, 6.56],
        reference: ['Soudziloskaia et al. 2013', 'Elumeeva et al., 2011']
    },
    'Polytrichum strictum': {
        height: [0.048, 0.051],
        dry_mass: [2.56, 3.0],
        bulk_density: [53.8, 58.3],
        max_water_content: [nan, 5.57],
        reference: ['Soudziloskaia et al. 2013', 'Elumeeva et al., 2011']
    },
    'Polytrichum juniperinum': {
        height: [0.047, 0.0539],
        dry_mass: [1.6, 0.546],
        bulk_density: [33.3, 10.13],
        max_water_content: [4.91, 6.6],
        reference: ['Elumeeva et al., 2011', 'Michel et al., 2013']
    },
    'Polytrichum commune': {
        height: [0.072, 0.074, 0.048, 0.073],
        dry_mass: [0.661, 0.795, 1.28, 2.68],
        bulk_density: [9.15, 10.77, 26.9, 41.35],
        max_water_content: [nan, 4.68, nan, 6.12],
        reference: ['Stuiver et al. 2014', 'Michel et al. 2012', 'Soudziloskaia et al. 2011', 'Lett et al. 2017']
    },
    'P.commune and P.schreberi': {
        height: [0.044],
        dry_mass: [0.566],
        bulk_density: [12.77],
        reference: ['Michel et al. 2012']
    },
    'Ptilidum ciliate': {
        height: [0.023, 0.023, 0.036],
        dry_mass: [0.6, 0.94, 1.72],
        bulk_density: [24.8, 40.5, 45.83],
        max_water_content: [12.29, nan, 9.18],
        reference: ['Elumeeva et al. 2011', 'Soudziloskaia et al. 2011']
    },
    'Racomitrium fasciculare': {
        height: [0.036],
        dry_mass: [1.4],
        bulk_density: [39.7],
        max_water_content: [11.61],
        reference: ['Elumeeva et al. 2011']
    },
    'Racomitrium lanuginosum': {
        height: [0.065],
        dry_mass: [4.95],
        bulk_density: [76.6],
        max_water_content: [6.8],
        reference: ['Soudziloskaia et al. 2013']
    },
    'Racomitrium pruinosum': {
        height: [0.0776],
        dry_mass: [0.224],
        bulk_density: [2.89],
        max_water_content: [10.9],
        reference: ['Michel et al. 2013']
    },
    'Racomitrium tonuginosum': {
        height: [0.044],
        dry_mass: [3.3],
        bulk_density: [75.3],
        max_water_content: [6.80],
        reference: ['Elumeeva et al. 2011']
    },
    'Rhytidium rugosum':{
        height: [0.04],
        dry_mass: [1.71],
        bulk_density: [42.7],
        reference: ['Soudziloskaia et al. 2013']
    },
    'Sphagnum capillifolium': {
        height: [0.058],
        dry_mass: [2.43],
        bulk_density: [46.66],
        max_water_content: [15.63],
        reference: ['Lett et al. 2017']
    },
    'Sphagnum fuscum': {
        height: [0.044, 0.046, 0.057],
        dry_mass: [1.58, 1.6, 1.40],
        bulk_density: [35.1, 34.9, 27.09],
        max_water_content: [nan, 23.73, 17.79],
        reference: ['Soudziloskaia et al. 2013', 'Elumeeva et al. 2011']
    },
    'Sphagnum girgensohnii': {
        height: [0.064],
        dry_mass: [0.592],
        bulk_density: [9.28],
        max_water_content: [17.18],
        reference: ['Stuiver et al. 2014']
    },
    'Sphagnum lindbergii': {
        height: [0.065],
        dry_mass: [1.8],
        bulk_density: [27.1],
        max_water_content: [18.12],
        reference: ['Elumeeva et al. 2011']
    },
    'Sphagnum riporium': {
        height: [0.076],
        dry_mass: [1.2],
        bulk_density: [16.1],
        max_water_content: [16.86],
        reference: ['Elumeeva et al. 2011']
    },
    'Sphagnum russowii': {
        height: [0.066],
        dry_mass: [1.4],
        bulk_density: [20.6],
        max_water_content: [24.35],
        reference: ['Elumeeva et al., 2011']
    },
    'Tomentypnum nitens': {
        height: [0.051, 0.051, 0.04],
        dry_mass: [1.41, 1.3, 0.474],
        bulk_density: [27.8, 25.0, 11.76],
        max_water_content: [nan, 16.73, nan],
        reference: ['Soudziloskaia et al. 2013', 'Elumeeva et al., 2011', 'Michel et al. 2012']
    },
    'T.nitens and C.stygium' : {
        height: [0.044],
        dry_mass: [0.372],
        bulk_density: [8.45],
        reference: ['Michel et al. 2012']
    },
    'Warnstorfia pseudostaminea': {
        height: [0.055],
        dry_mass: [1.8],
        bulk_density: [33.0],
        max_water_content: [17.69],
        reference: ['Elumeeva et al., 2011']
    }
}