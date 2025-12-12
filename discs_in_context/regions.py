"""
Preset region configurations for Milky Way star-forming regions.
"""

REGIONS = {
    'taurus': {
        'galactic': {
            'l_min': 165.6565,
            'l_max': 180.6565,
            'b_min': -23.4071,
            'b_max': -8.4071
        },
        'icrs': {
            'ra_min': 4.0,
            'ra_max': 5.1667,
            'dec_min': 16.0,
            'dec_max': 32.0
        }
    },
    'lupus': {
        'galactic': {   # e.g. Figure 1 of Tachihara+1996
            'l_min': 334.0,
            'l_max': 344.0,
            'b_min': 8.0,
            'b_max': 18.0
        },
        'icrs': {
            'ra_min': 15.5,
            'ra_max': 16.5,
            'dec_min': -45.0,
            'dec_max': -35.0
        }
    },
    'orion': {
        'galactic': {
            'l_min': 199.0,
            'l_max': 219.0,
            'b_min': -29.5,
            'b_max': -9.5
        },
        'icrs': {
            'ra_min': 5.0,
            'ra_max': 6.0,
            'dec_min': -5.0,
            'dec_max': 5.0
        }
    },
    'ophiuchus': {
        'galactic': {
            'l_min': 348.0,
            'l_max': 358.0,
            'b_min': 12.0,
            'b_max': 22.0
        },
        'icrs': {
            'ra_min': 16.0,
            'ra_max': 17.0,
            'dec_min': -25.0,
            'dec_max': -15.0
        }
    },
    'chamaeleon': {
        'galactic': {
            'l_min': 295.0,
            'l_max': 305.0,
            'b_min': -21.0,
            'b_max': -11.0
        },
        'icrs': {
            'ra_min': 11.0,
            'ra_max': 13.0,
            'dec_min': -80.0,
            'dec_max': -75.0
        }
    }
}

