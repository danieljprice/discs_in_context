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
    },
    # Additional Gould Belt star-forming regions
    'serpens': {
        'galactic': {
            'l_min': 28.0,
            'l_max': 36.0,
            'b_min': 2.0,
            'b_max': 8.0
        },
        'icrs': {
            'ra_min': 18.0,
            'ra_max': 18.8,
            'dec_min': -5.0,
            'dec_max': 5.0
        }
    },
    'aquila': {
        'galactic': {
            'l_min': 20.0,
            'l_max': 40.0,
            'b_min': -5.0,
            'b_max': 10.0
        },
        'icrs': {
            'ra_min': 19.0,
            'ra_max': 20.5,
            'dec_min': -10.0,
            'dec_max': 10.0
        }
    },
    'cepheus': {
        'galactic': {
            'l_min': 95.0,
            'l_max': 115.0,
            'b_min': 10.0,
            'b_max': 25.0
        },
        'icrs': {
            'ra_min': 21.0,
            'ra_max': 23.5,
            'dec_min': 60.0,
            'dec_max': 80.0
        }
    },
    'perseus': {
        'galactic': {
            # Approximate box around the Perseus molecular cloud
            'l_min': 155.0,
            'l_max': 165.0,
            'b_min': -25.0,
            'b_max': -15.0
        },
        'icrs': {
            # Roughly 3–4h in RA, 25–40 deg in Dec
            'ra_min': 3.0,
            'ra_max': 4.0,
            'dec_min': 25.0,
            'dec_max': 40.0
        }
    },
    'scocen': {
        'galactic': {
            # Broad box covering the Sco-Cen OB association (UCL, LCC, Upper Sco)
            'l_min': 285.0,
            'l_max': 360.0,
            'b_min': -15.0,
            'b_max': 30.0
        },
        'icrs': {
            'ra_min': 13.0,
            'ra_max': 17.0,
            'dec_min': -60.0,
            'dec_max': 0.0
        }
    },
    'corona australis': {
        'galactic': {
            # Compact box around the Corona Australis molecular cloud
            'l_min': 0.0,
            'l_max': 15.0,
            'b_min': -24.0,
            'b_max': -10.0
        },
        'icrs': {
            'ra_min': 18.5,
            'ra_max': 20.0,
            'dec_min': -45.0,
            'dec_max': -30.0
        }
    },
    'corona borealis': {
        'galactic': {
            # Rough high-latitude box around Corona Borealis
            'l_min': 20.0,
            'l_max': 50.0,
            'b_min': 20.0,
            'b_max': 50.0
        },
        'icrs': {
            'ra_min': 15.0,
            'ra_max': 17.0,
            'dec_min': 20.0,
            'dec_max': 40.0
        }
    },
    'crux': {
        'galactic': {
            # Small box around Crux on the Galactic plane
            'l_min': 295.0,
            'l_max': 310.0,
            'b_min': -10.0,
            'b_max': 10.0
        },
        'icrs': {
            'ra_min': 11.0,
            'ra_max': 13.0,
            'dec_min': -65.0,
            'dec_max': -50.0
        }
    },
    'allsky': {
        'galactic': {
            'l_min': 0.0,
            'l_max': 360.0,
            'b_min': -90.0,
            'b_max': 90.0
        },
        'icrs': {
            'ra_min': 0.0,
            'ra_max': 24.0,
            'dec_min': -90.0,
            'dec_max': 90.0
        }
    }
}

