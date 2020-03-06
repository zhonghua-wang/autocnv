DEFAULT_SCORE = {
    '1A': 0, '1B': -0.6,
    '2A': 1, '2B': 0,
    '2C': -1, '2C-1': 0.9, '2C-2': 0,
    '2D': -1, '2D-1': 0, '2D-2': 0.9, '2D-3': 0.3, '2D-4': 0.9,
    '2E': 0, '2F': -1, '2G': 0, '2H': 0, '2I': 0, '2J': 0, '2K': 0.45,
    '3A': 0, '3B': 0.45, '3C': 0.9,
    '4O': -1,
    'PVS1': 0.9, 'PVS1_S': 0.45, 'PVS1_M': 0.3, 'PVS1_P': 0.15, 'PVS1_U': 0
}

try:
    from local_settings import *
except ImportError:
    pass

try:
    from acit.local_settings import *
except ImportError:
    pass
