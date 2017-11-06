"""
All best parameters from parameter sweeps
"""

open2_mm = {"Enzymes": {
    "cds": {"k": 248.299, "kinetics": "michaelis_menten", "v": 4.7169},
    "dagk": {"k": 13.0596, "kinetics": "michaelis_menten", "v": 8.5278},
    "ip3ptase": {"k": 935.7779, "kinetics": "michaelis_menten", "v": 7.2798},
    "laza": {"k": 19.1597, "kinetics": "michaelis_menten", "v": 9.9025},
    "p4tase": {"k": 317.6421, "kinetics": "michaelis_menten", "v": 18.9616},
    "p5tase": {"k": 676.8825, "kinetics": "michaelis_menten", "v": 4.3239},
    "patp": {"k": 629.2997, "kinetics": "michaelis_menten", "v": 0.0118},
    "pi4k": {"k": 80.8653, "kinetics": "michaelis_menten", "v": 0.5888},
    "pip5k": {"k": 29.3552, "kinetics": "michaelis_menten", "v": 1.9867},
    "pis": {"k": 3.3836, "kinetics": "michaelis_menten", "v": 9.7055},
    "pitp": {"k": 0.2273, "kinetics": "michaelis_menten", "v": 9.9013},
    "plc": {"k": 21.7648, "kinetics": "michaelis_menten", "v": 1.5667},
    "sink": {"k": 18.1918, "kinetics": "michaelis_menten", "v": 6.6503},
    "source": {"k": 0.3656, "kinetics": "michaelis_menten", "v": 17.1201}},
    "Mutants": {"0.1": {"LAZA": {"DAG": 0.9957, "PA": 2.4998},
                        "RDGA": {"DAG": 1.0002, "PA": 0.9461}}}}

open2_mass_action = {
    "Enzymes": {"cds": {"k": 0.8679, "kinetics": "mass_action", "v": 1},
                "dagk": {"k": 125.7593, "kinetics": "mass_action", "v": 1},
                "ip3ptase": {"k": 10, "kinetics": "mass_action", "v": 1},
                "laza": {"k": 36.9801, "kinetics": "mass_action", "v": 1},
                "p4tase": {"k": 10, "kinetics": "mass_action", "v": 1},
                "p5tase": {"k": 10, "kinetics": "mass_action", "v": 1},
                "patp": {"k": 0.0286, "kinetics": "mass_action", "v": 1},
                "pi4k": {"k": 0.1227, "kinetics": "mass_action", "v": 1},
                "pip5k": {"k": 2.4389, "kinetics": "mass_action", "v": 1},
                "pis": {"k": 121.9473, "kinetics": "mass_action", "v": 1},
                "pitp": {"k": 20.0771, "kinetics": "mass_action", "v": 1},
                "plc": {"k": 2.4389, "kinetics": "mass_action", "v": 1},
                "sink": {"k": 15.1462, "kinetics": "mass_action", "v": 1},
                "source": {"k": 45.8352, "kinetics": "mass_action", "v": 1}},
    "Mutants": {"0.1": {"LAZA": {"DAG": 0.9461, "PA": 2.3611},
                        "RDGA": {"DAG": 1.0058, "PA": 0.8542}}}}
