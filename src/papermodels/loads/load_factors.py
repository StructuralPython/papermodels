def nbcc_2020_combos():
    nbcc_2020_combinations = {
        "LC1": {"D": 1.4},
        "LC2a": {"D": 1.25, "L": 1.5},
        "LC2b": {"D": 1.25, "L": 1.5, "S": 1.0},
        "LC2c": {"D": 1.25, "L": 1.5, "W": 0.4},
        "LC2d": {"D": 0.9, "L": 1.5},
        "LC2e": {"D": 0.9, "L": 1.5, "S": 1.0},
        "LC2f": {"D": 0.9, "L": 1.5, "W": 0.4},
        "LC3a": {"D": 1.25, "S": 1.5},
        "LC3b": {"D": 1.25, "S": 1.5, "L": 1.0},
        "LC3c": {"D": 1.25, "S": 1.5, "W": 0.4},
        "LC3d": {"D": 0.9, "S": 1.5},
        "LC3e": {"D": 0.9, "S": 1.5, "L": 1.0},
        "LC3f": {"D": 0.9, "S": 1.5, "W": 0.4},
        "LC4a": {"D": 1.25, "W": 1.4},
        "LC4b": {"D": 1.25, "W": 1.4, "L": 0.5},
        "LC4c": {"D": 1.25, "W": 1.4, "S": 0.5},
        "LC4d": {"D": 0.9, "W": 1.4},
        "LC4e": {"D": 0.9, "W": 1.4, "L": 0.5},
        "LC4f": {"D": 0.9, "W": 1.4, "S": 0.5},
        "LC5a": {"D": 1.0, "E": 1.0},
        "LC5b": {"D": 1.0, "E": 1.0, "L": 0.5, "S": 0.25},
    }
    return nbcc_2020_combinations


def nbcc_2015_combos():
    nbcc_2015_combinations = {
        "LC1": {"D": 1.4},
        "LC2a": {"D": 1.25, "L": 1.5},
        "LC2b": {"D": 1.25, "L": 1.5, "S": 1.0},
        "LC2c": {"D": 1.25, "L": 1.5, "W": 0.4},
        "LC2d": {"D": 0.9, "L": 1.5},
        "LC2e": {"D": 0.9, "L": 1.5, "S": 1.0},
        "LC2f": {"D": 0.9, "L": 1.5, "W": 0.4},
        "LC3a": {"D": 1.25, "S": 1.5},
        "LC3b": {"D": 1.25, "S": 1.5, "L": 1.0},
        "LC3c": {"D": 1.25, "S": 1.5, "W": 0.4},
        "LC3d": {"D": 0.9, "S": 1.5},
        "LC3e": {"D": 0.9, "S": 1.5, "L": 1.0},
        "LC3f": {"D": 0.9, "S": 1.5, "W": 0.4},
        "LC4a": {"D": 1.25, "W": 1.4},
        "LC4b": {"D": 1.25, "W": 1.4, "L": 0.5},
        "LC4c": {"D": 1.25, "W": 1.4, "S": 0.5},
        "LC4d": {"D": 0.9, "W": 1.4},
        "LC4e": {"D": 0.9, "W": 1.4, "L": 0.5},
        "LC4f": {"D": 0.9, "W": 1.4, "S": 0.5},
        "LC5a": {"D": 1.0, "E": 1.0},
        "LC5b": {"D": 1.0, "E": 1.0, "L": 0.5, "S": 0.25},
    }
    return nbcc_2015_combinations


def nbcc_2010_combos():
    nbcc_2010_combinations = {
        "LC1": {"D": 1.4},
        "LC2a": {"D": 1.25, "L": 1.5},
        "LC2b": {"D": 1.25, "L": 1.5, "S": 0.5},
        "LC2c": {"D": 1.25, "L": 1.5, "W": 0.4},
        "LC2d": {"D": 0.9, "L": 1.5},
        "LC2e": {"D": 0.9, "L": 1.5, "S": 0.5},
        "LC2f": {"D": 0.9, "L": 1.5, "W": 0.4},
        "LC3a": {"D": 1.25, "S": 1.5},
        "LC3b": {"D": 1.25, "S": 1.5, "L": 0.5},
        "LC3c": {"D": 1.25, "S": 1.5, "W": 0.4},
        "LC3d": {"D": 0.9, "S": 1.5},
        "LC3e": {"D": 0.9, "S": 1.5, "L": 0.5},
        "LC3f": {"D": 0.9, "S": 1.5, "S": 0.5},
        "LC4a": {"D": 1.25, "W": 1.4},
        "LC4b": {"D": 1.25, "W": 1.4, "L": 0.5},
        "LC4c": {"D": 1.25, "W": 1.4, "S": 0.5},
        "LC4d": {"D": 0.9, "W": 1.4},
        "LC4e": {"D": 0.9, "W": 1.4, "L": 0.5},
        "LC4f": {"D": 0.9, "W": 1.4, "S": 0.5},
        "LC5a": {"D": 1.0, "E": 1.0},
        "LC5b": {"D": 1.0, "E": 1.0, "L": 0.5, "S": 0.25},
    }
    return nbcc_2010_combinations


def factor_load(
    D_load: float = 0.0,
    D: float = 0.0,
    L_load: float = 0.0,
    L: float = 0.0,
    S_load: float = 0.0,
    S: float = 0.0,
    W_load: float = 0.0,
    W: float = 0.0,
    E_load: float = 0.0,
    E: float = 0.0,
) -> float:
    """
    Returns the factored load given the load components (e.g. 'D_load', 'L_load', etc.)
    and the load factors (e.g. 'D', 'L', etc.) provided.
    """
    factored_load = D_load * D + L_load * L + S_load * S + W_load * W + E_load * E
    return factored_load


def max_factored_load(
    loads: dict,
    load_combos: dict,
) -> float:
    """
    Returns the max factored load from the given 'loads' and 'load_combos'.
    'loads' is a dict with keys from {"D_load", "L_load", "S_load", "W_load", "E_load"}
        and float values representing load magnitudes. +ve/-ve is dependent on the designer's
        sign convention.
    'load_combos' is a dict of str keys (representing the "names" of each load combination)
        with dict values representing load components and their factors, e.g.
        {"D": 1.25, "L": 1.5, "S": 1.0}
    """
    factored_loads = []
    for load_combo in load_combos.values():
        factored_load = factor_load(**loads, **load_combo)
        factored_loads.append(factored_load)
    return max(factored_loads)


def min_factored_load(
    loads: dict,
    load_combos: dict,
) -> float:
    """
    Returns the max factored load from the given 'loads' and 'load_combos'.
    'loads' is a dict with keys from {"D_load", "L_load", "S_load", "W_load", "E_load"}
        and float values representing load magnitudes. +ve/-ve is dependent on the designer's
        sign convention.
    'load_combos' is a dict of str keys (representing the "names" of each load combination)
        with dict values representing load components and their factors, e.g.
        {"D": 1.25, "L": 1.5, "S": 1.0}
    """
    factored_loads = []
    for load_combo in load_combos.values():
        factored_load = factor_load(**loads, **load_combo)
        factored_loads.append(factored_load)
    return min(factored_loads)


def envelope_max(result_arrays: dict) -> list[list[float]]:
    """
    Returns the maximum factored array across all factored result arrays in 'result_arrays'.

    'result_arrays': a dict of factored result arrays for an action on a specific framing member,
        keyed by load combo name. The result array values are a Nx2 array where the x-coordinates
        are in index 0 and the y-coordinates are in index 1
    """
    stacked_results = []
    for result_array in result_arrays.values():
        stacked_results.append(result_array[1])
        x_array = result_array[0]

    enveloped = []
    for idx, _ in enumerate(x_array):
        result_elems = []
        for result_array in stacked_results:
            result_elems.append(result_array[idx])
        enveloped.append(max(result_elems))
    return [x_array, enveloped]


# Fancy, alternate implementation below...(used in min, afterwards)
# def envelope_max(result_arrays: dict) -> list[list[float]]:
#     stacked_results = [result_array[1] for result_array in result_arrays.values()]
#     enveloped = [max(results_at_station) for results_at_station in zip(*stacked_results)]
#     first_key = list(result_arrays.keys())[0]
#     return [result_arrays[first_key][0], enveloped]


def envelope_min(result_arrays: dict) -> list[list[float]]:
    """
    Returns the minimum factored array across all factored result arrays in 'result_arrays'.

    'result_arrays': a dict of factored result arrays for an action on a specific framing member,
        keyed by load combo name. The result array values are a Nx2 array where the x-coordinates
        are in index 0 and the y-coordinates are in index 1
    """
    stacked_results = [result_array[1] for result_array in result_arrays.values()]
    enveloped = [
        min(results_at_station) for results_at_station in zip(*stacked_results)
    ]
    first_key = list(result_arrays.keys())[0]
    return [result_arrays[first_key][0], enveloped]
