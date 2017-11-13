from biology.component import Enzyme
from constants.namespace import *


def gfs(concentrations: list, feedback_substrate: str):
    """
    Get Feedback Substrate
    :param concentrations: current concentrations of all elements
    :param feedback_substrate: Substrate which is giving feedback
    :return: amount
    """
    con_list = [L_PMPI, L_PI4P, L_PIP2, L_DAG, L_PMPA, L_ERPA, L_CDPDAG,
                L_ERPI]
    try:
        con_index = con_list.index(feedback_substrate)
        return concentrations[con_index]
    except ValueError:
        return 1


def rwf(enzyme: Enzyme, substrate: float, con: list):
    """
    React with Feedback
    :param enzyme: Enzyme
    :param substrate: Substrate
    :param con: ALl Concentrations
    :return: After reaction value
    """
    return enzyme.react_with_feedback(substrate,
                                      gfs(con, enzyme.feedback_substrate))


def get_equations(c: list, time: tuple, *args) -> list:
    assert len(c) == 8, "You should provide all concentrations"
    enzyme_list = args[0]  # type: dict

    pitp = enzyme_list.get(E_PITP)  # type: Enzyme
    pip5k = enzyme_list.get(E_PIP5K)  # type: Enzyme
    plc = enzyme_list.get(E_PLC)  # type: Enzyme
    pi4k = enzyme_list.get(E_PI4K)  # type: Enzyme
    dagk = enzyme_list.get(E_DAGK)  # type: Enzyme
    laza = enzyme_list.get(E_LAZA)  # type: Enzyme
    patp = enzyme_list.get(E_PATP)  # type: Enzyme
    cds = enzyme_list.get(E_CDS)  # type: Enzyme
    pis = enzyme_list.get(E_PIS)  # type: Enzyme
    sink = enzyme_list.get(E_SINK)  # type: Enzyme
    source = enzyme_list.get(E_SOURCE)  # type: Enzyme

    pmpi, pi4p, pip2, dag, pmpa, erpa, cdpdag, erpi = c

    d_pmpi = rwf(pitp, erpi, c) - rwf(pi4k, pmpi, c)
    d_pi4p = rwf(pi4k, pmpi, c) - rwf(pip5k, pi4p, c)
    d_pip2 = rwf(pip5k, pi4p, c) - rwf(plc, pip2, c)
    d_dag = rwf(plc, pip2, c) - rwf(dagk, dag, c) + rwf(laza, pmpa, c) - rwf(
        sink, dag, c)
    d_pmpa = rwf(dagk, dag, c) - rwf(laza, pmpa, c) - rwf(patp, pmpa, c)
    d_erpa = rwf(patp, pmpa, c) - rwf(cds, erpa, c) + rwf(source, None,
                                                          c)  # equal to k value of source
    d_cdpdag = rwf(cds, erpa, c) - rwf(pis, cdpdag, c)
    d_erpi = rwf(pis, cdpdag, c) - rwf(pitp, erpi, c)

    return [d_pmpi, d_pi4p, d_pip2, d_dag, d_pmpa, d_erpa, d_cdpdag, d_erpi]
