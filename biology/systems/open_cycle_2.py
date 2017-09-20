from biology.component import Enzyme
from constants.namespace import *


def get_equations(concentrations: list, time: tuple, *args) -> list:
    assert len(concentrations) == 8, "You should provide all concentrations"
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

    pmpi, pi4p, pip2, dag, pmpa, erpa, cdpdag, erpi = concentrations

    d_pmpi = pitp.react_with(erpi) - pi4k.react_with(pmpi)
    d_pi4p = pi4k.react_with(pmpi) - pip5k.react_with(pi4p)
    d_pip2 = pip5k.react_with(pi4p) - plc.react_with(pip2)
    d_dag = plc.react_with(pip2) - dagk.react_with(dag) + laza.react_with(
        pmpa) - sink.react_with(dag)
    d_pmpa = dagk.react_with(dag) - laza.react_with(pmpa) - patp.react_with(
        pmpa)
    d_erpa = patp.react_with(pmpa) - cds.react_with(erpa) + source.react_with(
        None)  # equal to k value of source
    d_cdpdag = cds.react_with(erpa) - pis.react_with(cdpdag)
    d_erpi = pis.react_with(cdpdag) - pitp.react_with(erpi)

    return [d_pmpi, d_pi4p, d_pip2, d_dag, d_pmpa, d_erpa, d_cdpdag, d_erpi]
