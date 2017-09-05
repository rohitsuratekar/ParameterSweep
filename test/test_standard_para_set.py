import matplotlib.pylab as plt

from analysis.helper import *


def get_ratios(out, print_screen=False):
    pi = out[0] + out[7]
    pi4p = out[1] / pi
    pip2 = out[2] / pi
    dag = out[3] / pi
    pa = (out[4] + out[5]) / pi
    cdpdag = out[6] / pi
    if print_screen:
        print("=================")
        print("WT Output Ratios wrt Total PI")
        print("PIP2: %f, PI4P: %f, DAG: %f, PA: %f, CDPDAG: %f" % (pip2, pi4p, dag, pa, cdpdag))
    return pip2, pi4p, dag, pa, cdpdag


def mutant_ratios(wt, mt):
    print("=================")
    print("Mutant Output Ratios in given mutant")
    print("DAG: %f, PA: %f" % (mt[2] / wt[2], mt[3] / wt[3]))


with open("test_para.txt", "r") as f:
    for line in f:
        l = line.split(":", 1)
        data = json.loads(l[1].strip())
        enzymes = {}
        for key in data['Enzymes'].keys():
            enzymes[key] = Enzyme.make_with_values(key, data['Enzymes'][key])

        initial_con = get_random_concentrations(100, S_OPEN_2)

        end_time = 1000
        output = get_concentration_profile(S_OPEN_2, initial_con, enzymes, end_time, 10000)

        enzymes[E_LAZA].k *= 0.1

        output_mt = get_concentration_profile(S_OPEN_2, initial_con, enzymes, end_time, 10000)

        mutant_ratios(get_ratios(output[-1], True), get_ratios(output_mt[-1]))

        time = np.linspace(0, end_time, 10000)
        plt.plot(time, output)
        # plt.show()
