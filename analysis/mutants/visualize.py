"""
Visualize scatter plot of between different mutant.
Need mutant sweep file created with "mutant_check.py" script
"""

from collections import defaultdict

from matplotlib import pylab as plt

from analysis.helper import *


class DataObject:
    def __init__(self, enzymes, mutants, original):
        self.original = original
        self.enzymes = enzymes
        self.mutants = mutants
        self.expressions = mutants.keys()
        self.laza = {}
        self.rdga = {}
        self.laza_pa = {}
        self.laza_dag = {}
        self.rdga_pa = {}
        self.rdga_dag = {}
        self.error = {}
        for k in self.expressions:
            self.laza[k] = mutants[k]['LAZA']
            self.rdga[k] = mutants[k]['RDGA']
            self.laza_pa[k] = self.laza[k]['PA']
            self.laza_dag[k] = self.laza[k]['DAG']
            self.rdga_pa[k] = self.rdga[k]['PA']
            self.rdga_dag[k] = self.rdga[k]['DAG']
            self.error[k] = abs(self.laza_dag[k] - 1) + abs(
                self.laza_pa[k] - 2.5) / 2.5 + abs(self.rdga_dag[k] - 1) + abs(
                self.rdga_pa[k] - 1)


def get_parameter_set(filename) -> list:
    all_para = []
    with open(filename, "r") as f:
        for line in f:
            raw = extract_enz_from_log(line)
            para = {}
            for key in raw.keys():
                para[key] = Enzyme.make_with_values(key, raw[key])

            all_para.append(
                DataObject(para, json.loads(line.split(":", 1)[1])[
                    "Mutants"], line))

    return all_para


def print_enzymes(all_enzymes):
    for e in all_enzymes:
        print(
            "%s\t%s\t%s" % (all_enzymes[e].properties.get('name').upper(),
                            all_enzymes[
                                e].properties.get('v'),
                            all_enzymes[e].properties.get('k')))


def save_plot(rdga, laza, lowest_para, is_dag, system):
    plt.clf()
    for k in rdga.keys():
        if is_dag:
            val = 'DAG'
            plt.xlim(-10, 50)
            plt.ylim(0, 5)
            c = '#3F51B5'
            save_name = system + "_" + k + '_dag.png'
            wt = 1
        else:
            c = '#F44336'
            val = 'PA_{total}'
            plt.xlim(0, 4)
            plt.ylim(0, 10)
            save_name = system + "_" + k + '_pa.png'
            wt = 2.4

        plt.scatter(rdga[k], laza[k], marker='o', color=c, alpha=0.5)
        plt.axhline(y=wt, color='grey', linestyle='--')
        plt.axvline(x=1, color='grey', linestyle='--')
        plt.title(system)
        plt.annotate('$%s/PI_{total}$' % val, xy=(0.8, 0.95),
                     xycoords='axes fraction')
        plt.xlabel('Ratio in $rdgA^{3}$')
        plt.ylabel('Ratio in $laza^{22}$')

        if is_dag:
            plt.scatter(lowest_para.rdga_dag[k], lowest_para.laza_dag[k],
                        marker='*',
                        color='k')
        else:
            plt.scatter(lowest_para.rdga_pa[k], lowest_para.laza_pa[k],
                        marker='*',
                        color='k')

        plt.savefig(save_name, format='png', dpi=300)


def visualize(filename, system):
    para = get_parameter_set(filename)
    rdga_dag = defaultdict(list)
    rdga_pa = defaultdict(list)
    laza_dag = defaultdict(list)
    laza_pa = defaultdict(list)
    lowest_para = None
    top_para = []
    for p in para:
        for k in p.expressions:
            rdga_dag[k].append(p.rdga_dag[k])
            rdga_pa[k].append(p.rdga_pa[k])
            laza_dag[k].append(p.laza_dag[k])
            laza_pa[k].append(p.laza_pa[k])
            if lowest_para is None:
                top_para.append(p)
                lowest_para = p
            else:
                if p.error[k] < lowest_para.error[k]:
                    lowest_para = p
                    top_para.insert(0, p)

    print_enzymes(lowest_para.enzymes)
    print(lowest_para.original)

    save_plot(rdga_dag, laza_dag, lowest_para, True, system)
    save_plot(rdga_pa, laza_pa, lowest_para, False, system)
    para_print = 0
    with open("top_para.txt", "w") as f:
        for p in top_para:
            print(p.original.strip(), file=f, end='\n')
            para_print += 1
            if para_print > 10:
                break


def plot_histograms(values, mutant, lipid):
    plt.clf()
    for k in values.keys():
        plt.hist(values[k], alpha=0.5, label=k + " x WT")

    ratio = '$[DAG]/PI_{total}$'

    if lipid == "PA":
        ratio = '$[PA]_{total}/PI_{total}$'

    if lipid == "PA" and mutant == "$laza^{22}$":
        plt.axvline(x=2.5, color='k', linestyle='--')
    else:
        plt.axvline(x=1, color='k', linestyle='--')

    plt.xlabel('%s ratio in %s/WT' % (ratio, mutant))
    plt.ylabel('Frequency (on log scale)')
    plt.yscale('log')
    plt.legend(title="Expression level", loc='center left',
               bbox_to_anchor=(1, 0.5))
    file_name = "sensitivity_" + lipid + "_" + mutant + ".png"
    file_name = file_name.replace("^", "").replace("$", "").replace("{", "") \
        .replace("}", "")
    plt.savefig(file_name, format='png', dpi=300, bbox_inches='tight')


def visualize_sensitivity(filename, system):
    para = get_parameter_set(filename)
    rdga_dag = defaultdict(list)
    rdga_pa = defaultdict(list)
    laza_dag = defaultdict(list)
    laza_pa = defaultdict(list)
    for p in para:
        for k in p.expressions:
            rdga_dag[k].append(p.rdga_dag[k])
            rdga_pa[k].append(p.rdga_pa[k])
            laza_dag[k].append(p.laza_dag[k])
            laza_pa[k].append(p.laza_pa[k])

    plot_histograms(rdga_dag, "$rdgA^3$", "DAG")
    plot_histograms(rdga_pa, "$rdgA^3$", "PA")
    plot_histograms(laza_dag, "$laza^{22}$", "DAG")
    plot_histograms(laza_pa, "$laza^{22}$", "PA")


def get_pa_ratio(wt, mt):
    return ((mt[4] + mt[5]) / (mt[0] + mt[7])) / (
        (wt[4] + wt[5]) / (wt[0] + wt[7]))


def get_dag_ratio(wt, mt):
    return (mt[3] / (mt[0] + mt[7])) / (wt[3] / (wt[0] + wt[7]))


def color_y_axis(ax, color):
    """Color your axes."""
    for t in ax.get_yticklabels():
        t.set_color(color)
    return None


def single_para_sensitivity(system: str):
    c1 = '#F44336'
    c2 = '#3F51B5'
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    for para in get_parameter_set("analysis/mutants/top_para.txt"):
        initial_con = get_random_concentrations(total_lipid_concentration,
                                                system)
        wt_output = get_concentration_profile(system, initial_con,
                                              para.enzymes,
                                              ode_end_time, ode_slices)

        expression = np.arange(0.1, 1.0, 0.03)
        ratio_pa = []
        ratio_dag = []
        for i in expression:
            para.enzymes[E_LAZA].v *= i
            output = get_concentration_profile(system, initial_con,
                                               para.enzymes,
                                               ode_end_time, ode_slices)
            para.enzymes[E_LAZA].v *= (1 / i)

            ratio_pa.append(get_pa_ratio(wt_output[-1], output[-1]))
            ratio_dag.append(get_dag_ratio(wt_output[-1], output[-1]))

        ax1.plot(expression, ratio_pa, color=c1)

        ax2.plot(expression, ratio_dag, color=c2)
        ax2.set_ylabel('$DAG/PI_{total}$')

    ax1.set_xlabel("$LAZA$ activity (wrt WT)")
    ax1.set_ylabel('$PA_{total}/PI_{total}$')
    color_y_axis(ax1, c1)
    color_y_axis(ax2, c2)

    plt.savefig("para_sensitivity.png", format='png', dpi=300,
                bbox_inches='tight')
