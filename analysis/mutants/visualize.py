"""
Visualize scatter plot of between different mutant.
Need mutant sweep file created with "mutant_check.py" script
"""

from collections import defaultdict
from random import shuffle

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


def save_plot(rdga, laza, lowest_para, is_dag, system, allowed_error):
    plt.clf()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for k in rdga.keys():
        new_rdga = []
        new_laza = []
        for i in range(len(allowed_error[k])):
            if allowed_error[k][i].count(True) == len(allowed_error[k][i]):
                new_rdga.append(rdga[k][i])
                new_laza.append(laza[k][i])

        if is_dag:
            val = 'DAG'
            ax1.set_xlim(-10, 50)
            ax1.set_ylim(0, 5)
            c = '#3F51B5'
            save_name = system + "_" + k + '_dag.png'
            wt = 1
        else:
            c = '#F44336'
            val = 'PA_{total}'
            ax1.set_xlim(0, 4)
            ax1.set_ylim(0, 10)
            save_name = system + "_" + k + '_pa.png'
            wt = 2.4

        ax1.axhline(y=wt, color='grey', linestyle='--')
        ax1.axvline(x=1, color='grey', linestyle='--')

        ax1.scatter(rdga[k], laza[k], marker='o', color=c, alpha=0.5)
        ax1.scatter(new_rdga, new_laza, marker='o', color='g', alpha=0.5)
        ax1.set_title(system)
        ax1.annotate('$%s/PI_{total}$' % val, xy=(0.8, 0.95),
                     xycoords='axes fraction')
        ax1.set_xlabel('Ratio in $rdgA^{3}$')
        ax1.set_ylabel('Ratio in $laza^{22}$')

        if is_dag:
            ax1.scatter(lowest_para.rdga_dag[k], lowest_para.laza_dag[k],
                        marker='*',
                        color='k')
        else:
            ax1.scatter(lowest_para.rdga_pa[k], lowest_para.laza_pa[k],
                        marker='*',
                        color='k')

        fig1.savefig(save_name, format='png', dpi=300)


def get_boolean_value(val, expected, error) -> bool:
    return expected * (1 - error) < val <= expected * (1 + error)


def visualize(filename, system, expression_level: float, allowed_error=0.15):
    para = get_parameter_set(filename)
    rdga_dag = defaultdict(list)
    rdga_pa = defaultdict(list)
    laza_dag = defaultdict(list)
    laza_pa = defaultdict(list)
    in_allowed_error = defaultdict(list)
    lowest_para = None
    for p in para:
        for k in p.expressions:
            rdga_dag[k].append(p.rdga_dag[k])
            rdga_pa[k].append(p.rdga_pa[k])
            laza_dag[k].append(p.laza_dag[k])
            laza_pa[k].append(p.laza_pa[k])
            in_allowed_error[k].append(
                [get_boolean_value(p.rdga_dag[k], 1, allowed_error),
                 get_boolean_value(p.rdga_pa[k], 1, allowed_error),
                 get_boolean_value(p.laza_dag[k], 1, allowed_error),
                 get_boolean_value(p.laza_pa[k], 2.5, allowed_error)])
            if lowest_para is None:
                lowest_para = p
            else:
                if p.error[k] < lowest_para.error[k]:
                    lowest_para = p

    print_enzymes(lowest_para.enzymes)
    print(lowest_para.original)

    save_plot(rdga_dag, laza_dag, lowest_para, True, system, in_allowed_error)
    save_plot(rdga_pa, laza_pa, lowest_para, False, system, in_allowed_error)
    para_print = 0
    para.sort(key=lambda x: x.error[str(expression_level)])
    first_para = para[:500]
    shuffle(first_para)
    with open("top_para.txt", "w") as f:
        for p in first_para:
            print(p.original.strip(), file=f, end='\n')
            para_print += 1
            if para_print > 10:
                break


def plot_histograms(values, mutant, lipid):
    plt.clf()
    if lipid == "PA" and mutant == "$laza^{22}$":
        plt.axvline(x=2.5, color='k', linestyle='--')
    else:
        plt.axvline(x=1, color='k', linestyle='--')

    for k in values.keys():
        plt.hist(values[k], alpha=0.5, label=k + " x WT")

    ratio = '$[DAG]/PI_{total}$'

    if lipid == "PA":
        ratio = '$[PA]_{total}/PI_{total}$'

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
    fig, axs = plt.subplots(3, 3)
    axs = axs.ravel()
    para_set = get_parameter_set("analysis/mutants/top_para.txt")
    plot_number = 0
    colors = ["#FFEBEE", "#FFCDD2", "#EF9A9A", "#E57373", "#EF5350", "#F44336",
              "#E53935", "#D32F2F", "#C62828", "#B71C1C"]

    for expression in np.linspace(0.01, 0.9, num=9):
        ratio_pa = []
        ratio_dag = []
        for para in para_set:
            initial_con = get_random_concentrations(total_lipid_concentration,
                                                    system)
            wt_output = get_concentration_profile(system, initial_con,
                                                  para.enzymes,
                                                  ode_end_time, ode_slices)

            para.enzymes[E_LAZA].v *= expression
            output = get_concentration_profile(system, initial_con,
                                               para.enzymes,
                                               ode_end_time, ode_slices)
            para.enzymes[E_LAZA].v *= (1 / expression)

            ratio_pa.append(get_pa_ratio(wt_output[-1], output[-1]))
            ratio_dag.append(get_dag_ratio(wt_output[-1], output[-1]))

        axs[plot_number].scatter(ratio_dag, ratio_pa,
                                 color=colors[plot_number + 1])
        axs[plot_number].axhline(y=2.5, linestyle="--", color="grey",
                                 alpha=0.5)
        axs[plot_number].axvline(x=1, linestyle="--", color="grey",
                                 alpha=0.5)
        axs[plot_number].set_title("%.2f x WT" % expression)
        color_y_axis(axs[plot_number], "grey")
        plot_number += 1

    fig.tight_layout()
    plt.savefig("para_sensitivity.png", format='png', dpi=300,
                bbox_inches='tight')
    plt.show()
