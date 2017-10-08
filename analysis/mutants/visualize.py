from matplotlib import pylab as plt

from analysis.helper import *


class DataObject:
    def __init__(self, enzymes, mutants, original):
        self.original = original
        self.enzymes = enzymes
        self.mutants = mutants
        self.laza = mutants['LAZA']
        self.rdga = mutants['RDGA']
        self.laza_pa = self.laza['PA']
        self.laza_dag = self.laza['DAG']
        self.rdga_pa = self.rdga['PA']
        self.rdga_dag = self.rdga['DAG']
        self.error = abs(self.laza_dag - 1) + abs(
            self.laza_pa - 2.5) / 2.5 + abs(self.rdga_dag - 1) + abs(
            self.rdga_pa - 1)


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
    if is_dag:
        val = 'DAG'
        plt.xlim(-10, 50)
        plt.ylim(0, 5)
        c = '#3F51B5'
        save_name = system + '_dag.png'
        wt = 1
    else:
        c = '#F44336'
        val = 'PA_{total}'
        plt.xlim(0, 4)
        plt.ylim(0, 10)
        save_name = system + '_pa.png'
        wt = 2.4

    plt.scatter(rdga, laza, marker='o', color=c, alpha=0.5)
    plt.axhline(y=wt, color='grey', linestyle='--')
    plt.axvline(x=1, color='grey', linestyle='--')
    plt.title(system)
    plt.annotate('$%s/PI_{total}$' % val, xy=(0.8, 0.95),
                 xycoords='axes fraction')
    plt.xlabel('Ratio in $rdgA^{3}$')
    plt.ylabel('Ratio in $laza^{22}$')

    if is_dag:
        plt.scatter(lowest_para.rdga_dag, lowest_para.laza_dag, marker='*',
                    color='k')
    else:
        plt.scatter(lowest_para.rdga_pa, lowest_para.laza_pa, marker='*',
                    color='k')

    plt.savefig(save_name, format='png', dpi=300)


def visualize(filename, system):
    para = get_parameter_set(filename)
    rdga_dag = []
    rdga_pa = []
    laza_dag = []
    laza_pa = []
    lowest_para = None
    for p in para:
        rdga_dag.append(p.rdga_dag)
        rdga_pa.append(p.rdga_pa)
        laza_dag.append(p.laza_dag)
        laza_pa.append(p.laza_pa)
        if lowest_para is None:
            lowest_para = p
        else:
            if p.error < lowest_para.error:
                lowest_para = p

    print_enzymes(lowest_para.enzymes)
    print(lowest_para.original)

    save_plot(rdga_dag, laza_dag, lowest_para, True, system)
    save_plot(rdga_pa, laza_pa, lowest_para, False, system)
