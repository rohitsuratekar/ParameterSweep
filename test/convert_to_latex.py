import re
from collections import defaultdict

from analysis.helper import *
from biology.component import Enzyme


def get_parameter_set(filename) -> list:
    all_para = []
    with open(filename, "r") as f:
        for line in f:
            raw = extract_enz_from_log(line)
            para = {}
            for key in raw.keys():
                para[key] = Enzyme.make_with_values(key, raw[key])
            all_para.append(para)

    return all_para


def print_table(table: defaultdict):
    print_string = "\\begin{table}\n"
    print_string += table["header"] + "}\hline\n"
    print_string += r"&&\multicolumn{%d}{c|}{\textbf{Parameter Sets}}\\\hline" % \
                    (table["header"].count("c") - 2)
    print_string += "\n"
    for k in table:
        if k is not "header":
            table[k] = table[k][:-1]
            if k[-1] == "v":
                print_string += table[k] + r"\\\cline{2-%d}" % (table[
                                                                    k].count(
                    "&") + 1) + "\n"
            else:
                print_string += table[k] + r"\\\hline" + "\n"

    print_string += "\end{tabular}\n\end{table}"
    print(print_string)


def convert_to_latex(filename):
    line_dict = defaultdict(str)
    header_put = False
    set_number = 0
    for e in get_parameter_set(filename):
        set_number += 1
        try:
            e.pop(E_IP3_PTASE)
            e.pop(E_P5TASE)
            e.pop(E_P4TASE)
        except KeyError:
            pass
        if not header_put:
            line_dict["header"] += r"\begin{tabular}{|c|c|"
            line_dict["title"] += "&& "
            header_put = True
        line_dict["header"] += "c|"
        line_dict["title"] += str(set_number) + "&"
        for enx in e:
            if len(line_dict[enx + "v"]) == 0:
                line_dict[
                    enx + "v"] += "\multirow{2}{*}{\\textbf{%s}} & v & " % enx
                line_dict[enx + "k"] += " & k &"
            line_dict[enx + "v"] += str(e[enx].v) + " &"
            line_dict[enx + "k"] += str(e[enx].k) + " &"

    print_table(line_dict)


def convert_from_latex(filename):
    all_v = defaultdict(list)
    all_k = defaultdict(list)
    last_key = ""
    with open(filename, "r") as f:
        for lines in f:
            if lines.strip().startswith("\multirow{2}{*}"):
                s = lines.strip().replace(r"\\\cline{2-8}", "")
                key = re.search(r'textbf{(.*?)}', s).group(1)
                last_key = key
                for v in s.split("&"):
                    try:
                        all_v[key].append(float(v.strip()))
                    except ValueError:
                        pass
            else:
                s = lines.replace(r"\\\hline", "").strip()
                for k in s.split("&"):
                    try:
                        all_k[last_key].append(float(k.strip()))
                    except ValueError:
                        pass

    para_sets = []
    for i in range(len(all_k[E_SOURCE])):
        s = {}
        for key in all_k.keys():
            e = Enzyme(key, v=all_v[key][i], k=all_k[key][i],
                       kinetics=MICHAELIS_MENTEN)
            s[key] = e
        para_sets.append(s)

    for key in para_sets[0].keys():
        key_string = key + "\t"
        for p in para_sets:
            key_string += str(p[key].k) + "\t"
        print(key_string)

    with open("extracted.txt", "w") as b:
        for p in para_sets:
            data = {}
            for value in p.values():
                data[value.name] = {
                    "v": round(value.v, 4),
                    "k": round(value.k, 4),
                    "kinetics": value.kinetics
                }
            save_values = {
                "Enzymes": data
            }
            print("xyz :" + json.dumps(save_values, sort_keys=True), file=b)
