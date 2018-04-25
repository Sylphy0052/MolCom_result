import matplotlib.pyplot as plt
from matplotlib import ticker
from enum import Enum
import os
import numpy as np
import sys

COLOR_LIST = ['r', 'g', 'b', 'm', 'c', 'y', 'k', 'r', 'g', 'b', 'm', 'c', 'y', 'k']
STYLE_LIST = ['-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--', '-.', ':']

class Mode(Enum):
    NORMAL = 0
    DC = 1
    RTO = 2
    DECOMPOSING = 3
    ADJUST = 4
    FEC = 5

class Xvalue(Enum):
    DISTANCE = 0
    DUPLICATION = 1
    DC = 2
    RWT = 3
    RTO = 4
    DECOMPOSING = 5
    NONE = 6
    NUMMESSAGE = 7
    DECOMPOSINGMETHOD = 8
    RATE = 9
    PACKETNUM = 10

class Yvalue(Enum):
    MEAN = 0
    MEDIAN = 1
    TXRXMEAN = 2
    JITTER = 3
    COLLISION = 4
    FAILURE = 5
    RWT = 6
    MOLECULENUM = 7
    DECOMPOSINGNUM = 8
    RETRANSMITFAILURE = 9


def checkMode(data_dict, all_file_list):
    for file_name in all_file_list:
        data = data_dict[file_name]
        if data.config["FEC"] is not None:
            return Mode.FEC
        if data.params["decomposing"]:
            return Mode.DECOMPOSING
        if data.params["info_adjust"] != 0:
            return Mode.ADJUST
        if "rto_type" in data.params.keys():
            return Mode.RTO
        if data.params["step_length"] != 1:
            return Mode.DC

    return Mode.NORMAL

def checkOption(data_dict, all_file_list):
    is_ptime = True
    is_coll = True
    for file_name in all_file_list:
        data = data_dict[file_name]
        if not data.is_ptime:
            is_ptime = False
        if not data.is_coll:
            is_coll = False
    return [is_ptime, is_coll]

def drawBarGraph(X, Y, X_labels, y_label, fig_name):
    for x, y in zip(X, Y):
        y = round(y, 1)
        plt.text(x, y, y, ha='center', va='bottom')

    plt.ylabel(y_label)
    plt.bar(X, Y, color=COLOR_LIST, tick_label=X_labels, width=0.5)

    plt.savefig(fig_name, dpi=90, bbox_inches="tight", pad_inches=0.0)
    plt.close('all')

def drawTwoGraph(X, Y1, Y2, ax1_label, ax2_label, X_label, Y1_label, Y2_label, location, fig_name):
    fig, ax1 = plt.subplots()
    ln1 = ax1.plot(X, Y1, color=COLOR_LIST[0], label=ax1_label)
    ax2 = ax1.twinx()
    ln2 = ax2.plot(X, Y2, color=COLOR_LIST[1], label=ax2_label, linestyle=STYLE_LIST[0])

    ax1.set_xlabel(X_label)
    ax1.set_ylabel(Y1_label)
    ax2.set_ylabel(Y2_label)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc=location)

    plt.grid(True)

    plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    plt.gca().ticklabel_format(style="sci",  axis="x",scilimits=(0,0))

    plt.savefig(fig_name, dpi=90, bbox_inches="tight", pad_inches=0.0)
    plt.close('all')

def drawBySimulation(data_dict, all_file_list):
    path = "./result_img"
    if not os.path.isdir(path):
        os.makedirs(path)
    for file_name in all_file_list:
        drawProbCumprobEachSimulation(data_dict[file_name])
        drawMeanEachSimulation(data_dict[file_name])
        drawAdjustGraphBySimulation(data_dict[file_name])
        drawCollisionGraphBySimulation(data_dict[file_name])
        drawRetransmissionGraphBySimulation(data_dict[file_name])

def drawRetransmissionGraphBySimulation(data):
    if not os.path.isfile(data.retransmission_file_name):
        return
    path = "./retransmission_graph"
    if not os.path.isdir(path):
        os.makedirs(path)

    steps = data.input_data.collisions

    retransmission_file_name = data.retransmission_file_name
    fig_name = "./retransmission_graph/" + retransmission_file_name.split(".txt")[0].split("/")[2] + ".png"
    if not os.path.isfile(retransmission_file_name):
        return

    retransmission_num = 0
    retransmissions = []
    start_pos = 0
    steps = data.input_data.retransmit_steps
    for step in steps:
        if step == 0:
            start_pos += 1
            continue
        retransmission_num += 1 / len(data.input_data.steps)
        retransmissions.append(retransmission_num)

    plt.plot(steps[start_pos:], retransmissions, color=COLOR_LIST[0], markersize="5")

    plt.xlabel('RTT (s)')
    plt.ylabel('Retransmission Num')

    plt.grid(True)
    plt.legend(loc='upper left')

    plt.savefig(fig_name)
    plt.close('all')

def drawCollisionGraphBySimulation(data):
    collision_file_name = data.collision_file_name
    fig_name = "./collision_graph/" + collision_file_name.split(".txt")[0].split("/")[2] + ".png"
    if not os.path.isfile(collision_file_name):
        return

    path = "./collision_graph"
    if not os.path.isdir(path):
        os.makedirs(path)

    steps = data.input_data.collisions
    collision_num = 0
    collisions = []
    for step in steps:
        collision_num += (1 / len(data.input_data.steps))
        collisions.append(collision_num)

    plt.plot(steps, collisions, color=COLOR_LIST[0], markersize="5")

    plt.xlabel('RTT (s)')
    plt.ylabel('Collision Num')

    plt.grid(True)
    plt.legend(loc='upper left')

    plt.savefig(fig_name)
    plt.close('all')

def drawAdjustGraphBySimulation(data):
    if not os.path.isfile(data.adjust_file_name):
        return
    path = "./adjust_graph"
    if not os.path.isdir(path):
        os.makedirs(path)

    steps = data.input_data.adjust_steps
    info_adjust_num = data.input_data.info_adjust_num
    ack_adjust_num = data.input_data.ack_adjust_num

    adjust_file_name = data.adjust_file_name
    fig_name = "./adjust_graph/" + adjust_file_name.split(".txt")[0].split("/")[2] + ".png"

    if not os.path.isfile(adjust_file_name):
        return

    plt.plot(steps, info_adjust_num, color=COLOR_LIST[0], label="INFO Num", markersize="5")
    plt.plot(steps, ack_adjust_num, color=COLOR_LIST[1], label="ACK Num", markersize="5")

    plt.xlabel('RTT (s)')
    plt.ylabel('Molecule Num')

    plt.grid(True)
    plt.legend(loc='upper left')

    plt.savefig(fig_name)
    plt.close('all')

def drawProbCumprobEachSimulation(data):
    fig_name = "./result_img/" + data.file_name + "_prob_cumprob.png"
    data.setProbCumprobGraph(fig_name)
    plot_data = data.plot_data

    X = [0]
    X.extend(x[1] + 1 for x in plot_data.plot_range)
    Y1 = [0]
    Y1.extend(plot_data.prob)
    Y2 = [0]
    Y2.extend(plot_data.cum_prob)

    drawTwoGraph(X, Y1, Y2, "Probability", "Cumulative Probability", "RTT (s)", "Probability of RTT", "Cumulative Probability of RTT", "right", fig_name)

def drawMeanEachSimulation(data):
    fig_name = "./result_img/" + data.file_name + "_mean.png"
    data.setMeanGraph(fig_name)
    X_labels = ["Mean"]
    Y = [float(data.params["mean"])]
    # if data.is_ptime:
    #     X_labels.append("TxRx Mean")
    #     Y.append(float(data.params["txrx_mean"]))
    # X_labels.append("Analytical")
    # Y.append(float(data.params["analytical_rtt"]))
    X = range(len(Y))

    drawBarGraph(X, Y, X_labels, "Mean RTT (s)", fig_name)

def rwtSort(X, Y, X_labels):
    XY = [{} for i in range(len(X_labels))]
    for i in range(len(Y)):
        x = X[i]
        y = Y[i]
        for j in range(len(X)):
            XY[i][x[j]] = y[j]
    retX = []
    retY = []
    for i in range(len(X_labels)):
        x = []
        y = []
        for k, v in sorted(XY[i].items()):
            x.append(k)
            y.append(v)
        retX.append(x)
        retY.append(y)
    return [retX, retY]

def drawLineGraph(X, Y, X_labels, x_val, y_val, fig_name, is_math):
    if x_val is Xvalue.RWT:
        X, Y = rwtSort(X, Y, X_labels)
    # elif len(X) != len(Y):
    #     X = [X for i in range(len(Y))]
    elif np.array(X).shape != np.array(Y).shape:
        X = [X for i in range(len(Y))]

    if x_val is Xvalue.RTO and y_val is Yvalue.RWT:
        # print(X, Y)
        X = [x[1:] for x in X]
        Y = [y[1:] for y in Y]
        # X_labels = X_labels[1:]
        # print(X, Y)

    # print("\n---X---\n")
    # for x in X:
    #     print(x)
    # print("\n---Y---\n")
    # for y in Y:
    #     print(y)
    for i in range(len(Y)):
        plt.plot(X[i], Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

    plt.xlabel(defineXlabel(x_val))
    plt.ylabel(defineYlabel(y_val))
    plt.grid(True)

    if x_val is Xvalue.DISTANCE:
        plt.xticks(X[0])
        plt.legend(loc='upper left')
    elif x_val is Xvalue.DUPLICATION:
        plt.xticks(X[0])
        plt.legend(loc='upper right')
    elif x_val is Xvalue.DC:
        x = [round(i, 1) for i in X[0]]
        plt.xticks(x)
        plt.legend(loc='upper right')
    elif x_val is Xvalue.RTO:
        x = ["RTO-{}".format(i) for i in X[0]]
        plt.xticks(X[0], x)
        plt.legend(loc='upper right')
    elif x_val is Xvalue.DECOMPOSINGMETHOD:
        x = ["Decomposing{}".format(i) for i in X[0]]
        plt.xticks(X[0], x)
        plt.legend(loc='upper right')
    elif x_val is Xvalue.PACKETNUM:
        plt.xticks(X[0])
        plt.legend(loc='upper right')
    elif x_val is Xvalue.RATE:
        plt.xticks(X[0])
        plt.legend(loc='upper right')
    else:
        plt.legend(loc='upper right')

    if y_val is Yvalue.FAILURE:
        plt.legend(loc='lower right')

    # if x_val is Xvalue.DUPLICATION and y_val is Yvalue.COLLISION:
    #     plt.legend(loc='upper left')

    if x_val is Xvalue.RTO and y_val is Yvalue.RETRANSMITFAILURE:
        plt.ylim([-5, 80])

    plt.locator_params(axis='y',min_n_ticks=4)
    if is_math:
        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

    plt.savefig(fig_name)
    plt.close('all')

def drawRegressionGraph(X, Y, X_labels, x_val, y_val, fig_name, is_math, dir_path):
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)

    if x_val is Xvalue.RWT:
        X, Y = rwtSort(X, Y, X_labels)
    elif np.array(X).shape != np.array(Y).shape:
        X = [X for i in range(len(Y))]

    for i in range(len(Y)):
        y = np.poly1d(np.polyfit(X[i], Y[i], len(X[i]) - 1))(X[i])
        y1 = np.poly1d(np.polyfit(X[i], Y[i], len(X[i])-1))

        f_c = ": $"
        for j in range(len(y1.c)):
            num_str = ""

            if len(y1.c) - j == 2:
                num_str = "{0:+g}x".format(y1.c[j], 3)
            elif len(y1.c) - j == 1:
                num_str = "{0:+g}".format(y1.c[j], 3)
            else:
                num_str = "{0:+g}x^{1}".format(y1.c[j], len(y1.c) - 1 - j)

            f_c += num_str

        f_c += '$'

        X_labels[i] += f_c
        X_labels[i] = repr(X_labels[i])
        plt.plot(X[i], Y[i], COLOR_LIST[i] + "o")
        plt.plot(X[i], y, color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])
        # plt.plot(X[i], Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

    plt.xlabel(defineXlabel(x_val))
    plt.ylabel(defineYlabel(y_val))
    plt.grid(True)
    if x_val is Xvalue.DISTANCE:
        plt.xticks(X[0])
    elif x_val is Xvalue.DUPLICATION:
        plt.xticks(X[0])
    elif x_val is Xvalue.DC:
        x = [round(i, 1) for i in X[0]]
        plt.xticks(x)
    elif x_val is Xvalue.RTO:
        x = ["RTO-{}".format(i) for i in X[0]]
        plt.xticks(X[0], x)
    else:
        plt.legend(loc='upper right')

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.locator_params(axis='y',min_n_ticks=4)
    if is_math:
        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

    plt.savefig(fig_name, bbox_inches='tight')
    plt.close('all')


def defineDirectoryPath(x_val, y_val, each_val):
    dir_path = "compare"
    dir_path += "_" + y_val.name
    dir_path += "_by_" + x_val.name
    if each_val != Xvalue.NONE:
        dir_path += "_each_" + each_val.name
    dir_path = dir_path.lower()
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
    return dir_path

def classifyDistance(data_dict, all_file_list):
    r_list = []
    for file_name in all_file_list:
        r = data_dict[file_name].params["r"]
        if not r in r_list:
            r_list.append(r)

    file_list = [[] for i in range(len(r_list))]
    for file_name in all_file_list:
        r = data_dict[file_name].params["r"]
        for i in range(len(r_list)):
            if r == r_list[i]:
                file_list[i].append(file_name)
    return file_list


def classifyDuplication(data_dict, all_file_list):
    d_list = []
    for file_name in all_file_list:
        d = data_dict[file_name].params["duplication"]
        if not d in d_list:
            d_list.append(d)

    file_list = [[] for i in range(len(d_list))]
    for file_name in all_file_list:
        d = data_dict[file_name].params["duplication"]
        for i in range(len(d_list)):
            if d == d_list[i]:
                file_list[i].append(file_name)
    return file_list

def classifyDistance(data_dict, all_file_list):
    r_list = []
    for file_name in all_file_list:
        r = data_dict[file_name].params["r"]
        if not r in r_list:
            r_list.append(r)

    file_list = [[] for i in range(len(r_list))]
    for file_name in all_file_list:
        r = data_dict[file_name].params["r"]
        for i in range(len(r_list)):
            if r == r_list[i]:
                file_list[i].append(file_name)
    return file_list

def classifyRate(data_dict, all_file_list):
    r_list = []
    for file_name in all_file_list:
        data = data_dict[file_name].params["FEC_rate"]
        if not data in r_list:
            r_list.append(data)

    file_list = [[] for i in range(len(r_list))]
    for file_name in all_file_list:
        data = data_dict[file_name].params["FEC_rate"]
        for i in range(len(r_list)):
            if data == r_list[i]:
                file_list[i].append(file_name)

    return file_list

def classifyPacket(data_dict, all_file_list):
    r_list = []
    for file_name in all_file_list:
        data = data_dict[file_name].params["FEC_packet"]
        if not data in r_list:
            r_list.append(data)

    file_list = [[] for i in range(len(r_list))]
    for file_name in all_file_list:
        data = data_dict[file_name].params["FEC_packet"]
        for i in range(len(r_list)):
            if data == r_list[i]:
                file_list[i].append(file_name)
    return file_list

def classifyFileList(data_dict, all_file_list, x_val, each_val):
    file_list = []
    is_distance = True
    # is_duplication = True
    is_rate = True
    is_packet = True

    if x_val is Xvalue.DISTANCE or each_val is Xvalue.DISTANCE:
        is_distance = False
    # if x_val is Xvalue.DUPLICATION or each_val is Xvalue.DUPLICATION:
    #     is_duplication = False

    if x_val is Xvalue.RATE or each_val is Xvalue.RATE:
        is_rate = False
    if x_val is Xvalue.PACKETNUM or each_val is Xvalue.PACKETNUM:
        is_packet = False

    # if not is_distance and not is_duplication:
    #     file_list = [all_file_list]
    # elif is_distance and not is_duplication:
    #     file_list = classifyDistance(data_dict, all_file_list)
    # elif not is_distance and is_duplication:
    #     file_list = classifyDuplication(data_dict, all_file_list)
    if is_distance:
        file_list = classifyDistance(data_dict, all_file_list)
    elif is_rate:
        file_list = classifyRate(data_dict, all_file_list)
        # for a in file_list:
        #     print(a)
    elif is_packet:
        file_list = classifyPacket(data_dict, all_file_list)

    return file_list

def defineFigureName(data_dict, file_list, x_val, each_val):
    fig_name = ""
    is_distance = True
    # is_duplication = True
    is_rate = True
    is_packet = True

    if x_val is Xvalue.DISTANCE or each_val is Xvalue.DISTANCE:
        is_distance = False
    # if x_val is Xvalue.DUPLICATION or each_val is Xvalue.DUPLICATION:
    #     is_duplication = False
    if x_val is Xvalue.RATE or each_val is Xvalue.RATE:
        is_rate = False
    if x_val is Xvalue.PACKETNUM or each_val is Xvalue.PACKETNUM:
        is_packet = False

    data = data_dict[file_list[0]]
    if is_distance:
        fig_name += "TxRx{}_".format(data.params["r"])
    # if is_duplication:
    #     fig_name += "{}-{}_".format(data.params["info_arq"], data.params["ack_arq"])
    if is_rate:
        fig_name +="Rate{}_".format(data.params["FEC_rate"])
    if is_packet:
        fig_name +="Packet{}_".format(data.params["FEC_packet"])

    fig_name += "{}-{}.png".format(data.params["info_type"], data.params["ack_type"])
    # fig_name += ".png"
    return fig_name

def defineXlabels(data_dict, file_list, each_val):
    X_labels = []
    if each_val is Xvalue.DUPLICATION:
        for file_name in file_list:
            data = data_dict[file_name]
            label = "SW-ARQ{}-{}".format(data.params["info_arq"], data.params["ack_arq"])
            if not label in X_labels:
                X_labels.append(label)

    return X_labels

def getX(data, x_val):
    if x_val is Xvalue.DISTANCE:
        return data.params["r"]
    elif x_val is Xvalue.DUPLICATION:
        return data.params["duplication"]
    elif x_val is Xvalue.DC:
        return data.params["step_length"]
    elif x_val is Xvalue.RWT:
        return data.params["retransmitWaitTime"]
    elif x_val is Xvalue.RTO:
        return int(data.params["rto_type"])
    elif x_val is Xvalue.DECOMPOSING:
        return data.params["decomposing"]
    elif x_val is Xvalue.NUMMESSAGE:
        return data.params["message_num"]
    elif x_val is Xvalue.DECOMPOSINGMETHOD:
        return data.params["decomposing"]
    elif x_val is Xvalue.RATE:
        return data.params["FEC_rate"]
    elif x_val is Xvalue.PACKETNUM:
        return data.params["FEC_packet"]
    else:
        return 0

def getY(data, y_val):
    if y_val is Yvalue.MEAN:
        return data.params["mean"]
    elif y_val is Yvalue.MEDIAN:
        return data.params["median"]
    elif y_val is Yvalue.TXRXMEAN:
        return data.params["txrx_mean"]
    elif y_val is Yvalue.JITTER:
        return data.params["std"]
    elif y_val is Yvalue.COLLISION:
        return data.params["coll"]
    elif y_val is Yvalue.FAILURE:
        return data.params["failure"]
    elif y_val is Yvalue.RWT:
        return data.params["retransmitWaitTime"]
    elif y_val is Yvalue.MOLECULENUM:
        info_num = int(data.params["last_info_num"])
        ack_num = int(data.params["last_ack_num"])
        return (info_num + ack_num) / 2
    elif y_val is Yvalue.DECOMPOSINGNUM:
        return data.params["decomposing_num"] / len(data.input_data.steps)
    elif y_val is Yvalue.RETRANSMITFAILURE:
        return data.params["retransmit_failure"]
    else:
        return 0

def getLabel(data, each_val):
    if each_val is Xvalue.DUPLICATION:
        return "SW-ARQ{}-{}".format(data.params["info_arq"], data.params["ack_arq"])
    elif each_val is Xvalue.DISTANCE:
        return "d={}".format(data.params["r"])
    elif each_val is Xvalue.RTO:
        return "RTO{}".format(data.params["rto_type"])
    elif each_val is Xvalue.DECOMPOSING:
        if data.params["decomposing"]:
            return "Decomposing"
        else:
            return "Non Decomposing"
    elif each_val is Xvalue.NUMMESSAGE:
        return "MessageNum={}".format(data.params["message_num"])
    elif each_val is Xvalue.DECOMPOSINGMETHOD:
        return "Decomposing{}".format(data.params["decomposing"])
    elif each_val is Xvalue.RATE:
        return "Rate{}".format(data.params["FEC_rate"])
    elif each_val is Xvalue.PACKETNUM:
        return "FEC Packet={}".format(data.params["FEC_packet"])
    else:
        return ""

def getXRWT(data_dict, file_list, x_val, each_val, X_labels):
    X = [[] for i in range(len(X_labels))]
    for file_name in file_list:
        data = data_dict[file_name]
        label = getLabel(data, each_val)
        for i in range(len(X_labels)):
            if label == X_labels[i]:
                X[i].append(getX(data_dict[file_name], x_val))
                break

    return X

def getInfo(data_dict, file_list, x_val, y_val, each_val):
    X = []
    X_labels = []
    for file_name in file_list:
        data = data_dict[file_name]
        label = getLabel(data, each_val)
        if not label in X_labels:
            X_labels.append(label)
        x = getX(data_dict[file_name], x_val)
        if not x in X:
            X.append(x)

    if x_val is Xvalue.RWT:
        X = getXRWT(data_dict, file_list, x_val, each_val, X_labels)
    Y = [[] for i in range(len(X_labels))]
    for file_name in file_list:
        data = data_dict[file_name]
        label = getLabel(data, each_val)
        for i in range(len(X_labels)):
            if label == X_labels[i]:
                Y[i].append(getY(data_dict[file_name], y_val))
                break

    return [X, Y, X_labels]

def defineXlabel(x_val):
    if x_val is Xvalue.DISTANCE:
        return "Tx-Rx distance (um)"
    elif x_val is Xvalue.DUPLICATION:
        return "Duplication level (n)"
    elif x_val is Xvalue.DC:
        return "Diffusion Coefficient"
    elif x_val is Xvalue.RWT:
        return "RTO (s)"
    elif x_val is Xvalue.RTO:
        return "RTO schemes"
    elif x_val is Xvalue.DECOMPOSINGMETHOD:
        return "Decomposing schemes"
    elif x_val is Xvalue.RATE:
        return "FEC Rate"
    elif x_val is Xvalue.PACKETNUM:
        return "FEC PacketNum"
    else:
        return ""

def defineYlabel(y_val):
    if y_val is Yvalue.MEAN:
        return "Mean RTT (s)"
    elif y_val is Yvalue.MEDIAN:
        return "Median RTT (s)"
    elif y_val is Yvalue.TXRXMEAN:
        return "Mean RTT (s)"
    elif y_val is Yvalue.JITTER:
        return "Jitter of RTT"
    elif y_val is Yvalue.COLLISION:
        return "Average number of collisions"
    elif y_val is Yvalue.FAILURE:
        return "Transmission failure rate"
    elif y_val is Yvalue.RWT:
        return "RTO (s)"
    elif y_val is Yvalue.MOLECULENUM:
        return "Molecule Num"
    elif y_val is Yvalue.DECOMPOSINGNUM:
        return "Decomposing Num"
    elif y_val is Yvalue.RETRANSMITFAILURE:
        return "Transmission failure rate"
    else:
        return 0

def drawCompareGraph(data_dict, all_file_list, x_val, y_val, each_val):
    dir_path = defineDirectoryPath(x_val, y_val, each_val)
    print(dir_path)
    classify_list = classifyFileList(data_dict, all_file_list, x_val, each_val)
    for file_list in classify_list:
        fig_name = dir_path + '/' + defineFigureName(data_dict, file_list, x_val, each_val)
        X, Y, X_labels = getInfo(data_dict, file_list, x_val, y_val, each_val)
        is_math = False
        if max(max(Y)) > 10 ** 5:
            is_math = True
        drawLineGraph(X, Y, X_labels, x_val, y_val, fig_name, is_math)
        reg_dir_path = "regression_" + dir_path
        fig_name = reg_dir_path + '/' + defineFigureName(data_dict, file_list, x_val, each_val)
        drawRegressionGraph(X, Y, X_labels, x_val, y_val, fig_name, is_math, reg_dir_path)

def drawNormalGraph(data_dict, all_file_list, is_ptime, is_coll):
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MEAN, Xvalue.DUPLICATION)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DUPLICATION, Yvalue.MEAN, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MEDIAN, Xvalue.DUPLICATION)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DUPLICATION, Yvalue.MEDIAN, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.JITTER, Xvalue.DUPLICATION)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DUPLICATION, Yvalue.JITTER, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.FAILURE, Xvalue.DUPLICATION)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DUPLICATION, Yvalue.FAILURE, Xvalue.DISTANCE)

    if is_coll:
        drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.COLLISION, Xvalue.DUPLICATION)
        drawCompareGraph(data_dict, all_file_list, Xvalue.DUPLICATION, Yvalue.COLLISION, Xvalue.DISTANCE)

    drawEnergyByDuplicationEachDistance()

def drawDCGraph(data_dict, all_file_list, is_ptime, is_coll):
    drawCompareGraph(data_dict, all_file_list, Xvalue.DC, Yvalue.MEDIAN, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DC, Yvalue.MEAN, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DC, Yvalue.JITTER, Xvalue.DISTANCE)

def drawRTOGraph(data_dict, all_file_list, is_ptime, is_coll):
    drawCompareGraph(data_dict, all_file_list, Xvalue.RTO, Yvalue.RWT, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.RTO, Yvalue.MEDIAN, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.RTO, Yvalue.MEAN, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.RTO, Yvalue.JITTER, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MEDIAN, Xvalue.RTO)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MEAN, Xvalue.RTO)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.JITTER, Xvalue.RTO)
    drawCompareGraph(data_dict, all_file_list, Xvalue.RWT, Yvalue.MEDIAN, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.RWT, Yvalue.MEAN, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.RWT, Yvalue.JITTER, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.RTO, Yvalue.RETRANSMITFAILURE, Xvalue.DISTANCE)

def checkDecomposingMethod(data_dict, all_file_list):
    decomposings = []
    for file_name in all_file_list:
        data = data_dict[file_name]
        decomposing_method = data.params["decomposing"]
        if not decomposing_method in decomposings:
            decomposings.append(decomposing_method)

    return len(decomposings) != 1

def drawDecomposingGraph(data_dict, all_file_list, is_ptime, is_coll):
    if checkDecomposingMethod(data_dict, all_file_list):
        drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MEAN, Xvalue.DECOMPOSINGMETHOD)
        drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MEDIAN, Xvalue.DECOMPOSINGMETHOD)
        drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.JITTER, Xvalue.DECOMPOSINGMETHOD)
        drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.DECOMPOSINGNUM, Xvalue.DECOMPOSINGMETHOD)
        drawCompareGraph(data_dict, all_file_list, Xvalue.DECOMPOSINGMETHOD, Yvalue.MEAN, Xvalue.DISTANCE)
        drawCompareGraph(data_dict, all_file_list, Xvalue.DECOMPOSINGMETHOD, Yvalue.MEDIAN, Xvalue.DISTANCE)
        drawCompareGraph(data_dict, all_file_list, Xvalue.DECOMPOSINGMETHOD, Yvalue.JITTER, Xvalue.DISTANCE)
        drawCompareGraph(data_dict, all_file_list, Xvalue.DECOMPOSINGMETHOD, Yvalue.DECOMPOSINGNUM, Xvalue.DISTANCE)
    else:
        drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MEAN, Xvalue.DECOMPOSING)
        drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MEDIAN, Xvalue.DECOMPOSING)
        drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.JITTER, Xvalue.DECOMPOSING)

def drawAdjustGraph(data_dict, all_file_list, is_ptime, is_coll):
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MOLECULENUM, Xvalue.RTO)

def drawFECGraph(data_dict, all_file_list, is_ptime, is_coll):
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MEAN, Xvalue.RATE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MEDIAN, Xvalue.RATE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.JITTER, Xvalue.RATE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MEAN, Xvalue.PACKETNUM)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.MEDIAN, Xvalue.PACKETNUM)
    drawCompareGraph(data_dict, all_file_list, Xvalue.DISTANCE, Yvalue.JITTER, Xvalue.PACKETNUM)

    drawCompareGraph(data_dict, all_file_list, Xvalue.RATE, Yvalue.MEAN, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.RATE, Yvalue.MEDIAN, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.RATE, Yvalue.JITTER, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.RATE, Yvalue.MEAN, Xvalue.PACKETNUM)
    drawCompareGraph(data_dict, all_file_list, Xvalue.RATE, Yvalue.MEDIAN, Xvalue.PACKETNUM)
    drawCompareGraph(data_dict, all_file_list, Xvalue.RATE, Yvalue.JITTER, Xvalue.PACKETNUM)

    drawCompareGraph(data_dict, all_file_list, Xvalue.PACKETNUM, Yvalue.MEAN, Xvalue.RATE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.PACKETNUM, Yvalue.MEDIAN, Xvalue.RATE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.PACKETNUM, Yvalue.JITTER, Xvalue.RATE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.PACKETNUM, Yvalue.MEAN, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.PACKETNUM, Yvalue.MEDIAN, Xvalue.DISTANCE)
    drawCompareGraph(data_dict, all_file_list, Xvalue.PACKETNUM, Yvalue.JITTER, Xvalue.DISTANCE)



def drawEnergyByDuplicationEachDistance():
    dir_path = "./compare_energy_by_duplication_each_distance/"
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)

    fig_name = dir_path + "/PASSIVE-PASSIVE.png"

    X = [1, 10, 20, 30, 50, 100]
    Y = []
    base_pairs = 52000
    for i in range(len(X)):
        # D=0.5 -> base pairs = 52000
        Y.append(X[i] * 45.6 * 2 * base_pairs)

    plt.plot(X, Y, color=COLOR_LIST[0], marker='o', markersize="5")
    plt.xlabel('Duplication level (n)')
    plt.ylabel('Energy consumption (kJ/mol)')
    plt.xticks(X, [str(x) for x in X])
    plt.grid(True)
    plt.legend(loc='upper right')

    plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

    plt.savefig(fig_name)
    plt.close('all')

# distance
def classifyData(data_dict, all_file_list):
    return_list = [[] for i in range(2)]
    for name in all_file_list:
        data = data_dict[name]
        r = int(data.params["r"])
        if r == 30:
            return_list[0].append(name)
        elif r == 50:
            return_list[1].append(name)

    return return_list

# rate
def classifyRateData(data_dict, file_list):
    return_list = [[] for i in range(2)]
    for name in file_list:
        data = data_dict[name]
        num = int(data.params["numRetransmissions"])
        if num == 0:
            return_list[0].append(name)
        else:
            return_list[1].append(name)

    return return_list

def getSpecificInfo(data_dict, file_list, y_val):
    X = [[] for i in range(2)]
    Y = [[] for i in range(2)]
    X_labels = []
    i = 0

    for names in file_list:
        data = data_dict[names[0]]
        pre_data = data.params["numRetransmissions"]
        if pre_data == 0:
            X_labels.append("FEC")
        else:
            X_labels.append("FEC and SW-ARQ")

        for name in names:
            X[i].append(int(data_dict[name].params["FEC_rate"] * 100))
            Y[i].append(getY(data_dict[name], y_val))
        i += 1

    return [X, Y, X_labels]

def drawSpecificGraph(data_dict, all_file_list, y_val):
    dir_path = "figures"
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
    x_val = Xvalue.RATE
    file_list = classifyData(data_dict, all_file_list)

    for fileNames in file_list:
        data = data_dict[fileNames[0]]
        r = data.params["r"]
        fig_name = dir_path + '/' + "compare_distance{}_{}3.png".format(str(r), y_val.name.lower())
        files = classifyRateData(data_dict, fileNames)

        X, Y, X_labels = getSpecificInfo(data_dict, files, y_val)
        is_math = False
        if max(max(Y)) > 10 ** 5:
            is_math = True
        drawLineGraph(X, Y, X_labels, x_val, y_val, fig_name, is_math)

def drawGraph(data_dict, all_file_list):
    # Each simulation graph
    # drawBySimulation(data_dict, all_file_list)

    drawSpecificGraph(data_dict, all_file_list, Yvalue.MEDIAN)
    drawSpecificGraph(data_dict, all_file_list, Yvalue.JITTER)
    drawSpecificGraph(data_dict, all_file_list, Yvalue.FAILURE)

    # mode = checkMode(data_dict, all_file_list)
    # is_ptime, is_coll = checkOption(data_dict, all_file_list)
    # if mode is Mode.NORMAL:
    #     print("Normal")
    #     drawNormalGraph(data_dict, all_file_list, is_ptime, is_coll)
    # elif mode is Mode.DC:
    #     print("DC")
    #     drawDCGraph(data_dict, all_file_list, is_ptime, is_coll)
    # elif mode is Mode.RTO:
    #     print("RTO")
    #     drawRTOGraph(data_dict, all_file_list, is_ptime, is_coll)
    # elif mode is Mode.DECOMPOSING:
    #     print("Decomposing")
    #     drawDecomposingGraph(data_dict, all_file_list, is_ptime, is_coll)
    # elif mode is Mode.ADJUST:
    #     print("Adjust")
    #     drawAdjustGraph(data_dict, all_file_list, is_ptime, is_coll)
    #
    # elif mode is Mode.FEC:
    #     drawFECGraph(data_dict, all_file_list, is_ptime, is_coll)
