import re
from enum import Enum
import linecache
# import sys
import os
import numpy as np
from statistics import mean, variance, stdev, median
import math
from natsort import natsorted

class SimData:
    def __init__(self, config_file_name):
        self.config_file_name = config_file_name
        self.parseConfig()
        self.input_data = InputData(self)
        self.output_data = OutputData(self)
        self.plot_data = PlotData(self)
        self.regData()
        self.defineFileName()

        self.is_ptime = False
        self.is_coll = True

    def setMeanGraph(self, fig_name):
        self.mean_fig_name = fig_name

    def setProbCumprobGraph(self, fig_name):
        self.prob_cumprob_fig_name = fig_name

    def defineFileName(self):
        # self.file_name = "TxRx{}_ARQ{}-{}_{}-{}".format(self.params["r"], self.params["info_arq"], self.params["ack_arq"], self.params["info_type"], self.params["ack_type"])
        self.file_name = "TxRx{}_ARQ{}-{}_{}-{}_RTO{}_DC{}_FEC{}-{}".format(self.params["r"], self.params["info_arq"], self.params["ack_arq"], self.params["info_type"], self.params["ack_type"], self.params["rto_type"], str(self.params["step_length"]).replace('.', ''), self.params["FEC_packet"], str(self.params["FEC_rate"]).replace('.', ''))
        self.collision_file_name = "./result/collision_batch_" + self.config["outputFile"]
        self.decomposing_file_name = "./result/decomposing_batch_" + self.config["outputFile"]
        self.adjust_file_name = "./result/adjust_batch_" + self.config["outputFile"]
        self.retransmission_file_name = "./result/retransmission_batch_" + self.config["outputFile"]

        # if "rto_type" in self.params.keys():
        #     self.file_name += "_RTO{}".format(self.params["rto_type"])
        # if self.params["decomposing"]:
        #     self.file_name += "_decomposing"
        #
        # if self.params["adjust"]:
        #     self.file_name += "_adjust{}".format(self.params["adjust"])
        #
        # self.file_name += "_{}message".format(self.params["message_num"])

    def regData(self):
        self.params = {}
        # d, stepLength, rate
        self.params["r"] = int(self.output_data.analytical_model.r)
        for mp in self.config["moleculeParams"]:
            if mp.type_of_molecule == MoleculeType["INFO"]:
                self.params["info_arq"] = mp.duplication
                self.params["duplication"] = mp.duplication
                self.params["info_type"] = mp.type_of_movement.name
            elif mp.type_of_molecule == MoleculeType["ACK"]:
                self.params["ack_arq"] = mp.duplication
                self.params["ack_type"] = mp.type_of_movement.name
        self.params["step_length"] = self.config["stepLengthX"]
        if "RTO" in self.config_file_name:
            self.params["rto_type"] = self.config_file_name.split("RTO")[1][0]
        if self.config["FEC"] is not None:
            self.params["FEC_rate"] = self.config["FEC"].rate
            self.params["FEC_packet"] = self.config["FEC"].require_packet
        self.params["mean"] = self.output_data.mean
        self.params["median"] = self.output_data.med
        self.params["std"] = self.output_data.std
        self.params["failure"] = self.output_data.transmission_failure_rate

        # self.params["r"] = int(self.output_data.analytical_model.r)
        # for mp in self.config["moleculeParams"]:
        #     if mp.type_of_molecule == MoleculeType["INFO"]:
        #         self.params["info_arq"] = mp.duplication
        #         self.params["info_type"] = mp.type_of_movement.name
        #         self.params["info_adjust"] = mp.adaptive_change_number
        #     elif mp.type_of_molecule == MoleculeType["ACK"]:
        #         self.params["ack_arq"] = mp.duplication
        #         self.params["ack_type"] = mp.type_of_movement.name
        #         self.params["ack_adjust"] = mp.adaptive_change_number
        # self.params["step_length"] = self.config["stepLengthX"]
        # self.params["decomposing"] = self.config["decomposing"]
        # self.params["duplication"] = self.params["info_arq"]
        # self.params["mean"] = self.output_data.mean
        # self.params["analytical_rtt"] = self.output_data.analytical_model.rtt
        # self.params["median"] = self.output_data.med
        # self.params["std"] = self.output_data.std
        # self.params["retransmitWaitTime"] = self.config["retransmitWaitTime"]
        # self.params["txrx_mean"] = self.output_data.txrx_mean
        # self.params["failure"] = self.output_data.transmission_failure_rate
        # self.params["adjust"] = self.params["info_adjust"]
        # # self.params["decomposing_num"] = self.input_data.decomposing_sum
        # self.params["message_num"] = self.config["numMessages"]
        # if "RTO" in self.config_file_name:
        #     self.params["rto_type"] = self.config_file_name.split("RTO")[1][0]
        #
        # self.is_ptime = self.input_data.is_ptime
        # self.is_coll = self.input_data.is_coll
        #
        # if self.config["FEC"] is not None:
        #     self.params["FEC_rate"] = self.config["FEC"].rate
        #     self.params["FEC_packet"] = self.config["FEC"].require_packet
        #
        # if self.is_coll:
        #     self.params["coll"] = self.input_data.coll_avg
        #
        # if os.path.isfile("./result/adjust_batch_" + self.config["outputFile"]):
        #     self.params["last_info_num"] = self.input_data.last_info_num
        #     self.params["last_ack_num"] = self.input_data.last_ack_num
        #
        # if os.path.isfile("./result/retransmission_batch_" + self.config["outputFile"]):
        #     self.params["retransmit"] = self.input_data.retransmit_steps
        #     self.params["retransmit_failure"] = self.input_data.failure_count


    def parseConfig(self):
        self.config = {}
        self.config["moleculeParams"] = []
        self.config["microtubuleParams"] = []
        self.config_header = []

        with open(self.config_file_name, 'r') as f:
            for line in f:
                if line[0] == '*' or line[0] == '\n':
                    continue
                line = line.rstrip()
                key, val = line.split(' ', 1)
                if not key in self.config_header:
                    self.config_header.append(key)

                if key in ["transmitter", "receiver"]:
                    self.config[key] = NanoMachine(val)
                elif key in ["intermediateNode"]:
                    self.config[key] = IntermediateNode(val)
                elif key in ["moleculeParams"]:
                    self.config[key].append(MoleculeParams(val))
                elif key in ["microtubuleParams"]:
                    self.config[key].append(MicrotubuleParams(val))
                elif key in ["probDRail", "stepLengthX", "stepLengthY", "stepLengthZ", "packetStepLengthX", "packetStepLengthY", "packetStepLengthZ", "packetDiameter"]:
                    self.config[key] = float(val)
                elif key in ["outputFile"]:
                    self.config[key] = val
                elif key in ["FEC"]:
                    self.config[key] = FEC(val)
                else:
                    self.config[key] = int(val)

class InputData:
    def __init__(self, data):
        self.input_file_name = "./result/batch_" + data.config["outputFile"]
        self.parseFile()
        # self.is_ptime = True
        self.is_ptime = False
        self.is_coll = True
        # self.txrx_mean = 0

        # self.parseCollisionFile(data)
        # self.parseAdjustFile(data)
        self.parseRetransmissionFile(data)
        # self.parseWaitFile(data)

        # if self.is_ptime:
        #     self.calcAvg()
        # if self.is_coll:
        #     self.calcColl()

    def parseWaitFile(self, data):
        file_name = "./result/wait_batch_" + data.config["outputFile"]
        if not os.path.isfile(file_name):
            self.is_ptime = False
            return

        self.info_time = []
        self.info_num = []
        self.ack_time = []
        self.ack_num = []

        with open(file_name, 'r') as f:
            for line in f:
                datas = line.split(',')
                self.info_time.append(datas[0])
                self.info_num.append(datas[1])
                self.ack_time.append(datas[2])
                self.ack_num.append(datas[3])

    def parseCollisionFile(self, data):
        file_name = "./result/collision_batch_" + data.config["outputFile"]
        if not os.path.isfile(file_name):
            self.is_coll = False
            return

        self.collision_header = ["ACK/ACK", "ACK/INFO", "ACK/NOISE", "INFO/INFO", "INFO/NOISE"]
        self.collisions = []
        self.each_collision = []
        self.decomposing_num = []

        with open(file_name, 'r') as f:
            for line in f:
                if len(line.split(',')) == 3:
                    all_collision, collision_num, decomposing = line.split(',')
                    self.collisions.extend(all_collision.split('/'))
                    self.each_collision.append(collision_num.split('/'))
                    self.decomposing_num.append(decomposing)
                else:
                    all_collision, collision_num = line.split(',')
                    self.collisions.extend(all_collision.split('/'))
                    self.each_collision.append(collision_num.split('/'))

        self.collisions = natsorted(self.collisions)
        self.decomposing_sum = sum(list(map(int, self.decomposing_num)))

    def parseAdjustFile(self, data):
        file_name = "./result/adjust_batch_" + data.config["outputFile"]
        if not os.path.isfile(file_name):
            return

        self.adjust_steps = []
        self.info_adjust_num = []
        self.ack_adjust_num = []

        with open(file_name, 'r') as f:
            for line in f:
                datas = line.split(',')
                for data in datas:
                    step, info, ack = data.split('/')
                    self.adjust_steps.append(step)
                    self.info_adjust_num.append(info)
                    self.ack_adjust_num.append(ack)

        self.last_info_num = self.info_adjust_num[-1]
        self.last_ack_num = self.ack_adjust_num[-1]

    def parseRetransmissionFile(self, data):
        file_name = "./result/retransmission_batch_" + data.config["outputFile"]
        if not os.path.isfile(file_name):
            return

        length = len(linecache.getline(file_name, 1).split(','))
        self.retransmit_steps = []
        self.tx_retransmit_steps = []
        self.rx_retransmit_steps = []
        self.failure_count = 0

        with open(file_name, 'r') as f:
            for line in f:
                datas = line.split(',')
                if datas[0] == "F":
                    self.failure_count += 1

                self.retransmit_steps.extend(datas[1].split('/'))
                if len(datas) != 2:
                    self.tx_retransmit_steps.extend(datas[2].split('/'))
                    self.rx_retransmit_steps.extend(datas[3].split('/'))

        self.retransmit_steps = natsorted(self.retransmit_steps)
        self.failure_count = self.failure_count / len(self.steps)
        if length != 2:
            self.tx_retransmit_steps = natsorted(self.tx_retransmit_steps)
            self.rx_retransmit_steps = natsorted(self.rx_retransmit_steps)

    def calcAvg(self):
        info_sum_time = float(sum(self.info_time))
        info_sum_num = sum(self.info_num)
        ack_sum_time = float(sum(self.ack_time))
        ack_sum_num = sum(self.ack_num)
        self.txrx_mean = info_sum_time / info_sum_num + ack_sum_time / ack_sum_num

    def calcColl(self):
        self.coll_aa = []
        self.coll_ai = []
        self.coll_an = []
        self.coll_ii = []
        self.coll_in = []

        for collisions in self.each_collision:
            self.coll_aa.append(int(collisions[0]))
            self.coll_ai.append(int(collisions[1]))
            self.coll_an.append(int(collisions[2]))
            self.coll_ii.append(int(collisions[3]))
            self.coll_in.append(int(collisions[4]))

        self.sum_coll_aa = sum(self.coll_aa)
        self.sum_coll_ai = sum(self.coll_ai)
        self.sum_coll_an = sum(self.coll_an)
        self.sum_coll_ii = sum(self.coll_ii)
        self.sum_coll_in = sum(self.coll_in)
        self.sum_coll = self.sum_coll_aa + self.sum_coll_ai + self.sum_coll_an + self.sum_coll_ii + self.sum_coll_in
        self.coll_avg = self.sum_coll / len(self.steps)

    def parseFile(self):
        self.steps = []
        # length = len(linecache.getline(self.input_file_name, 1).split(','))
        # if length != 1:
        #     self.oldParseFile(length)
        #     return

        with open(self.input_file_name, 'r') as f:
            for line in f:
                self.steps.append(int(line))

    # def oldParseFile(self, length):
    #     self.steps = []
    #     self.coll_aa = []
    #     self.coll_ai = []
    #     self.coll_an = []
    #     self.coll_ii = []
    #     self.coll_in = []
    #     self.is_ptime = False
    #     with open(self.input_file_name, 'r') as f:
    #         if length == 5:
    #             for line in f:
    #                 data = line.split(',')
    #                 self.steps.append(int(data[0]))
    #                 # self.info_time.append(int(data[1]))
    #                 # self.info_num.append(int(data[2]))
    #                 # self.ack_time.append(int(data[3]))
    #                 # self.ack_num.append(int(data[4]))
    #                 # self.is_ptime = True
    #                 self.is_coll = False
    #
    #         if length == 6:
    #             for line in f:
    #                 line.rstrip('\n')
    #                 data = line.split(',')
    #                 self.steps.append(int(data[0]))
    #                 self.coll_aa.append(int(data[1]))
    #                 self.coll_ai.append(int(data[2]))
    #                 self.coll_an.append(int(data[3]))
    #                 self.coll_ii.append(int(data[4]))
    #                 self.coll_in.append(int(data[5]))
    #
    #             self.sum_coll_aa = sum(self.coll_aa)
    #             self.sum_coll_ai = sum(self.coll_ai)
    #             self.sum_coll_an = sum(self.coll_an)
    #             self.sum_coll_ii = sum(self.coll_ii)
    #             self.sum_coll_in = sum(self.coll_in)
    #             self.sum_coll = self.sum_coll_aa + self.sum_coll_ai + self.sum_coll_an + self.sum_coll_ii + self.sum_coll_in
    #             self.coll_avg = self.sum_coll / len(self.steps)
    #
    #         if length == 1:
    #             for line in f:
    #                 self.steps.append(int(line))
    #
    #             self.sum_coll_aa = 0
    #             self.sum_coll_ai = 0
    #             self.sum_coll_an = 0
    #             self.sum_coll_ii = 0
    #             self.sum_coll_in = 0
    #             self.sum_coll = self.sum_coll_aa + self.sum_coll_ai + self.sum_coll_an + self.sum_coll_ii + self.sum_coll_in
    #             self.coll_avg = self.sum_coll / len(self.steps)

            # elif self.is_ptime and self.is_coll:
            #     for line in f:
            #         data = line.split(',')
            #         self.steps.append(int(data[0]))
            #         self.info_time.append(int(data[1]))
            #         self.info_num.append(int(data[2]))
            #         self.ack_time.append(int(data[3]))
            #         self.ack_num.append(int(data[4]))
            #         self.coll_aa.append(int(data[5]))
            #         self.coll_ai.append(int(data[6]))
            #         self.coll_an.append(int(data[7]))
            #         self.coll_ii.append(int(data[8]))
            #         self.coll_in.append(int(data[9]))

class OutputData:
    def __init__(self, data):
        self.steps = np.sort(data.input_data.steps)
        self.count = len(self.steps)
        self.minimum = self.steps[0]
        self.maximum = self.steps[-1]
        self.range_num = int(self.maximum / 1000) * 10
        if self.range_num == 0:
            self.range_num = int(self.maximum / 100)
        self.var = np.var(self.steps)
        self.std = np.std(self.steps)
        self.med = median(self.steps)
        self.mean = mean(self.steps)
        # self.txrx_mean = data.input_data.txrx_mean
        self.analytical_model = AnalyticalModel(data.config)
        # self.threashold = self.med
        if self.analytical_model.r == 30:
            self.threashold = (17028 + self.med) / 2
        elif self.analytical_model.r == 50:
            self.threashold = (21837 + self.med) / 2
        elif self.analytical_model.r == 70:
            self.threashold = (22915 + self.med) / 2
        elif self.analytical_model.r == 90:
            self.threashold = (21160 + self.med) / 2

        self.calcTransmissionFailureRate()
        self.is_ptime = data.input_data.is_ptime

    def calcTransmissionFailureRate(self):
        length = len(self.steps)
        o = 0
        for step in self.steps:
            if step > self.threashold:
                break
            o += 1
        x = length - o
        self.transmission_failure_rate = x / length * 100

    def toArray(self):
        if self.is_ptime:
            return [self.count, self.minimum, self.maximum, self.range_num,
            self.var, self.std, self.med, self.mean, self.txrx_mean]
        else:
            return [self.count, self.minimum, self.maximum, self.range_num,
            self.var, self.std, self.med, self.mean]

class AnalyticalModel:
    def __init__(self, config):
        tx = config["transmitter"].center_position
        rx = config["receiver"].center_position
        step_length = config["stepLengthX"]
        if step_length in DNA_DICT.keys():
            self.D = DNA_DICT[step_length]["dc"]
        else:
            self.D = 0.5
        self.r = self.calcDistance(tx, rx)
        self.L = int(config["mediumDimensionX"] / 2)
        self.l = (int(config["receiver"].size) * 2 - 1) / 2

        self.info_winf = 0.0
        self.ack_winf = 0.0

        molecule_params = config["moleculeParams"]
        for molecule_param in molecule_params:
            # from transmitter to receiver
            if molecule_param.type_of_molecule == MoleculeType["INFO"]:
                if molecule_param.type_of_movement == MovementType["PASSIVE"]:
                    self.info_winf = self.calcPassiveRtt()
                else:
                    passive_winf = self.calcPassiveRtt()
                    self.info_winf = self.calcActiveRtt(passive_winf)

            # from receiver to transmitter
            elif molecule_param.type_of_molecule == MoleculeType["ACK"]:
                if molecule_param.type_of_movement == MovementType["PASSIVE"]:
                    self.ack_winf = self.calcPassiveRtt()
                else:
                    passive_winf = self.calcPassiveRtt()
                    self.ack_winf = self.calcActiveRtt(passive_winf)

        self.rtt = self.info_winf + self.ack_winf

    def calcPassiveRtt(self):
        winf = (self.r - self.l) * (2 * self.L ** 3 - self.l * self.r ** 2 - self.l ** 2 * self.r) / (2 * self.D * self.l * self.r)
        return winf

    def calcActiveRtt(self, passive_winf):
        V1 = 2 * 2 * self.r
        V = 4.0 / 3.0 * math.pi * self.L ** 3
        p = V1 / V
        va = 1
        vp = self.r / passive_winf
        ve = p * va + (1.0 - p) * vp
        winf = self.r / ve
        return winf

    def calcDistance(self, tx, rx):
        x = (tx.x - rx.x) ** 2
        y = (tx.y - rx.y) ** 2
        z = (tx.z - rx.z) ** 2
        return math.sqrt(x + y + z)

    def toArray(self):
        return [self.D, self.r, self.L, self.l, self.rtt]

class PlotData:
    def __init__(self, data):
        self.range_num = data.output_data.range_num

        self.makePlotData(data.output_data.steps)
        self.calcProb(len(data.output_data.steps))

    def calcProb(self, step_size):
        self.prob = []
        self.cum_prob = []
        prob = 0.0
        cum_prob = 0.0
        for i in range(len(self.plot_data)):
            head, tail = self.plot_range[i]
            prob = float(len(self.plot_data[i])) / step_size * 100.0
            cum_prob += prob
            self.prob.append(prob)
            self.cum_prob.append(cum_prob)

    def makePlotData(self, datas):
        self.plot_data = []
        self.plot_range = []
        data_arr = []
        head = 0
        tail = self.range_num
        self.plot_range.append([head, tail - 1])
        i = 0
        while True:
            if datas[i] < tail:
                data_arr.append(datas[i])
                i += 1
                if i == len(datas):
                    break
            else:
                self.plot_data.append(data_arr)
                head += self.range_num
                tail += self.range_num
                self.plot_range.append([head, tail - 1])
                data_arr = []

        self.plot_data.append(data_arr)

class Position:
    def __init__(self, args):
        self.x = int(args[0])
        self.y = int(args[1])
        self.z = int(args[2])

    def toString(self):
        return "({}, {}, {})".format(self.x, self.y, self.z)

class MicrotubuleParams:
    def __init__(self, val):
        args = [i for i in re.split(r"[,( )]", val) if i != '']
        self.start_position = Position(args[0:3])
        self.end_position = Position(args[3:6])

    def toString(self):
        return "{} {}".format(self.start_position.toString(), self.end_position.toString())

class MoleculeParams:
    def __init__(self, val):
        args = [i for i in re.split(r"[ ]", val) if i != '']
        self.duplication = int(args[0])
        self.type_of_molecule = MoleculeType[args[1]]
        self.size = 1
        if self.type_of_molecule != MoleculeType["NOISE"]:
            self.type_of_movement = MovementType[args[2]]
            self.adaptive_change_number = int(args[3])
            if len(args) == 5:
                self.size = float(args[4])
            else:
                self.size = float(1)
        else:
            if len(args) == 3:
                self.size = float(args[2])
            else:
                self.size = float(1)

    def toString(self):
        if self.type_of_molecule != MoleculeType["NOISE"]:
            return "{} {} {} {} {}".format(self.duplication, self.type_of_molecule.name, self.type_of_movement.name, self.adaptive_change_number, self.size)
        else:
            return "{} {} {}".format(self.duplication, self.type_of_molecule.name, self.size)

class IntermediateNode:
    def __init__(self, val):
        args = [i for i in re.split(r"[,( )]", val) if i != '']
        self.center_position = Position(args[0:3])
        self.size = int(args[3])
        self.info_release_position = Position(args[4:7])
        self.ack_release_position = Position(args[7:10])

    def toString(self):
        return "{} {} {} {}".format(self.center_position.toString(), self.size, self.info_release_position.toString(), self.ack_release_position.toString())

class NanoMachine:
    def __init__(self, val):
        args = [i for i in re.split(r"[,( )]", val) if i != '']
        self.center_position = Position(args[0:3])
        self.size = int(args[3])
        self.release_position = Position(args[4:7])

    def toString(self):
        return "{} {} {}".format(self.center_position.toString(), self.size, self.release_position.toString())

class FEC:
    def __init__(self, val):
        args = [i for i in re.split(r"[,( )]", val) if i != '']
        self.type = args[0]
        self.require_packet = int(args[1])
        self.rate = float(args[2])

class MoleculeType(Enum):
    INFO = 0
    ACK = 1
    NOISE = 2

class MovementType(Enum):
    PASSIVE = 0
    ACTIVE = 1

DNA_DICT = {
    11.31: {"diameter": 0.01, "dc": 63.944},
    9.19: {"diameter": 0.02, "dc": 42.187},
    7.46: {"diameter": 0.02, "dc": 27.833},
    6.06: {"diameter": 0.04, "dc": 18.363},
    4.92: {"diameter": 0.06, "dc": 12.115},
    4.00: {"diameter": 0.09, "dc": 7.993},
    3.25: {"diameter": 0.13, "dc": 5.273},
    2.64: {"diameter": 0.20, "dc": 3.479},
    2.14: {"diameter": 0.30, "dc": 2.295},
    1.18: {"diameter": 0.98, "dc": 0.695},
    0.67: {"diameter": 3.00, "dc": 0.227},
    2.20: {"diameter": 0.28, "dc": 2.416},
    1.42: {"diameter": 0.67, "dc": 1.013},
    10.58: {"diameter": 0.01, "dc": 55.931},
}
