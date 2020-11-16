""" Simple MATPOWER format parser
"""
import logging
import re
import numpy as np
import Market_obj as MO


def read_structure(market, file):
    """Read a MATPOWER data file and build market elements"""
    func = re.compile(r'function\s')
    mva = re.compile(r'\s*mpc.baseMVA\s*=\s*')
    bus = re.compile(r'\s*mpc.bus\s*=\s*\[?')
    gen = re.compile(r'\s*mpc.gen\s*=\s*\[')
    branch = re.compile(r'\s*mpc.branch\s*=\s*\[')
    area = re.compile(r'\s*mpc.areas\s*=\s*\[')
    gencost = re.compile(r'\s*mpc.gencost\s*=\s*\[')
    bus_name = re.compile(r'\s*mpc.bus_name\s*=\s*{')
    end = re.compile(r'\s*\];?')
    has_digit = re.compile(r'.*\d+\s*]?;?')

    ret = True
    field = None
    info = True

    base_mva = 100
    mpc = {
        'bus': [],
        'gen': [],
        'branch': [],
        'area': [],
        'gencost': [],
        'bus_name': [],
    }

    fid = open(file, 'r')

    for line in fid:
        line = line.strip().rstrip(';')
        if not line:
            continue
        elif func.search(line):  # skip function declaration
            continue
        elif len(line.split('%')[0]) == 0:
            if info is True:
                info = False
            else:
                continue
        elif mva.search(line):
            base_mva = float(line.split('=')[1])

        if not field:
            if bus.search(line):
                field = 'bus'
            elif gen.search(line):
                field = 'gen'
            elif branch.search(line):
                field = 'branch'
            elif area.search(line):
                field = 'area'
            elif gencost.search(line):
                field = 'gencost'
            elif bus_name.search(line):
                field = 'bus_name'
            else:
                continue
        elif end.search(line):
            field = None
            continue

        # parse mpc sections
        if field:
            if line.find('=') >= 0:
                line = line.split('=')[1]
            if line.find('[') >= 0:
                line = re.sub(r'\[', '', line)
            elif line.find('{') >= 0:
                line = re.sub(r'{', '', line)

            if line.find('\'') >= 0:  # bus_name
                line = line.split(';')
                data = [i.strip('\'').strip() for i in line]
                mpc['bus_name'].extend(data)
            else:
                if not has_digit.search(line):
                    continue
                line = line.split('%')[0].strip()
                line = line.split(';')
                for item in line:
                    if not has_digit.search(item):
                        continue
                    try:
                        data = np.array([float(val) for val in item.split()])
                    except Exception as e:
                        raise e
                    mpc[field].append(data)

    fid.close()

    # convert mpc to np array
    mpc_array = dict()
    for key, val in mpc.items():
        mpc_array[key] = np.array(val)

    # list of buses with slack gen
    sw = []

    bus_idx = 0
    for data in mpc['bus']:
        # idx  ty   pd   qd  gs  bs  area  vmag  vang  baseKV  zone  vmax  vmin
        # 0    1    2   3   4   5    6      7     8     9      10    11    12
        market_type = market.Type
        load = MO.Load(market_type)
        load.P = data[2]
        bus_idx += 1
        if data[1] == 3:
            market.sw = data[0]
        market.load.append(load)
    market.Nb = int(bus_idx)

    gen_idx = 0
    for data in mpc['gen']:
        # bus  pg  qg  qmax  qmin  vg  mbase  status  pmax  pmin  pc1  pc2
        #  0   1   2    3     4     5    6      7       8    9    10    11
        # qc1min  qc1max  qc2min  qc2max  ramp_agc  ramp_10  ramp_30  ramp_q
        #  12      13       14      15      16        17       18      19
        # apf
        #  20
        gen_idx += 1
        market_type = market.Type
        gen = MO.Genco(market_type)
        gen.bus = int(data[0])
        gen.status = int(data[7])
        gen.pmax = data[8]
        gen.pmin = data[9]
        market.genco.append(gen)
    market.Ng = int(gen_idx)

    gen_idx = 0
    for data in mpc['gencost']:
        # idx	startup	shutdown  n a	b	c
        market.genco[gen_idx].bids = data[4]
        gen_idx += 1

    line_idx = 0
    for data in mpc['branch']:
        # fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle
        #  0     1      2   3   4     5       6       7       8      9
        # status	angmin	angmax	Pf	Qf	Pt	Qt
        #   10       11       12    13  14  15  16
        market_type = market.Type
        line = MO.Line(market_type)
        line.fbus = int(data[0])
        line.tbus = int(data[1])
        line.r = data[2]
        line.x = data[3]
        line.b = data[4]
        if data[5] == 0:
            line.rating = 1000000000
        else:
            line.rating = data[5]
        line.status = data[10]
        market.Line.append(line)
        line_idx += 1
    market.Nl = line_idx
    return ret


def read_load(market, file):
    fid = open(file, 'r')
    Read_Type = re.compile(r'\s*Type\s*=\s*')
    Read_load_level = re.compile(r'\s*Load_level\s*=\s*')
    Read_factor= re.compile(r'\s*Factor\s*=\s*')
    load_level = []
    factor = []
    Title = " "
    for line in fid:
        line = line.strip().rstrip(';')
        if line == "" or line == "\n":
            continue
        if Read_Type.search(line):
            Title = "Type"
        if Read_load_level.search(line):
            Title = "Load level"
            continue
        if Read_factor.search(line):
            Title = "Factor"
            continue
        # Switch "Title"  case ....
        if line.find(']') >= 0:
            Title = " "
        if Title == "Type":
            Type = line.split('=')[1]
        elif Title == "Load level":
            load_level.append(int(line))
        elif Title == "Factor":
            factor.append(int(line))
    # Read the load profile into market objective
    N_T = load_level.__len__()
    for idx, ld in enumerate(market.load):
        par_fac = factor[idx]/sum(factor)
        ld.T_P = [load_level[i]*par_fac for i in range(N_T)]
    market.N_T = N_T
    market.load_level = load_level
