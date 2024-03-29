""" Simple MATPOWER format parser
"""
import logging
import re
import numpy as np
# import numpy.random.common
# import numpy.random.bounded_integers
# import numpy.random.entropy
import market_obj as MO


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
    coordinates = re.compile(r'\s*mpc.coordinates\s*=\s*\[')
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
        'coordinates': [],
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
            elif coordinates.search(line):
                field = 'coordinates'
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

    for data in mpc['bus_name']:
        market_type = market.Type
        Bus = MO.Bus(market_type)
        Bus.name = data
        market.bus.append(Bus)

    bus_idx = 0
    for data in mpc['coordinates']:
        market.bus[bus_idx].latitude = data[0]
        market.bus[bus_idx].longitude = data[1]
        bus_idx += 1

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
        gen.ramp_up = data[18]*2
        gen.ramp_down = data[18] * 2
        market.genco.append(gen)
    market.Ng = int(gen_idx)

    gen_idx = 0
    for data in mpc['gencost']:
        # idx	startup	shutdown  ~ a	b	c
        if data[3] == 2:
            # linear
            market.genco[gen_idx].bid_type = 2
            market.genco[gen_idx].bids = data[4]
        elif data[3] == 3:
            # quadratic
            market.genco[gen_idx].bid_type = 3
            market.genco[gen_idx].bids = [data[4], data[5], data[6],]
        market.genco[gen_idx].start_up = data[1]
        market.genco[gen_idx].shut_down = data[2]
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
    Read_load = re.compile(r'\s*Load\s*=\s*')
    Title = " "
    for line in fid:
        line = line.strip().rstrip(';')
        if Read_Type.search(line):
            Type = int(line.split('=')[1])
    fid.close()
    load = []
    load_level = []
    factor = []
    if Type == 1:
        fid = open(file, 'r')
        for line in fid:
            if Read_load.search(line):
                Title = "Load"
                continue
            if line.find(']') >= 0:
                Title = " "
            elif Title == "Load":
                load_txt = line.split()
                load_float = [float(load_txt[i]) for i in range(load_txt.__len__())]
                load.append(load_float)
        N_T = load[0].__len__()
        for idx, ld in enumerate(market.load):
            ld.T_P = load[idx]
        for i in range(N_T):
            load_level.append(sum([load[x][i] for x in range(market.load.__len__())]))
        market.N_T = N_T
        market.load_level = load_level
    if Type == 2:
        fid = open(file, 'r')
        for line in fid:
            line = line.strip().rstrip(';')
            if line == "" or line == "\n":
                continue
            if Read_load_level.search(line):
                Title = "Load level"
                continue
            if Read_factor.search(line):
                Title = "Factor"
                continue
            # Switch "Title"  case ....
            if line.find(']') >= 0:
                Title = " "
            elif Title == "Load level":
                load_level.append(float(line))
            elif Title == "Factor":
                factor.append(float(line))
        # Read the load profile into market objective
        N_T = load_level.__len__()
        for idx, ld in enumerate(market.load):
            par_fac = factor[idx]/sum(factor)
            ld.T_P = [load_level[i]*par_fac for i in range(N_T)]
        market.N_T = N_T
        market.load_level = load_level


def read_xgd(market, file):
    fid = open(file, 'r')
    Read_factor= re.compile(r'\s*xgd_table.data\s*=\s*')
    CommitKey = []
    MinUp = []
    MinDown = []
    Title = " "
    for line in fid:
        line = line.strip().rstrip(';')
        if line == "" or line == "\n":
            continue
        if Read_factor.search(line):
            Title = "xgd_table.data"
            continue
        # Switch "Title"  case ....
        if line.find(']') >= 0:
            Title = " "
        if Title == "xgd_table.data":
            line = line.split(' ')
            line = [x for x in line if x]
            CommitKey.append(int(line[0]))
            MinUp.append(int(line[2]))
            MinDown.append(int(line[3]))
    # Read the load profile into market objective
    for idx, gen in enumerate(market.genco):
        gen.commit_key = CommitKey[idx]
        gen.min_dowm = MinDown[idx]
        gen.min_up = MinUp[idx]
