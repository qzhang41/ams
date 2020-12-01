import numpy as np
import gurobipy as gb


def Make_Bdc(Market):
    Nb = Market.Nb
    Nl = Market.Nl
    Cft = np.zeros([Nl, Nb])
    Bf = np.zeros([Nl, Nb])
    Bbus = np.zeros([Nb, Nb])

    Lines = Market.Line
    for idx, l in enumerate(Lines):
        status = l.status
        Bf[idx, int(l.fbus - 1)] = status / l.x
        Bf[idx, int(l.tbus - 1)] = -status / l.x
        Cft[idx, int(l.fbus - 1)] = 1
        Cft[idx, int(l.tbus - 1)] = -1
    Bbus = np.matmul(np.transpose(Cft), Bf)
    Market.Bbus = Bbus
    Market.Bf = Bf
    Market.Cft = Cft


def Make_PTDF(Market):
    Nb = Market.Nb
    Nl = Market.Nl
    sw = Market.sw
    nosw = range(Nb)
    nosw = np.setdiff1d(nosw, sw - 1)
    Bf = Market.Bf
    Bbus = Market.Bbus
    PTDF = np.matmul(Bf[:, nosw], np.linalg.inv(Bbus[np.ix_(nosw, nosw)]))
    PTDF = np.insert(PTDF, int(sw - 1), 0, axis=1)
    Market.PTDF = PTDF


def Ecnomic_dispatch(Market):
    Make_Bdc(Market)
    Make_PTDF(Market)
    opt_model = gb.Model(str(Market.Type) + 'Ecnomic_dispatch')
    Pg = {}
    load_level = 0
    obj = 0
    gen_bus = np.zeros([len(Market.genco), 1])
    # add Pg cap
    for idx, gen in enumerate(Market.genco):
        gen_bus[idx] = gen.bus
        Pg[idx] = opt_model.addVar(name='Power generation' + str(idx), vtype=gb.GRB.CONTINUOUS, \
                                   ub=gen.pmax, lb=gen.pmin)
        cost = gen.bids
        obj += Pg[idx]*cost
    # add line flow cons
    line_flow = {}
    for line_idx, line in enumerate(Market.Line):
        line_flow[line_idx] = 0
        for bus_idx in range(Market.Nb):
            load = Market.load[bus_idx].P
            line_flow[line_idx] = line_flow[line_idx] + Market.PTDF[line_idx, bus_idx] * (-load)
            line_flow[line_idx] = line_flow[line_idx] + Market.PTDF[line_idx, bus_idx] * sum([Pg[x] for x in sum(np.where(gen_bus == bus_idx + 1))])
        opt_model.addConstr(line_flow[line_idx] <= line.rating, name='TC p' + str(line_idx))
        opt_model.addConstr(line_flow[line_idx] >= -line.rating, name='TC n' + str(line_idx))
    # add power balance
    load_level = sum([np.sum(Market.load[i].P) for i in range(Market.load.__len__())])
    gen_level = sum([np.sum(Pg[i]) for i in range(Market.Ng)])
    opt_model.addConstr(gen_level == load_level, name="balance")
    opt_model.setObjective(obj, gb.GRB.MINIMIZE)
    opt_model.optimize()
    #opt_model.write('math_model.lp')
    for idx, gen in enumerate(Market.genco):
        gen.opt_pg = Pg[idx].X
    for idx, line in enumerate(Market.Line):
        line.opt_fl = line_flow[idx].getValue()
    # LMP and dispatched settlements
    lamda = opt_model.getConstrByName('balance').Pi
    LMP = np.zeros([1, Market.Nb])
    for b, ld in enumerate(Market.load):
        LMP[0, b] = lamda
        for l, line in enumerate(Market.Line):
            ng = opt_model.getConstrByName('TC n' + str(l)).Pi
            po = opt_model.getConstrByName('TC p' + str(l)).Pi
            LMP[0, b] += Market.PTDF[l, b]*(ng-po)
    Market.LMP = LMP
    for gen in Market.genco:
        gen.Revenue = Market.LMP[0, int(gen.bus-1)]*gen.opt_pg
    for idx, ld in enumerate(Market.load):
        ld.Revenue = -Market.LMP[0, idx]*ld.P


def Multi_Ecnomic_dispatch(Market):
    Make_Bdc(Market)
    Make_PTDF(Market)
    for t in range(Market.N_T):
        opt_model = gb.Model(str(Market.Type) + 'Ecnomic_dispatch')
        Pg = {}
        load_level = 0
        obj = 0
        gen_bus = np.zeros([len(Market.genco), 1])
        # add Pg cap
        for idx, gen in enumerate(Market.genco):
            gen_bus[idx] = gen.bus
            Pg[idx] = opt_model.addVar(name='Power generation' + str(idx), vtype=gb.GRB.CONTINUOUS, \
                                       ub=gen.pmax, lb=gen.pmin)
            cost = gen.bids
            obj += Pg[idx]*cost
        # add line flow cons
        line_flow = {}
        for line_idx, line in enumerate(Market.Line):
            line_flow[line_idx] = 0
            for bus_idx in range(Market.Nb):
                load = Market.load[bus_idx].T_P[t]
                line_flow[line_idx] = line_flow[line_idx] + Market.PTDF[line_idx, bus_idx] * (-load)
                line_flow[line_idx] = line_flow[line_idx] + Market.PTDF[line_idx, bus_idx] *\
                                      sum([Pg[x] for x in sum(np.where(gen_bus == bus_idx + 1))])
            opt_model.addConstr(line_flow[line_idx] <= line.rating, name='TC p' + str(line_idx))
            opt_model.addConstr(line_flow[line_idx] >= -line.rating, name='TC n' + str(line_idx))
        # add power balance
        load_level = Market.load_level
        gen_level = sum([np.sum(Pg[i]) for i in range(Market.Ng)])
        opt_model.addConstr(gen_level == load_level[t], name="balance")
        opt_model.setObjective(obj, gb.GRB.MINIMIZE)
        opt_model.optimize()
        if opt_model.Status == 2:
            for idx, gen in enumerate(Market.genco):
                gen.T_Pg.append(Pg[idx].X)
            for idx, line in enumerate(Market.Line):
                line.T_Lf.append(line_flow[idx].getValue())
            # LMP and dispatched settlements
            lamda = opt_model.getConstrByName('balance').Pi
            LMP = np.zeros([1, Market.Nb])
            for b, ld in enumerate(Market.load):
                LMP[0, b] = lamda
                for l, line in enumerate(Market.Line):
                    ng = opt_model.getConstrByName('TC n' + str(l)).Pi
                    po = opt_model.getConstrByName('TC p' + str(l)).Pi
                    LMP[0, b] += Market.PTDF[l, b]*(ng-po)
            Market.T_LMP.append(LMP)
            for gen in Market.genco:
                gen.T_Revenue.append(LMP[0, int(gen.bus-1)]*gen.T_Pg[t])
            for idx, ld in enumerate(Market.load):
                ld.Revenue.append(-LMP[0, idx]*ld.P)
        else:
            for idx, gen in enumerate(Market.genco):
                gen.T_Pg.append([])
                gen.T_Revenue.append([])
            for idx, line in enumerate(Market.Line):
                line.T_Lf.append([])
            for idx, ld in enumerate(Market.load):
                ld.Revenue.append([])
            Market.T_LMP.append([])
        del opt_model


def Unit_commitment(Market):
    Make_Bdc(Market)
    Make_PTDF(Market)
    opt_model = gb.Model(str(Market.Type) + 'Unit Commitment')
    Pg = {}
    Status = {}
    Comm = {}
    obj = 0
    for t in range(Market.N_T):
        gen_bus = np.zeros([len(Market.genco), 1])
        # add Pg cap and obj
        for idx, gen in enumerate(Market.genco):
            gen_bus[idx] = gen.bus
            Pg[t, idx] = opt_model.addVar(name='Power generation' + str(idx), vtype=gb.GRB.CONTINUOUS)
            Status[t, idx] = opt_model.addVar(name='T:' + str(t) + 'Status' + str(idx), vtype=gb.GRB.BINARY)
            Comm[t, idx] = opt_model.addVar(name='T:' + str(t) + 'Comm' + str(idx), vtype=gb.GRB.BINARY)
            # start up variable
            if t == 0:
                opt_model.addConstr(Status[t, idx] - gen.status <= Comm[t, idx])
                opt_model.addConstr(Status[t, idx] >= Comm[t, idx])
            else:
                opt_model.addConstr(Status[t, idx] - Status[t-1, idx] <= Comm[t, idx])
                opt_model.addConstr(Status[t, idx] >= Comm[t, idx])
            # ramp
            if t >= 1:
                opt_model.addConstr(Pg[t, idx] - Pg[t-1, idx] <= 0,      #gen.T_ramp_up[t],
                                    name='T:' + str(t) + 'ramp_up' + str(idx))
                opt_model.addConstr(Pg[t, idx] - Pg[t-1, idx] >= 0,        #gen.T_ramp_down[t],
                                    name='T:' + str(t) + 'ramp_up' + str(idx))
            # power capacity
            opt_model.addConstr(Pg[t, idx] <= gen.pmax*Status[t, idx], name='T:' + str(t) + 'Capacity_max' + str(idx))
            opt_model.addConstr(Pg[t, idx] >= gen.pmin*Status[t, idx], name='T:' + str(t) + 'Capacity_min' + str(idx))
            # obj
            cost = gen.bids
            no_load_c = 0      #gen.T_no_load[t]
            start_up_c = 0     #gen.T_no_load[t]
            obj += Pg[t, idx]*cost+Status[t, idx]*no_load_c+Comm[t, idx]*start_up_c
        # add line flow cons
        line_flow = {}
        for line_idx, line in enumerate(Market.Line):
            line_flow[t, line_idx] = 0
            for bus_idx in range(Market.Nb):
                load = Market.load[bus_idx].T_P[t]
                line_flow[t, line_idx] = line_flow[t, line_idx] + Market.PTDF[line_idx, bus_idx] * (-load)
                line_flow[t, line_idx] = line_flow[t, line_idx] + Market.PTDF[line_idx, bus_idx] *\
                                      sum([Pg[t, x] for x in sum(np.where(gen_bus == bus_idx + 1))])
            opt_model.addConstr(line_flow[t, line_idx] <= line.rating, name='T:' + str(t) + 'TC p' + str(line_idx))
            opt_model.addConstr(line_flow[t, line_idx] >= -line.rating, name='T:' + str(t) + 'TC n' + str(line_idx))
        # add power balance
        load_level = Market.load_level
        gen_level = sum([np.sum(Pg[t, i]) for i in range(Market.Ng)])
        opt_model.addConstr(gen_level == load_level[t], name="balance")
    opt_model.setObjective(obj, gb.GRB.MINIMIZE)
    opt_model.optimize()

