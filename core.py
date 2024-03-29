import numpy as np
import gurobipy as gb
import Disturb as db


def make_Bdc(market):
    Nb = market.Nb
    Nl = market.Nl
    Cft = np.zeros([Nl, Nb])
    Bf = np.zeros([Nl, Nb])
    Bbus = np.zeros([Nb, Nb])

    Lines = market.Line
    for idx, l in enumerate(Lines):
        status = l.status
        Bf[idx, int(l.fbus - 1)] = status / l.x
        Bf[idx, int(l.tbus - 1)] = -status / l.x
        Cft[idx, int(l.fbus - 1)] = 1
        Cft[idx, int(l.tbus - 1)] = -1
    Bbus = np.matmul(np.transpose(Cft), Bf)
    market.Bbus = Bbus
    market.Bf = Bf
    market.Cft = Cft


def make_PTDF(market):
    Nb = market.Nb
    Nl = market.Nl
    sw = market.sw
    nosw = range(Nb)
    nosw = np.setdiff1d(nosw, sw - 1)
    Bf = market.Bf
    Bbus = market.Bbus
    PTDF = np.matmul(Bf[:, nosw], np.linalg.inv(Bbus[np.ix_(nosw, nosw)]))
    PTDF = np.insert(PTDF, int(sw - 1), 0, axis=1)
    market.PTDF = PTDF


def ecnomic_dispatch(market):
    make_Bdc(market)
    make_PTDF(market)
    opt_model = gb.Model(str(market.Type) + 'Ecnomic_dispatch')
    pg = {}
    load_level = 0
    obj = 0
    gen_bus = np.zeros([len(market.genco), 1])
    # add pg cap
    for idx, gen in enumerate(market.genco):
        gen_bus[idx] = gen.bus
        pg[idx] = opt_model.addVar(name='Power generation' + str(idx), vtype=gb.GRB.CONTINUOUS,
                                   ub=gen.pmax * gen.status, lb=gen.pmin * gen.status)
        if gen.bid_type == 2:
            cost = gen.bids
            obj += pg[idx] * cost
        elif gen.bid_type == 3:
            cost = gen.bids
            obj += (pg[idx] * pg[idx]) * cost[0] + pg[idx] * cost[1] + cost[2]
    # add line flow cons
    line_flow = {}
    for line_idx, line in enumerate(market.Line):
        line_flow[line_idx] = 0
        for bus_idx in range(market.Nb):
            load = market.load[bus_idx].P
            line_flow[line_idx] = line_flow[line_idx] + market.PTDF[line_idx, bus_idx] * (-load)
            line_flow[line_idx] = line_flow[line_idx] + market.PTDF[line_idx, bus_idx] * sum([pg[x] for x in sum(np.where(gen_bus == bus_idx + 1))])
        opt_model.addConstr(line_flow[line_idx] <= line.rating, name='TC p' + str(line_idx))
        opt_model.addConstr(line_flow[line_idx] >= -line.rating, name='TC n' + str(line_idx))
    # add power balance
    load_level = sum([np.sum(market.load[i].P) for i in range(market.load.__len__())])
    gen_level = sum([np.sum(pg[i]) for i in range(market.Ng)])
    opt_model.addConstr(gen_level == load_level, name="balance")
    opt_model.setObjective(obj, gb.GRB.MINIMIZE)
    opt_model.optimize()
    #opt_model.write('math_model.lp')
    for idx, gen in enumerate(market.genco):
        gen.opt_pg = pg[idx].X
    for idx, line in enumerate(market.Line):
        line.opt_fl = line_flow[idx].getValue()
    # LMP and dispatched settlements
    lamda = opt_model.getConstrByName('balance').Pi
    LMP = np.zeros([1, market.Nb])
    for b, ld in enumerate(market.load):
        LMP[0, b] = lamda
        for l, line in enumerate(market.Line):
            ng = opt_model.getConstrByName('TC n' + str(l)).Pi
            po = opt_model.getConstrByName('TC p' + str(l)).Pi
            LMP[0, b] += market.PTDF[l, b]*(ng-po)
    market.LMP = LMP
    for gen in market.genco:
        gen.revenue = market.LMP[0, int(gen.bus-1)]*gen.opt_pg
    for idx, ld in enumerate(market.load):
        ld.revenue = -market.LMP[0, idx]*ld.P
    del opt_model


def multi_ED(market):
    make_Bdc(market)
    make_PTDF(market)
    for t in range(market.N_T):
        opt_model = gb.Model(str(market.Type) + 'Ecnomic_dispatch')
        pg = {}
        load_level = 0
        obj = 0
        gen_bus = np.zeros([len(market.genco), 1])
        # add pg cap
        for idx, gen in enumerate(market.genco):
            gen_bus[idx] = gen.bus
            pg[idx] = opt_model.addVar(name='Power generation' + str(idx), vtype=gb.GRB.CONTINUOUS)
            if not gen.T_status:
                opt_model.addConstr(pg[idx] <= gen.pmax * gen.status, name='T_' + str(t) + 'Capacity_max' + str(idx))
                opt_model.addConstr(pg[idx] >= gen.pmin * gen.status, name='T_' + str(t) + 'Capacity_min' + str(idx))
            else:
                opt_model.addConstr(pg[idx] <= gen.pmax * gen.T_status[t], name='T_' + str(t) + 'Capacity_max' + str(idx))
                opt_model.addConstr(pg[idx] >= gen.pmin * gen.T_status[t], name='T_' + str(t) + 'Capacity_min' + str(idx))
            if gen.bid_type == 2:
                cost = gen.bids
                obj += pg[idx] * cost
            elif gen.bid_type == 3:
                cost = gen.bids
                obj += (pg[idx] * pg[idx]) * cost[0] + pg[idx] * cost[1] + cost[2]
        # add line flow cons
        line_flow = {}
        for line_idx, line in enumerate(market.Line):
            line_flow[line_idx] = 0
            for bus_idx in range(market.Nb):
                load = market.load[bus_idx].T_P[t]
                line_flow[line_idx] = line_flow[line_idx] + market.PTDF[line_idx, bus_idx] * (-load)
                line_flow[line_idx] = line_flow[line_idx] + market.PTDF[line_idx, bus_idx] *\
                                      sum([pg[x] for x in sum(np.where(gen_bus == bus_idx + 1))])
            opt_model.addConstr(line_flow[line_idx] <= line.rating, name='TC p' + str(line_idx))
            opt_model.addConstr(line_flow[line_idx] >= -line.rating, name='TC n' + str(line_idx))
        # add power balance
        load_level = market.load_level
        gen_level = sum([np.sum(pg[i]) for i in range(market.Ng)])
        opt_model.addConstr(gen_level == load_level[t], name="balance")
        opt_model.setObjective(obj, gb.GRB.MINIMIZE)
        opt_model.optimize()
        if opt_model.Status == 2:
            for idx, gen in enumerate(market.genco):
                gen.T_pg.append(pg[idx].X)
            for idx, line in enumerate(market.Line):
                line.T_Lf.append(line_flow[idx].getValue())
            # LMP and dispatched settlements
            lamda = opt_model.getConstrByName('balance').Pi
            LMP = np.zeros([1, market.Nb])
            for b, ld in enumerate(market.load):
                LMP[0, b] = lamda
                for l, line in enumerate(market.Line):
                    ng = opt_model.getConstrByName('TC n' + str(l)).Pi
                    po = opt_model.getConstrByName('TC p' + str(l)).Pi
                    LMP[0, b] += market.PTDF[l, b]*(ng-po)
            market.LMP = LMP
            market.T_LMP.append(LMP)
            for gen in market.genco:
                gen.T_revenue.append(LMP[0, int(gen.bus-1)]*gen.T_pg[t])
            for idx, ld in enumerate(market.load):
                ld.T_revenue.append(-LMP[0, idx]*ld.P)
        else:
            for idx, gen in enumerate(market.genco):
                gen.T_pg.append([])
                gen.T_revenue.append([])
            for idx, line in enumerate(market.Line):
                line.T_Lf.append([])
            for idx, ld in enumerate(market.load):
                ld.T_revenue.append([])
            market.T_LMP.append([])
            market.LMP = []
        del opt_model


def unit_commitment(market):
    make_Bdc(market)
    make_PTDF(market)
    opt_model = gb.Model(str(market.Type) + 'Unit Commitment')
    pg = {}
    status = {}
    comm_up = {}
    comm_down = {}
    obj = 0
    for t in range(market.N_T):
        gen_bus = np.zeros([len(market.genco), 1])
        # add pg cap and obj
        for idx, gen in enumerate(market.genco):
            gen_bus[idx] = gen.bus
            pg[t, idx] = opt_model.addVar(name='T_' + str(t) + 'P_g' + str(idx), vtype=gb.GRB.CONTINUOUS)
            status[t, idx] = opt_model.addVar(name='T_' + str(t) + 'Status' + str(idx), vtype=gb.GRB.BINARY)
            comm_up[t, idx] = opt_model.addVar(name='T_' + str(t) + 'comm_up' + str(idx), vtype=gb.GRB.BINARY)
            comm_down[t, idx] = opt_model.addVar(name='T_' + str(t) + 'comm_down' + str(idx), vtype=gb.GRB.BINARY)
            # commit key: if 1 then 0/1; if 2 then 1 Must_run_unit
            if gen.commit_key == 2:
                opt_model.addConstr(status[t, idx] == 1)
            elif gen.commit_key == 3:
                opt_model.addConstr(status[t, idx] == 0)
            # min down time
            if t <= gen.min_dowm-1:
                if gen.status == 0:
                    opt_model.addConstr(status[t, idx] == 0)
            elif t >= gen.min_dowm:
                if t != 0:
                    opt_model.addConstr((sum([comm_up[t-i, idx] for i in range(gen.min_dowm-1)])) <= 1-status[t-gen.min_dowm, idx])
            if t <= gen.min_up-2:
                if gen.status == 1:
                    opt_model.addConstr(status[t, idx] == 1)
            elif t >= gen.min_up-1:
                if t != 0:
                    opt_model.addConstr((sum([comm_up[t-i, idx] for i in range(gen.min_up-1)])) <= status[t, idx])
            # state transition constraint
            if t == 0:
                opt_model.addConstr(status[t, idx] - gen.status == comm_up[t, idx]-comm_down[t, idx])
            else:
                opt_model.addConstr(status[t, idx] - status[t-1, idx] == comm_up[t, idx]-comm_down[t, idx])
            # ramp
            if t >= 1:
                opt_model.addConstr(pg[t, idx] - pg[t-1, idx] <= gen.ramp_up,      #gen.T_ramp_up[t],
                                    name='T_' + str(t) + 'ramp_up' + str(idx))
                opt_model.addConstr(pg[t, idx] - pg[t-1, idx] >= -gen.ramp_down,        #gen.T_ramp_down[t],
                                    name='T_' + str(t) + 'ramp_down' + str(idx))
            # power capacity
            opt_model.addConstr(pg[t, idx] <= gen.pmax*status[t, idx], name='T_' + str(t) + 'Capacity_max' + str(idx))
            opt_model.addConstr(pg[t, idx] >= gen.pmin*status[t, idx], name='T_' + str(t) + 'Capacity_min' + str(idx))
            # obj
            cost = gen.bids
            no_load_c = 0      #gen.T_no_load[t]
            start_up_c = 0     #gen.T_no_load[t]
            obj += pg[t, idx]*cost+status[t, idx]*no_load_c+comm_up[t, idx]*gen.start_up+comm_down[t, idx]*gen.shut_down
        # add line flow cons
        line_flow = {}
        for line_idx, line in enumerate(market.Line):
            line_flow[t, line_idx] = 0
            for bus_idx in range(market.Nb):
                load = market.load[bus_idx].T_P[t]
                line_flow[t, line_idx] = line_flow[t, line_idx] + market.PTDF[line_idx, bus_idx] * (-load)
                line_flow[t, line_idx] = line_flow[t, line_idx] + market.PTDF[line_idx, bus_idx] *\
                                      sum([pg[t, x] for x in sum(np.where(gen_bus == bus_idx + 1))])
            opt_model.addConstr(line_flow[t, line_idx] <= line.rating, name='T_' + str(t) + 'TC_p' + str(line_idx))
            opt_model.addConstr(line_flow[t, line_idx] >= -line.rating, name='T_' + str(t) + 'TC_n' + str(line_idx))
        # add power balance
        load_level = market.load_level
        gen_level = sum([np.sum(pg[t, i]) for i in range(market.Ng)])
        opt_model.addConstr(gen_level == load_level[t], name='T_' + str(t) + 'balance')
    opt_model.update()
    opt_model.setObjective(obj, gb.GRB.MINIMIZE)
    opt_model.optimize()
    for idx, gen in enumerate(market.genco):
        for t in range(market.N_T):
            gen.T_status.append(int(status[t, idx].x))
    market.UC_result = 1


def real_time(market):
    for t in range(market.N_T):
        db.ex_ante_attack(market)
        for bus_idx in range(market.Nb):
            market.load[bus_idx].P = market.load[bus_idx].T_P[t]
        ecnomic_dispatch(market)
        P_cog_list = []
        N_cog_list = []
        for line_idx, line in enumerate(market.Line):
            if abs(line.opt_fl - line.rating) <= 0.00001:
                P_cog_list.append(line_idx)
            if abs(line.opt_fl+line.rating) <= 0.00001:
                N_cog_list.append(line_idx)
        db.ex_post_attack(market)
        d_pg_down = -2
        d_pg_up = 0.01
        opt_model = gb.Model(str(market.Type) + 'Ex_post_pricing')
        pg = {}
        obj = 0
        gen_bus = np.zeros([len(market.genco), 1])
        # add pg cap
        for idx, gen in enumerate(market.genco):
            gen_bus[idx] = gen.bus
            if gen.opt_pg >= -d_pg_down:
                pg[idx] = opt_model.addVar(name='Pg' + str(idx), vtype=gb.GRB.CONTINUOUS, lb=d_pg_down, ub=d_pg_up)
            else:
                pg[idx] = opt_model.addVar(name='Pg' + str(idx), vtype=gb.GRB.CONTINUOUS, lb=-gen.opt_pg, ub=d_pg_up)
            if gen.bid_type == 2:
                cost = gen.bids
                obj += pg[idx] * cost
            elif gen.bid_type == 3:
                cost = gen.bids
                obj += (pg[idx] * pg[idx]) * cost[0] + pg[idx] * cost[1] + cost[2]
        # add line flow cons
        line_flow = {}
        for line_idx, line in enumerate(market.Line):
            if line_idx in P_cog_list:
                line_flow[line_idx] = 0
                for bus_idx in range(market.Nb):
                    line_flow[line_idx] = line_flow[line_idx] + market.PTDF[line_idx, bus_idx] *\
                                          sum([pg[x] for x in sum(np.where(gen_bus == bus_idx + 1))])
                opt_model.addConstr(line_flow[line_idx] <= 0, name='TC p' + str(line_idx))
            if line_idx in N_cog_list:
                line_flow[line_idx] = 0
                for bus_idx in range(market.Nb):
                    line_flow[line_idx] = line_flow[line_idx] + market.PTDF[line_idx, bus_idx] * \
                                          sum([pg[x] for x in sum(np.where(gen_bus == bus_idx + 1))])
                opt_model.addConstr(line_flow[line_idx] >= 0, name='TC n' + str(line_idx))
        # add power balance
        gen_level = sum([np.sum(pg[i]) for i in range(market.Ng)])
        opt_model.addConstr(gen_level == 0, name="balance")
        opt_model.setObjective(obj, gb.GRB.MINIMIZE)
        opt_model.optimize()
        if opt_model.Status == 2:
            # LMP and dispatched settlements
            lamda = opt_model.getConstrByName('balance').Pi
            LMP = np.zeros([1, market.Nb])
            for b, ld in enumerate(market.load):
                LMP[0, b] = lamda
                po = 0
                ng = 0
                for l, line in enumerate(market.Line):
                    if l in P_cog_list:
                        po = opt_model.getConstrByName('TC p' + str(l)).Pi
                        LMP[0, b] += market.PTDF[l, b] * (-po)
                    if l in N_cog_list:
                        ng = opt_model.getConstrByName('TC n' + str(l)).Pi
                        LMP[0, b] += market.PTDF[l, b] * ng
            market.LMP = LMP
            market.RT_LMP.append(LMP)
            for gen in market.genco:
                gen.T_revenue.append(LMP[0, int(gen.bus-1)]*gen.T_pg[t])
            for idx, ld in enumerate(market.load):
                ld.T_revenue.append(-LMP[0, idx]*ld.P)
        else:
            for idx, gen in enumerate(market.genco):
                gen.T_pg.append([])
                gen.T_revenue.append([])
            for idx, line in enumerate(market.Line):
                line.T_Lf.append([])
            for idx, ld in enumerate(market.load):
                ld.T_revenue.append([])
            market.RT_LMP.append([])
            market.LMP = []
        del opt_model



