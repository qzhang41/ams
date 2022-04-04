import numpy as np
import gurobipy as gb
import Disturb as db
import data_output as ot
import multiprocessing


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

# def cons_lineflow(market,):
#     lineflow = np.dot(market.PTDF[i,:], np.matrix.transpose(G-np.matrix(L))) for i in range(market.Line.__len__())
#     return lineflow


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
    G = np.ndarray.astype(np.zeros([1, market.Nb]), object)
    L = np.array([market.load[idx].P for idx in range(market.Nb)])
    obj = 0
    # add pg cap
    for idx, gen in enumerate(market.genco):
        pg[idx] = opt_model.addVar(name='Power generation' + str(idx), vtype=gb.GRB.CONTINUOUS,
                                   ub=gen.pmax * gen.status, lb=gen.pmin * gen.status)
        opt_model.update()
        G[0, gen.bus-1] = G[0, gen.bus-1] + pg[idx]
        if gen.bid_type == 2:
            cost = gen.bids[0]
            obj += pg[idx] * cost + gen.bids[1]
        elif gen.bid_type == 3:
            cost = gen.bids
            obj += (pg[idx] * pg[idx]) * cost[0] + pg[idx] * cost[1] + cost[2]
    # add line flow cons
    # line_flow = [np.dot(market.PTDF[i,:], np.matrix.transpose(G-np.matrix(L))) for i in range(market.Line.__len__())]
    # G_L = np.matrix.transpose(G-np.matrix(L))
    # line_flow = market.PTDF * G_L
    G = np.matrix.transpose(np.matrix(G))
    L = np.matrix.transpose(L)
    line_flow = market.PTDF * (G-L)
    pos_con = [opt_model.addConstrs(line_flow[line_idx, 0] <= market.Line[line_idx].rating for line_idx in range(market.Line.__len__()))]
    neg_con = [opt_model.addConstrs(line_flow[line_idx, 0] >= -market.Line[line_idx].rating for line_idx in range(market.Line.__len__()))]
    opt_model.update()
    # add power balance
    load_level = sum([np.sum(market.load[i].P) for i in range(market.load.__len__())])
    gen_level = sum([np.sum(pg[i]) for i in range(market.Ng)])
    opt_model.addConstr(gen_level == load_level, name="balance")
    opt_model.setObjective(obj, gb.GRB.MINIMIZE)
    opt_model.optimize()
    # LMP and dispatched settlements
    for i in range(market.genco.__len__()):
        market.genco[i].opt_pg = pg[i].X
    for i in range(market.Line.__len__()):
        market.Line[i].opt_fl = line_flow[i, 0].getValue()
    lamda = opt_model.getConstrByName('balance').Pi
    LMP = np.ones([1, market.Nb]) * lamda
    gama = np.matrix([(abs(neg_con[0][l].Pi)-abs(pos_con[0][l].Pi)) for l in range(market.Line.__len__())])
    market.gama = gama.tolist()[0]
    d_LMP = gama * np.matrix(market.PTDF)
    LMP = LMP + d_LMP
    market.LMP = LMP
    del opt_model
    if market.Type == 'ED':
        if market.output == 1:
            ot.Out_to_CSV(market)
        elif market.output == 2:
            ot.Out_to_plot(market)

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
                cost = gen.bids[0]
                obj += pg[idx] * cost + gen.bids[1]
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
    opt_model.params.NonConvex = 2
    pg = {}
    status = {}
    comm_up = {}
    comm_down = {}
    obj = 0
    for t in range(market.N_T):
        gen_bus = np.zeros([len(market.genco), 1])
        # add pg cap and obj
        G = np.ndarray.astype(np.zeros([1, market.Nb]), object)
        L = np.array([market.load[idx].T_P[t] for idx in range(market.Nb)])
        for idx, gen in enumerate(market.genco):
            gen_bus[idx] = gen.bus
            pg[t, idx] = opt_model.addVar(name='T_' + str(t) + 'P_g' + str(idx), vtype=gb.GRB.CONTINUOUS)
            G[0, gen.bus - 1] = G[0, gen.bus - 1] + pg[t, idx]
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
            if gen.bid_type == 2:
                cost = gen.bids
                obj += (pg[t, idx] * cost[0] + cost[1])*status[t, idx] + status[t, idx]*no_load_c+comm_up[t, idx]*gen.start_up+comm_down[t, idx]*gen.shut_down
            elif gen.bid_type == 3:
                cost = gen.bids
                pg2 = opt_model.addVar(lb=0,vtype=gb.GRB.CONTINUOUS)
                opt_model.addConstr(pg2 == pg[t, idx]*pg[t, idx])
                obj += (pg2 * cost[0] + pg[t, idx] * cost[1] + cost[2])*status[t, idx]+status[t, idx]*no_load_c+comm_up[t, idx]*gen.start_up+comm_down[t, idx]*gen.shut_down
        # add line flow cons
        line_flow = market.PTDF * np.matrix.transpose(G - np.matrix(L))
        pos_con = [opt_model.addConstrs(
            line_flow[line_idx, 0] <= market.Line[line_idx].rating for line_idx in range(market.Line.__len__()))]
        neg_con = [opt_model.addConstrs(
            line_flow[line_idx, 0] >= -market.Line[line_idx].rating for line_idx in range(market.Line.__len__()))]
        # add power balance
        load_level = market.load_level
        gen_level = sum([np.sum(pg[t, i]) for i in range(market.Ng)])
        opt_model.addConstr(gen_level == load_level[t], name='T_' + str(t) + 'balance')
    opt_model.update()
    opt_model.setObjective(obj, gb.GRB.MINIMIZE)
    opt_model.optimize()
    for idx, gen in enumerate(market.genco):
        for t in range(market.N_T):
            gen.T_status.append(round(status[t, idx].x))
            gen.T_pg.append(pg[t,idx].X)
    if market.output == 1:
        ot.Out_to_CSV(market)
    elif market.output == 2:
        ot.Out_to_plot(market)

def real_time(market):
    for t in range(market.N_T):
        if market.attack == 1:
            if t == market.attack.t:
                db.ex_ante_attack(market)
        for bus_idx in range(market.Nb):
            market.load[bus_idx].P = market.load[bus_idx].T_P[t]
        ecnomic_dispatch(market)
        market.P_cog_list = []
        market.N_cog_list = []
        for line_idx, line in enumerate(market.Line):
            if abs(line.opt_fl - line.rating) <= 0.00001:
                market.P_cog_list.append(line_idx)
            if abs(line.opt_fl+line.rating) <= 0.00001:
                market.N_cog_list.append(line_idx)
        if market.attack == 1:
            if t == market.attack.t:
                db.ex_post_attack(market)
        P_cog_list = market.P_cog_list
        N_cog_list = market.N_cog_list
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
                cost = gen.bids[0]
                obj += pg[idx] * cost + gen.bids[1]
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
        else:
            market.RT_LMP.append([])
            market.LMP = []
        del opt_model
    if market.output == 1:
        ot.Out_to_CSV(market)
    elif market.output == 2:
        ot.Out_to_plot(market)

def CLL_collector(market):
    LMP_C = []
    Next_load_increase = 10000
    while Next_load_increase >= 1:
        ecnomic_dispatch(market)
        Next_load_increase = CLL(market)
        patern_f = [market.load[i].P for i in range(market.Nb)]
        patern_f = [patern_f[i] / sum(patern_f) for i in range(patern_f.__len__())]
        for i in range(market.load.__len__()):
            market.load[i].P += float(Next_load_increase) * patern_f[i]
        LMP_C.append(market.LMP)

def CLL(market):
    kc = 0.00001
    patern_f = [market.load[i].P for i in range(market.Nb)]
    patern_f = [patern_f[i]/sum(patern_f) for i in range(patern_f.__len__())]
    All_G = list(range(market.genco.__len__()))
    NG = []
    UL = []
    for i in range(market.genco.__len__()):
        if market.genco[i].opt_pg == market.genco[i].pmax or market.genco[i].opt_pg == market.genco[i].pmin:
            NG.append(i)
    for i in range(market.Line.__len__()):
        # if -market.Line[i].rating + kc <= market.Line[i].opt_fl <= market.Line[i].rating - kc:
        #     UL.append(i)
        if market.gama[i] == 0:
            UL.append(i)
    MG = list(set(All_G)-set(NG))
    CL = list(set(list(range(market.Nl))) - set(UL))
    MG_bus = [market.genco[i].bus-1 for i in MG]
    A00 = np.ones([1, MG.__len__()])
    A01 = np.zeros([1, UL.__len__()])
    A11 = np.matrix([market.PTDF[i, MG_bus] for i in CL])
    A12 = np.zeros([CL.__len__(), UL.__len__()])
    A21 = np.matrix([market.PTDF[i, MG_bus] for i in UL])
    A22 = np.identity(UL.__len__())
    A0 = np.concatenate((A00, A01), 1)
    A1 = np.concatenate((A11, A12), 1)
    A2 = np.concatenate((A21, A22), 1)
    A = np.concatenate((A0, A1), 0)
    A = np.concatenate((A, A2), 0)
    q0 = np.ones([1,market.Nb])
    q1 = np.matrix([market.PTDF[i, :] for i in CL])
    q2 = np.matrix([market.PTDF[i, :] for i in UL])
    q = np.concatenate((q0, q1), 0)
    q = np.concatenate((q, q2), 0)
    Q = np.linalg.inv(A)*q*np.matrix(patern_f).transpose()
    load_inc = []
    for i in range(MG.__len__()):
        if Q[i] >= 0:
            dis = market.genco[MG[i]].pmax - market.genco[MG[i]].opt_pg
            load_inc.append(dis/Q[i])
        if Q[i] <= -0.00000001:
            dis = market.genco[MG[i]].opt_pg - market.genco[MG[i]].pmin
            load_inc.append(dis/abs(Q[i]))
    for i in range(UL.__len__()):
        idx = i+MG.__len__()
        if Q[idx] >= 0:
            dis = market.Line[UL[i]].rating - market.Line[UL[i]].opt_fl
            load_inc.append(dis/Q[idx])
        if Q[idx] <= -0.00000001:
            dis = market.Line[UL[i]].opt_fl + market.Line[UL[i]].rating
            load_inc.append(dis/abs(Q[idx]))
    Next_load_increase = min(load_inc)
    Next_binding = load_inc.index(Next_load_increase)
    return Next_load_increase