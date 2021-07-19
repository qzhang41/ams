import numpy as np


class Market:
    def __init__(self, market_type):
        self.Type = market_type
        self.bus = []
        self.load = []
        self.Ng = 0
        self.Nl = 0
        self.Nb = 0
        self.sw = 0
        self.genco = []
        self.Line = []
        self.Bbus = []
        self.Bf = []
        self.Cft = []
        self.PTDF = []
        self.LMP = []
        self.RT_LMP = []
        self.Load_profile_flg = False
        self.T_LMP = []
        self.load_level = []
        self.action = []
        self.UC_result = 0
        self.dime = False
        self.DA_rs = []
        self.RT_LMP = []
        self.output = 0
        self.p_unit = -1
        self.disturb = False
        self.attack = []

class Bus:
    def __init__(self, market_type):
        self.Type = market_type
        self.name = ""
        self.latitude = 0
        self.longitude = 0


class Line:
    def __init__(self, market_type):
        self.Type = market_type
        self.fbus = 0
        self.tbus = 0
        self.rating = 0
        self.status = 0
        self.r = 0
        self.x = 0
        self.b = 0
        self.opt_fl = 0
        self.T_Lf = []


class Genco:
    def __init__(self, market_type):
        self.Type = market_type
        self.bus = 0
        self.bids = []
        self.bid_type = 2
        self.pmax = 0
        self.pmin = 0
        self.status = 0
        self.opt_pg = 0
        self.revenue = 0
        self.start_up = 0
        self.shut_down = 0
        self.ramp_up = 0
        self.ramp_down = 0
        self.ramp_up_price = 0
        self.ramp_down_price = 0
        self.commit_key = 0
        self.min_dowm = 0
        self.min_up = 0
        self.T_no_load = []
        self.T_pg = []
        self.T_status = []
        self.T_revenue = []


class Load:
    def __init__(self, market_type):
        self.Type = market_type
        self.P = 0
        self.revenue = []
        self.T_P = []
        self.T_revenue = []

class Attack:
    def __init__(self, market):
        self.Lr = []
        self.Load = []
        self.Cog_p = []
        self.Cog_n = []
        self.Bidding = []
        self.t = -1
