class Market:
    def __init__(self, market_type):
        self.Type = market_type
        self.load = []
        self.Ng = 0
        self.Nl = 0
        self.Nb = 0
        self.GSF = []
        self.sw = 0
        self.genco = []
        self.Line = []
        self.Bbus = []
        self.Bf = []
        self.Cft = []
        self.PTDF = []
        self.LMP = []
        self.Load_profile_flg = False
        self.T_LMP = []
        self.load_level = []


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
        self.bids = 0
        self.pmax = 0
        self.pmin = 0
        self.status = 0
        self.opt_pg = 0
        self.Revenue = 0
        self.T_Pg = []
        self.T_Revenue = []


class Load:
    def __init__(self, market_type):
        self.Type = market_type
        self.P = 0
        self.Revenue = []
        self.T_P = []
        self.T_Revenue = []
