class Market:
    def __init__(self, market_type):
        self.Type = market_type
    load = []
    Ng = 0
    Nl = 0
    Nb = 0
    GSF = []
    sw = 0
    genco = []
    Line = []
    Bbus = []
    Bf = []
    Cft = []
    PTDF = []
    LMP = []
    Load_profile_flg = False
    T_LMP = []
    load_level = []

class Line:
    def __init__(self, market_type):
        self.Type = market_type
    fbus = 0
    tbus = 0
    rating = 0
    status = 0
    r = 0
    x = 0
    b = 0
    Opt_fl = 0
    T_lf = []


class Genco:
    def __init__(self, market_type):
        self.Type = market_type
    bus = 0
    bids = 0
    pmax = 0
    pmin = 0
    status = 0
    opt_pg = 0
    Revenue = 0
    # from load profile
    T_Pg = []
    T_Revenue = []



class Load:
    def __init__(self, market_type):
        self.Type = market_type
    P = 0
    Revenue = []
    # from load profile
    T_P = []
    T_Revenue = []