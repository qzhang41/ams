

def load_error(market):
    return


def ex_ante_attack(market):
    for i in range(market.attack.Lr.__len__()):
        idx = int(market.attack.Lr[i][0])
        value = market.attack.Lr[i][1]
        market.Line[idx].rating += value
    for i in range(market.attack.Load.__len__()):
        idx = int(market.attack.Load[i][0])
        value = market.attack.Load[i][1]
        market.load[idx].P += value
    for i in range(market.attack.Bidding.__len__()):
        idx = int(market.attack.Bidding[i][0])
        value = market.attack.Bidding[i][1]
        market.genco[idx].bids[0] += value
    return


def ex_post_attack(market):
    market.N_cog_list = market.attack.Cog_n
    market.P_cog_list = market.attack.Cog_p
    return
