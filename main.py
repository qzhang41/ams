# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
# Qiwei Zhang
import argparse
import Input_Parse
import Market_obj as MO
import Core


def market_gateway(market, Topo):
    Input_Parse.read_structure(market, Topo)

def Load_forecast(market, load_profile):
    Input_Parse.read_load(market, load_profile)
    market.Load_profile_flg = True


if __name__ == '__main__':
    # input parser
    parser = argparse.ArgumentParser(description='LTB Market Module')
    parser.add_argument('-ED', '--Economic Dispatch', help='Perform Economic dispatch')
    parser.add_argument('-L', '--Load Profile', help='Input Load Profile')
    parser.add_argument('-UC', '--Unit Commitment', help='Perform Unit commitment')
    parser.add_argument('-xgd', '--UC gen data', help='Input gen Profile')
    parser.add_argument('-WM', '--Whole Market-clearing', help='Perform Whole Market-clearing')
    args = parser.parse_args()
    args = vars(args)
    # Bid/Post system
    market = MO.Market('RT')
    if bool(args['Economic Dispatch']):
        market_gateway(market, args['Economic Dispatch'])
        if bool(args['Load Profile']):
            Load_forecast(market, args['Load Profile'])
            Core.Multi_Ecnomic_dispatch(market)
        else:
            Core.Ecnomic_dispatch(market)
    if bool(args['Unit Commitment']):
        market_gateway(market, args['Unit Commitment'])
        if bool(args['UC gen data']):
            Input_Parse.read_xgd(market, args['UC gen data'])
        if bool(args['Load Profile']):
            Load_forecast(market, args['Load Profile'])
            Core.Unit_commitment(market)
        else:
            Core.Ecnomic_dispatch(market)


