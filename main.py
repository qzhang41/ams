# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
# Qiwei Zhang
import Input_Parse
import Market_obj as MO
import Core
import sys


def market_gateway(market):
    file_structure = sys.argv[1]
    Input_Parse.read_structure(market, file_structure)


def Load_forecast(market):
    if sys.argv.__len__() > 2:
        load_profile = sys.argv[2]
        Input_Parse.read_load(market, load_profile)
        market.Load_profile_flg = True


if __name__ == '__main__':
    # Bid/Post system
    market = MO.Market('RT')
    market_gateway(market)
    Load_forecast(market)
    # Dispatch
    if market.Load_profile_flg:
        Core.Multi_Ecnomic_dispatch(market)
    else:
        Core.Ecnomic_dispatch(market)