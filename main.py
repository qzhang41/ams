# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
# Qiwei Zhang
import Input_Parse
import Market_obj as MO
import Core
import sys


def market_gateway(market):
    file = sys.argv[1]
    Input_Parse.read(market, file)


if __name__ == '__main__':
    market = MO.Market('RT')
    market_gateway(market)
    Core.Make_Bdc(market)
    Core.Make_PTDF(market)
    Core.Ecnomic_dispatch(market)