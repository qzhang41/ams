# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
# Qiwei Zhang
import argparse
import input_Parse
import market_obj as MO
import core
import streaming


if __name__ == '__main__':
    # input parser
    parser = argparse.ArgumentParser(description='LTB Market Module')
    parser.add_argument('-ED', '--Economic Dispatch', help='Perform Economic dispatch')
    parser.add_argument('-L', '--Load Profile', help='Input Load Profile')
    parser.add_argument('-UC', '--Unit Commitment', help='Perform Unit commitment')
    parser.add_argument('-XGD', '--UC gen data', help='Input gen Profile')
    parser.add_argument('-DM', '--Day ahead market', help='Perform Whole UC+ED')
    parser.add_argument('-DI', '--Dime', help='Send data to dime')
    args = parser.parse_args()
    args = vars(args)
    # Bid/Post system
    market = MO.Market('run UC ED')

    if bool(args['Dime']):
        market.dime = True
        port = int(args['Dime'])
        market.streaming = streaming.Streaming(market, port)

    if bool(args['Economic Dispatch']):
        input_Parse.read_structure(market, args['Economic Dispatch'])
        market.streaming.send_init()
        if bool(args['Load Profile']):
            input_Parse.read_load(market, args['Load Profile'])
            market.Load_profile_flg = True
            core.multi_ED(market)
        else:
            core.ecnomic_dispatch(market)
    if bool(args['Unit Commitment']):
        input_Parse.read_structure(market, args['Unit Commitment'])
        market.streaming.send_init()
        try:
            input_Parse.read_xgd(market, args['UC gen data'])
            input_Parse.read_load(market, args['Load Profile'])
        except:
            print("Error: fail to load xgd and load file")
        else:
            print("xgd and load file loaded")
            market.Load_profile_flg = True
            core.unit_commitment(market)
    if bool(args['Day ahead market']):
        input_Parse.read_structure(market, args['Day ahead market'])
        market.streaming.send_init()
        try:
            input_Parse.read_xgd(market, args['UC gen data'])
            input_Parse.read_load(market, args['Load Profile'])
        except:
            print("Error: fail to load xgd and load file")
        else:
            print("xgd and load file loaded")
            market.Load_profile_flg = True
            core.unit_commitment(market)
            core.multi_ED(market)
            market.streaming.finalize(market.streaming)
