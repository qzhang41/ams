# Python based power market clearing model
# Author: Qiwei Zhang, Prof. Fangxing Li


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
    parser.add_argument('-RT', '--Real time market', help='Perform Whole Ex_ante+Ex_post')
    parser.add_argument('-OUT', '--output', help='to csv or plot')
    parser.add_argument('-PU', '--plot unit', help='plot unit')
    parser.add_argument('-DI', '--Dime', help='Send data to dime')
    parser.add_argument('-ATT', '--Attack', help='Real-time attack')
    args = parser.parse_args()
    args = vars(args)
    # Bid/Post system
    if bool(args['output']):
        output = args['output']
        if bool(args['plot unit']):
            plot_unit = args['plot unit']
    else:
        output = 0
    if bool(args['Economic Dispatch']):
        market = MO.Market('ED')
        market.output = int(output)
        input_Parse.read_structure(market, args['Economic Dispatch'])
        if market.dime:
            market.streaming.send_init()
        if bool(args['Load Profile']):
            input_Parse.read_load(market, args['Load Profile'])
            market.Load_profile_flg = True
            core.multi_ED(market)
        else:
            core.ecnomic_dispatch(market)
            core.CLL(market)
    if bool(args['Unit Commitment']):
        market = MO.Market('UC')
        market.output = int(output)
        market.p_unit = int(plot_unit)
        input_Parse.read_structure(market, args['Unit Commitment'])
        if market.dime:
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
        market = MO.Market('DA')
        market.output = int(output)
        market.p_unit = int(plot_unit)
        input_Parse.read_structure(market, args['Day ahead market'])
        if market.dime:
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
            # market.streaming.finalize(market.streaming)
    if bool(args['Real time market']):
        market = MO.Market('RT')
        market.output = int(output)
        input_Parse.read_structure(market, args['Real time market'])
        input_Parse.read_load(market, args['Load Profile'])
        market.Load_profile_flg = True
        if bool(args['Real time market']):
            market.attack = MO.Attack(market)
            input_Parse.read_attack(market, args['Attack'])
        core.real_time(market)

    if bool(args['Dime']):
        market.dime = True
        port = int(args['Dime'])
        market.streaming = streaming.Streaming(market, port)
