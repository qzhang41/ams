import pandas as pd
import numpy as np

def Out_to_CSV(market):
    type = market.Type
    path = 'C:/Users/qzhang41/PycharmProjects/Pmarket/Results.xlsx'
    writer = pd.ExcelWriter(path, engine='xlsxwriter')
    if type == 'ED':
        disp = np.transpose(np.matrix([market.genco[i].opt_pg for i in range(market.genco.__len__())]))
        frame_LMP = pd.DataFrame(np.transpose(market.LMP), columns=['LMP'])
        frame_LMP.to_excel(writer, header=True, sheet_name='LMP')
        frame_DP = pd.DataFrame(disp, columns=['Dispatch'])
        frame_DP.to_excel(writer, header=True, sheet_name='Dispatch')
        writer.save()
        writer.close()
    if type == 'UC':
        for i in range(market.genco.__len__()):
            if i == 0:
                T_pg = np.matrix(market.genco[i].T_pg)
                T_status = np.matrix(market.genco[i].T_status)
            else:
                T_pg = np.concatenate((T_pg,np.matrix(market.genco[i].T_pg)))
                T_status = np.concatenate((T_status, np.matrix(market.genco[i].T_status)))
        header = list(range(0, market.N_T))
        header = [str(header[i]) for i in range(market.N_T)]
        frame_pg = pd.DataFrame(T_pg, columns=header)
        frame_pg.to_excel(writer, header=True, sheet_name='T_pg')
        frame_st = pd.DataFrame(T_status, columns=header)
        frame_st.to_excel(writer, header=True, sheet_name='T_status')
        writer.save()
        writer.close()