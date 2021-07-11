import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

def Out_to_CSV(market):
    type = market.Type
    if type == 'ED':
        path = 'C:/Users/qzhang41/PycharmProjects/Pmarket/Results_ED.xlsx'
        writer = pd.ExcelWriter(path, engine='xlsxwriter')
        disp = np.transpose(np.matrix([market.genco[i].opt_pg for i in range(market.genco.__len__())]))
        frame_LMP = pd.DataFrame(np.transpose(market.LMP), columns=['LMP'])
        frame_LMP.to_excel(writer, header=True, sheet_name='LMP')
        frame_DP = pd.DataFrame(disp, columns=['Dispatch'])
        frame_DP.to_excel(writer, header=True, sheet_name='Dispatch')
        writer.save()
        writer.close()
    if type == 'UC':
        path = 'C:/Users/qzhang41/PycharmProjects/Pmarket/Results_UC.xlsx'
        writer = pd.ExcelWriter(path, engine='xlsxwriter')
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
    if type == 'RT':
        path = 'C:/Users/qzhang41/PycharmProjects/Pmarket/Results_RT.xlsx'
        writer = pd.ExcelWriter(path, engine='xlsxwriter')
        frame_LMP = pd.DataFrame(np.transpose(market.RT_LMP), columns=['LMP'])
        frame_LMP.to_excel(writer, header=True, sheet_name='LMP')
        writer.save()
        writer.close()

def Out_to_plot(market):
    type = market.Type
    if type == 'ED':
        y = np.matrix.tolist(market.LMP)[0]
        x = list(range(market.Nb))
        plt.plot(x, y, color='r', marker='o', linestyle='dashed')
        plt.ylabel('LMP ($)')
        plt.xlabel('Bus')
        plt.title('LMP plot')
        plt.show()
    if type == 'UC':
        x = list(range(market.N_T))
        y1 = market.genco[market.p_unit].T_pg
        y2 = market.genco[market.p_unit].T_status
        plt.subplot(1,2,1)
        plt.plot(x, y1, color='r', marker='o', linestyle='dashed')
        plt.ylabel('Generation (MW)')
        plt.xlabel('Time period')
        plt.title('Dispatch')
        plt.subplot(1,2,2)
        plt.plot(x, y2, color='r', marker='o', linestyle='dashed')
        plt.ylabel('Status (0/1)')
        plt.xlabel('Time period')
        plt.title('Unit Status')
        plt.show()
    if type == 'RT':
        for i in range(market.N_T):
            x = [0, 1, 2]
            y = [0, 1, 2]
            r = random.random()
            b = random.random()
            g = random.random()
            color = (r, g, b)
            y = np.matrix.tolist(market.RT_LMP[i])[0]
            x = list(range(market.Nb))
            x = [int(x[i]+1) for i in range(market.Nb)]
            plt.plot(x, y, color=color, marker='o', linestyle='dashed', label='LMP at t'+str(i+1))
        plt.ylabel('LMP ($)')
        plt.xlabel('Bus')
        plt.title('LMP plot')
        plt.legend()
        plt.show()
