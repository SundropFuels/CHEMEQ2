"""Processes datafiles created by the CEQ2 code and uses them to calculate various features"""

if __name__ == "__main__":

    temps = [1200,1300,1400,1500,1600]
    press = ["5.0"]
    fns = []
    for T in temps:
        for P in press:
            fns.append("methane_crack_hotwall_%s_%s.csv" % (T, P))

    for fn in fns:
        df = pd.from_csv(fn)
        df['C2H2_yield'] = df['C2H2']*2.0/df['CH4'][0]
        df['C2H4_yield'] = df['C2H4']*2.0/df['CH4'][0]
        df['CO_yield'] = df['CO']/df['CH4'][0]
        df['CO2_yield'] = df['CO2']/df['CH4'][0]
        df['C6H6_yield'] = df['C6H6']*6.0/df['CH4'][0]
        df['BIN1AB_yield'] = (df['BIN1A']+df['BIN1B'])*20.0/df['CH4'][0]
        df.to_csv("proc_" + fn)



