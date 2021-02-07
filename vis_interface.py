from dime import DimeClient


def send_LMP(market):
    d = DimeClient("tcp", "localhost", 5000)
    d.join("python")

    d["LMP"] = market.T_LMP
    d.send("python", "LMP")
    del d["LMP"]

    d.sync()
    print(d["LMP"])
