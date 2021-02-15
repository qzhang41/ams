import logging
from time import sleep
import numpy as np

logger = logging.getLogger(__name__)

try:
    from dime import DimeClient
except ImportError:
    logger.debug("Dime import failed.")


class Streaming:

    def __init__(self, market, port):
        self.params_built = False
        self.market = market
        self.dimec = None
        self.SysParam = dict()
        self.LMP = dict()
        self.port = port
        self.recipent = "Vis"

    def connect(self):
        self.dimec = DimeClient("tcp", "localhost", self.port)
        self.dimec.join("market")

    def build_init(self):
        bus_name = []
        latitude = []
        longitude = []
        for bus in self.market.bus:
            bus_name.append(bus.name)
            latitude.append(bus.latitude)
            longitude.append(bus.longitude)
        self.dimec["Bus"] = bus_name
        self.dimec["latitude"] = np.array(latitude)
        self.dimec["longitude"] = np.array(longitude)

    def send_init(self):
        if not self.market.dime:
            return
        self.connect()
        if not self.params_built:
            self.build_init()
            self.params_built = True
        t_sleep = 0.05
        self.dimec2 = DimeClient("tcp", "localhost", self.port)
        self.dimec2.join("Vis")
        sleep(t_sleep)
        self.dimec.send(self.recipent, "Bus")
        sleep(t_sleep)
        self.dimec2.sync()
        self.dimec.send(self.recipent, "latitude")
        sleep(t_sleep)
        self.dimec2.sync()
        self.dimec.send(self.recipent, "longitude")
        sleep(t_sleep)
        self.dimec2.sync()
        sleep(t_sleep)

    def build_LMP(self):
        self.dimec["LMP"] = np.array(self.market.LMP)

    def send_LMP(self):
        if not self.market.dime:
            return
        self.build_LMP()
        sleep(0.05)
        self.dimec.send(self.recipent, "LMP")
        sleep(0.05)
        self.dimec2.sync()
        print(self.dimec2["LMP"])

    @staticmethod
    def finalize(self):
        """
        Send ``DONE`` signal when simulation completes
        :return: None
        """
        if not self.market.dime:
            return
        self.dimec.broadcast_r(DONE=1)
        self.dimec.close()
        self.dimec2.close()