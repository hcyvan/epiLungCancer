import os
import yaml
from pathlib import Path


class Config:
    def __init__(self, config_path):
        self.config_abs_path = os.path.abspath(config_path)
        with open(self.config_abs_path, 'r') as file:
            data = yaml.safe_load(file)
            self.DataDir = Path(data.get('DataDir'))
            self.DataRaw = self.DataDir / 'raw'
            self.DataInter = self.DataDir / 'intermediate'
            self.DataResult = self.DataDir / 'result'

    def __str__(self):
        rets = []
        for k, v in vars(self).items():
            if 'Password' in k:
                v = "**********"
            rets.append("{}:{}".format(k, v))
        return "; ".join(rets)


CONFIG = Config(os.path.expanduser('./config.yaml'))
