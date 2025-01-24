import pandas as pd
from pycode.config import CONFIG


class Sample:
    groups = ['CTL', 'LUAD', 'LUSC', 'LCC', 'SCLC','Healthy']

    def __init__(self):
        df = pd.read_excel(CONFIG.DataRaw / 'SupplementaryData.xlsx', sheet_name='sample')
        self.table = pd.concat([df[df['Group'] == g] for g in self.groups], axis=0)

    def get_samples(self, group=None):
        ret = self.table[['SampleNameTissue', 'Group']]
        if group:
            ret = ret[ret['Group'] == group]
        return ret

    def get_samples_(self):
        return self.table[['SampleNameTissue','SampleNameCfDNA', 'Group', 'Modeling','WGBS_data', 'WGBS_cfDNA']]

    def old_name2new(self, old_name):
        data_from = pd.DataFrame({'WGBS_data': old_name})
        data_to = pd.merge(data_from, self.table, on='WGBS_data', how='left')
        data_to = data_to[['WGBS_data', 'SampleName', 'Group']]
        return data_to['SampleName'].to_list()
