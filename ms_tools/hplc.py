from rainbow.agilent.chemstation import parse_ch
import pandas as pd
import numpy as np
from glob import glob
import os

class ch_sample:
    def __init__(self, path: str):
        """Load data from .ch file

        Args:
            path (str): Path to .ch file
        """
        self.path = path
        
        ch = parse_ch(path)
        
        self._ch = ch
        self.time = ch.xlabels
        self.riu = ch.data.flatten()
        self.metadata = ch.metadata
        self.date = ch.metadata['date']
        self.sample = ch.metadata['notebook']
        self.data = np.array([self.time, self.riu])
        
    def __repr__(self) -> str:
        return f'.ch file "{self.path}" for sample "{self.sample}"'

class ch_directory:
    def __init__(self, path: str):
        """Load all .ch files from a directory
        
        Args:
            path (str): Path to directory
        """
        
        if not os.path.isdir(path):
            raise TypeError(f'{path} must be a directory.')
        
        files = glob(os.path.join(path, '**/*.ch'), recursive=True)
        
        if not files:
            raise FileNotFoundError(f'{path} does not contain any .ch files within it.')
        
        self.path = path
        self.files = files
        self._n_files = len(files)
        self.chs = []
        self.df = None
        
        self._load_ch_data()
        self._generate_df()
        
    def _load_ch_data(self):
        "Load all .ch files in `path`"
        for path in self.files:
            ch = ch_sample(path)
            self.chs.append(ch)
            
    def _generate_df(self):
        "Compile data into a single longform dataframe"
        data = [pd.DataFrame(
                    {'Sample': ch.sample,
                    'time': ch.time, 
                    'riu': ch.riu,
                    'date': ch.date}
                    ) for ch in self.chs]
        
        self.df = pd.concat(data)
        
    def __repr__(self):
        return 'Data for {self._n_files} *.ch files in "{self.path}"'