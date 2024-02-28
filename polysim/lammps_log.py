import numpy

class ThermoData:
    ''' Holds thermo data from a LAMMPS run step. '''
    def __init__(self, columns):
        self.columns = columns.split()
        self.data = []

    def add_row(self, thermo):
        ''' Adds a row of thermo data to this step. '''
        x = [float(s) for s in thermo.split()]
        if len(x) == len(self.columns):
            self.data.append(x)

    def __getitem__(self, name):
        ''' '''
        if name in self.columns:
            return self.data[:, self.columns.index(name)]
        else:
            print(f'Missing thermo data {name}')
            return None

def read_log(path):
    ''' Reads a LAMMPS log file '''

    fid = open(path)
    steps = []
    while True:
        for line in fid:
            if 'Step' in line:
                steps.append(ThermoData(line))
                break
        else:
            break
        for line in fid:
            if 'Loop time' in line:
                print(f'Read {len(steps[-1].data)} log steps')
                break
            else:
                steps[-1].add_row(line)
        else:
            break

    for step in steps:
        step.data = numpy.array(step.data)
    return steps

