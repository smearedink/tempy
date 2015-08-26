from collections import OrderedDict
from _settings import string_pars

# Only used by functions here, so defined here rather than in classes.py
class Par:

    def __init__(self,name,value,error=None,fit=None):
        self.name = name
        self.value = value
        self.error = error
        self.fit = fit
    def __repr__(self):
        return "Name: " + self.name + " Value: " + str(self.value) + \
               " Error: " + str(self.error) + " Fit: " + str(self.fit)

def read_parfile(parfn):
    pars = OrderedDict()
    parf = open(parfn,'r')
    for line in parf.readlines():
        split = line.split()
        if split[0] == 'C': #its a comment
            continue
        if len(split) == 1: #has no value
            continue

        value = split[1] if split[0] in string_pars else float(split[1].replace('D', 'E'))

        if len(split) == 2:
            pars[split[0]] = Par(split[0],value)
        elif len(split) == 3:
            if (split[2] == '0') or (split[2] == '1'):
                pars[split[0]] = Par(split[0],value,fit=int(split[2]))
            else:
                pars[split[0]] = Par(split[0],value,error=float(split[2]))
        elif len(split) == 4:
            pars[split[0]] = Par(split[0],value,error=float(split[3].replace('D', 'E')), \
                                 fit=split[2])
    return pars

def write_parfile(pars, outnm):
    with open (outnm, 'w') as f:
        for item in pars:
            if pars[item].fit!=None:
                f.write( item +' '+str(pars[item].value)+'  '+ str(pars[item].fit)+ '\n')
            else:
                f.write( item +' '+str(pars[item].value)+'  '+ '\n')
