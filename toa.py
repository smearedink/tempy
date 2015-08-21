from astropy import time as _atime

class TOA:
    """
    A class representing a single pulsar time of arrival for use with TEMPO.

    A TOA object contains the following information:
      MJD ----- Precise time of arrival, stored in astropy format
      err ----- Error on MJD in microseconds
      freq ---- Radio frequency this arrival time is associated with
      obs ----- Observatory (TEMPO's observatory code, '@' is barycentre)
      dm_corr - Optional offset in dispersion measure for this TOA
      label --- Optional label for this TOA

    Initialization:
      (1) toa = TOA(MJDi, MJDf, err, freq, obs, dm_corr, label)
        MJDi and MJDf are the integer and float (following the decimal point)
        parts of the MJD, respectively
      (2) toa = TOA.from_princeton_format(toa_str)
        toa_str is a line from a Princeton-format TEMPO .tim file

    Methods:
      to_princeton_format()
        Return a string formatted as a line from a Princeton-format
        TEMPO .tim file
    """
    def __init__(self, MJDi, MJDf, err, freq, obs, dm_corr=0, label=''):
        self.MJD = _atime.Time(MJDi, MJDf, format='mjd')
        self.err = float(err)
        self.freq = float(freq)
        self.obs = obs
        self.dm_corr = float(dm_corr)
        self.label = label

    def __repr__(self):
        return "<TOA: %d.%s>" % (int(self.MJD.mjd),
                                 ("%.13f" % (self.MJD.jd2 % 1))[2:])

    def __lt__(self, other):
        return self.MJD < other.MJD
    def __le__(self, other):
        return self.MJD <= other.MJD
    def __gt__(self, other):
        return self.MJD > other.MJD
    def __ge__(self, other):
        return self.MJD >= other.MJD

    @classmethod
    def from_princeton_format(cls, toa_str):
        obs = toa_str[0:1]
        label = toa_str[2:15].strip()
        freq = float(toa_str[15:24])
        mjd_str = toa_str[24:44].split('.')
        MJDi = int(mjd_str[0])
        MJDf = float('.' + mjd_str[1])
        err = float(toa_str[44:53])
        dm_corr_str = toa_str[68:78]
        if not dm_corr_str.strip():
            return cls(MJDi, MJDf, err, freq, obs, 0, label)
        else:
            return cls(MJDi, MJDf, err, freq, obs, float(dm_corr_str), label)

    def to_princeton_format(self):
        toa = "%5d"%int(self.MJD.mjd) + ("%.13f" % (self.MJD.jd2 % 1))[1:]
        if self.dm_corr != 0:
            return self.obs+" %13s %8.3f %s %8.3f              %9.4f" % \
                   (self.label, self.freq, toa, self.err, self.dm_corr)
        else:
            return self.obs+" %13s %8.3f %s %8.3f" % \
                   (self.label, self.freq, toa, self.err)



class TOAset:
    """
    A class representing a set of TOAs used for timing in TEMPO, typically
    stored in and read from a .tim file.

    A TOAset object contains the following information:
      TOAs -------- List of TOAs in the order they appear in a .tim file
      jump_ranges - List of tuples of the form (start, end) where
                    'start' and 'end' are indices referring to sections of the
                    'TOAs' list bracketed by TEMPO JUMP statements
                    (TOAS[start:stop] would be the list of those bracketed)
      phase_wraps - Dictionary where each key is an index in the 'TOAs' list
                    before which a TEMPO PHASE statement occurs, and the value
                    is the argument to the PHASE statement
      mode -------- If set, the argument to a .tim file MODE statement
      track ------- If set, the argument to a .ti file TRACK statement

    Initialization:
      (1) toa_set = TOAset(TOAs, jump_ranges, phase_wraps, mode, track)
      (2) toa_set = TOAset.from_princeton_file(fname)
        fname is the path to a TEMPO .tim file

    Methods:
      get_TOAs_from_jump_range(index)
        Returns a list of the TOAs from the jump range identified by the tuple
        at position 'index' in the 'jump_ranges' list
      jump_statement_before(index)
        Returns 'True' if a JUMP statement appears before the TOA at index
      jump_statement_after(index)
        Returns 'True' if a JUMP statement appears after the TOA at index
      get_nTOAs()
        Returns the total number of TOAs in this TOAset
      get_span()
        Returns the total number of days spanned by this TOAset
      to_princeton_file(fname)
        If fname provided, outputs TOAset to a new .tim file in Princeton
        format; otherwise, prints to screen in Princeton format
    """
    def __init__(self, list_of_TOAs=[], jump_ranges=[], phase_wraps={},
                 mode=None, track=None):
        self.TOAs = list_of_TOAs
        self.jump_ranges = jump_ranges
        self.phase_wraps = phase_wraps
        self.mode = mode
        self.track = track

    def __repr__(self):
        return "<TOAset: %d TOAs over %.1f days>" %\
          (self.get_nTOAs(), self.get_span())

    @classmethod
    def from_princeton_file(cls, fname):
        TOAs = []
        jump_ranges = []
        jump = False
        phase_wraps = {}
        mode = None
        track = None
        with open(fname, 'r') as f:
            for line in f.readlines():
                if line[0].upper() != 'C':
                    if line.strip() == "JUMP":
                        jump = not jump
                        if jump:
                            jump_ranges.append([len(TOAs), 0])
                        else:
                            jump_ranges[-1][1] = len(TOAs)
                            jump_ranges[-1] = tuple(jump_ranges[-1])
                    elif line.strip()[:4] == "MODE":
                        mode = int(line.split()[1])
                    elif line.strip()[:5] == "TRACK":
                        track = int(line.split()[1])
                    elif line.strip()[:5] == "PHASE":
                        phase_wrap = int(line.strip()[5:])
                        phase_wraps[len(TOAs)] = phase_wrap
                    elif len(line.strip()) >= 50:
                        TOAs.append(TOA.from_princeton_format(line.strip()))
        return cls(TOAs, jump_ranges, phase_wraps, mode, track)

    def get_TOAs_from_jump_range(self, index):
        return self.TOAs[self.jump_ranges[index][0]:self.jump_ranges[index][1]]

    def jump_statement_before(self, index):
        return index in [r[0] for r in self.jump_ranges]

    def jump_statement_after(self, index):
        return index+1 in [r[1] for r in self.jump_ranges]

    def get_nTOAs(self):
        return len(self.TOAs)

    def get_span(self):
        mjds = [t.MJD for t in self.TOAs]
        first_mjd = min(mjds)
        last_mjd = max(mjds)
        return (last_mjd - first_mjd).jd

    def to_princeton_file(self, fname=None):
        """
        If filename not provided, formatted TOAs are printed to screen
        """
        lines = []
        if self.mode is not None:
            lines.append("MODE %d" % self.mode)
        if self.track is not None:
            lines.append("TRACK %d" % self.track)
        for ii in range(self.get_nTOAs()):
            if ii in self.phase_wraps:
                if self.phase_wraps[ii] > 0:
                    ph_arg = "+" + str(self.phase_wraps[ii])
                else:
                    ph_arg = str(self.phase_wraps[ii])
                lines.append("PHASE %s" % ph_arg)
            if self.jump_statement_before(ii):
                lines.append("JUMP")
            lines.append(self.TOAs[ii].to_princeton_format())
            if self.jump_statement_after(ii):
                lines.append("JUMP")
        if fname is None:
            for line in lines:
                print line
        else:
            with open(fname, 'w') as f:
                for line in lines:
                    f.write(line + "\n")


string_pars = ['PSR','PSRJ','EPHEM','CLK','BINARY','JUMP',\
               'UNITS','TZRSITE','TIMEEPH','T2CMETHOD',\
               'CORRECT_TROPOSPHERE','PLANET_SHAPIRO','DILATEFREQ',\
               'RAJ', 'DECJ', 'NITS']

from collections import OrderedDict

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
