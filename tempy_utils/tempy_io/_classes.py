from astropy import time as _atime

class TOA:
    """
    A class representing a single pulsar time of arrival for use with TEMPO.
    A TOA object contains the following information:
      MJD ---------- Precise time of arrival, stored in astropy format
      err ---------- Error on MJD in microseconds
      freq --------- Radio frequency this arrival time is associated with
      obs ---------- Observatory (TEMPO's observatory code, '@' is barycentre)
      dm_corr ------ Optional offset in dispersion measure for this TOA
      label -------- Optional label for this TOA
      phase_offset - Optional offset in phase for this TOA

    Initialization:
      (1) toa = TOA(MJDi, MJDf, err, freq, obs, dm_corr, label, phase_offset)
        MJDi and MJDf are the integer and float (following the decimal point)
        parts of the MJD, respectively
      (2) toa = TOA.from_[princeton/parkes/ITOA]_format(toa_str)
        toa_str is a line from a TEMPO .tim file

    Methods:
      to_princeton_format()
        Return a string formatted as a line from a Princeton-format
        TEMPO .tim file
    """
    def __init__(self, MJDi, MJDf, err, freq, obs, dm_corr=0, label='', phase_offset=0):
        self.MJD = _atime.Time(MJDi, MJDf, format='mjd')
        self.err = float(err)
        self.freq = float(freq)
        self.obs = obs
        self.dm_corr = dm_corr
        if self.dm_corr is not None:
            self.dm_corr = float(self.dm_corr)
        self.label = label
        self.phase_offset = phase_offset

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
        if self.dm_corr is not None:
            if self.dm_corr != 0:
                return self.obs+" %13s %8.3f %s %8.3f              %9.4f" % \
                       (self.label, self.freq, toa, self.err, self.dm_corr)
            else:
                return self.obs+" %13s %8.3f %s %8.3f" % \
                       (self.label, self.freq, toa, self.err)
        else:
            return self.obs+" %13s %8.3f %s %8.3f" % \
                   (self.label, self.freq, toa, self.err)

    @classmethod
    def from_parkes_format(cls, toa_str):
        obs = toa_str[79:80]
        label = toa_str[2:16].strip()
        freq = float(toa_str[25:34])
        mjd_str = toa_str[34:55].split('.')
        MJDi = int(mjd_str[0])
        MJDf = float('.' + mjd_str[1])
        phase_offset_str = toa_str[55:63]
        err = float(toa_str[63:71])
        dm_corr=None
        if not phase_offset_str.strip():
            return cls(MJDi, MJDf, err, freq, obs, dm_corr, label, 0)
        else:
            return cls(MJDi, MJDf, err, freq, obs, dm_corr, label,
                       float(phase_offset_str))

    def to_parkes_format(self):
        toa = "%5d"%int(self.MJD.mjd) + ("%.13f" % (self.MJD.jd2 % 1))[1:]
        return "%12s              %8.3f  %s%7.6f  %8.3f      " % \
          (self.label, self.freq, toa, self.phase_offset, self.err)+self.obs

    @classmethod
    def from_ITOA_format(cls, toa_str):
        obs = toa_str[57:59]
        freq = float(toa_str[34:45])
        mjd_str = toa_str[9:28].split('.')
        MJDi = int(mjd_str[0])
        MJDf = float('.' + mjd_str[1])
        err = float(toa_str[28:34])
        dm_corr_str = toa_str[45:55]
        phase_offset = None
        if not dm_corr_str.strip():
            return cls(MJDi, MJDf, err, freq, obs, 0, '', phase_offset)
        else:
            return cls(MJDi, MJDf, err, freq, obs, float(dm_corr_str), '',
                       phase_offset)

    def to_ITOA_format(self):
        toa = "%5d"%int(self.MJD.mjd) + ("%.13f" % (self.MJD.jd2 % 1))[1:]
        # TODO: needs to have spacing fixed
        return " %8s %s %8.3f %8.3f %9.4f   " % \
          (self.label, toa, self.err, self.freq, self.dm_corr) + self.obs
    
    @classmethod
    def from_Tempo2_format(cls, toa_str):
        split_toa_str = toa_str.split()[1:]
        if (len(split_toa_str) == 5) and (len(split_toa_str[-1]) == 1):
            split_toa_str = split_toa_str[1:]
        freq = float(split_toa_str[0])
        mjd_str = split_toa_str[1].split('.')
        MJDi = int(mjd_str[0])
        MJDf = float('.' + mjd_str[1])
        err = float(split_toa_str[2])
        obs = split_toa_str[3]
        dm_corr = None
        phase_offset = None
        return cls(MJDi, MJDf, err, freq, obs, dm_corr, '', phase_offset)

    def to_Tempo2_format(self):
        toa = "%5d"%int(self.MJD.mjd) + ("%.13f" % (self.MJD.jd2 % 1))[1:]
        return " NOT %8.3f %s %8.3f " % \
          (self.freq, toa, self.err)+self.obs
        

class TOAset:
    def __init__(self, list_of_TOAs=[], jumped=False):
        self.TOAs = list(list_of_TOAs)
        self.jumped = jumped

    def append_TOA(self, toa):
        self.TOAs.append(toa)

    def get_nTOAs(self):
        return len(self.TOAs)



class TOAfile:
    """
    A class representing many TOAs used for timing in TEMPO, typically
    stored in and read from a .tim file.
    """
    def __init__(self, TOAsets=[], phase_wraps={}, mode=None, track=None,
                 default_format='princeton'):
        self.TOAsets = TOAsets
        # phase_wraps keys are tuple (TOAsets index, TOAset.TOAs index)
        # and values are phase wrap value (eg, -2)
        self.phase_wraps = phase_wraps
        self.mode = mode
        self.track = track
        self.default_format = default_format
        self.sorted_TOAs = []

    @classmethod
    def from_tim_file(cls, fname):
        TOAsets = [TOAset(jumped=False)]
        TOAset_index = 0
        jump = False
        phase_wraps = {}
        mode = None
        track = None
        with open(fname, 'r') as f:
            for line in f.readlines():
                if line[0].upper() != 'C':
                    use_Tempo2 = False
                    if line.strip() == "JUMP":
                        jump = not jump
                        if TOAsets[-1].get_nTOAs() == 0:
                            TOAsets[-1].jumped = jump
                        else:
                            TOAsets.append(TOAset(jumped=jump))
                            TOAset_index += 1
                    elif line.strip()[:4] == "MODE":
                        mode = int(line.split()[1])
                    elif line.strip()[:5] == "TRACK":
                        track = int(line.split()[1])
                    elif line.strip()[:5] == "PHASE":
                        phase_wrap = int(line.strip()[5:])
                        phase_wraps[(TOAset_index, TOAsets[-1].get_nTOAs())] \
                          = phase_wrap
                    elif line.strip()[:6] == "FORMAT":
                        if line.split()[1] == "1":
                            default_format = 'Tempo2'
                            use_Tempo2 = True
                    elif len(line) >= 20:
                        err_str = "Unsupported file format. Format must be " \
                                  "princeton, parkes, Tempo2, or ITOA."
                        if use_Tempo2:
                            TOAsets[-1].append_TOA(TOA.from_Tempo2_format(line))
                        else:
                            try:
                                toa = TOA.from_princeton_format(line)
                                default_format = 'princeton'
                            except:
                                try:
                                    toa = TOA.from_parkes_format(line)
                                    default_format = 'parkes'
                                except:
                                    try:
                                        toa = TOA.from_ITOA_format(line)
                                        default_format = 'ITOA'
                                    except:
                                        print err_str
                                        exit(1)
                            TOAsets[-1].append_TOA(toa)
        return cls(TOAsets, phase_wraps, mode, track, default_format)

    def get_nTOAsets(self):
        return len(self.TOAsets)

    def get_nTOAs(self, before_TOAset_index=None):
        """
        If before_TOAset_index is not None, this will return the number of TOAs
        before that set.  Otherwise it returns the total number of TOAs.
        """
        if before_TOAset_index is None:
            return sum([tset.get_nTOAs() for tset in self.TOAsets])
        else:
            return sum([tset.get_nTOAs() for tset in
                       self.TOAsets[:before_TOAset_index]])

    def get_position_of_TOA(self, n):
        TOAset_index = 1
        while self.get_nTOAs(TOAset_index) <= n:
            TOAset_index += 1
        return (TOAset_index-1, n-self.get_nTOAs(TOAset_index-1))
        
    def get_ordered_index_of_position(self, pos):
        """
        pos is a tuple: (index of TOAset, index of TOA in TOAset)
        
        Returns what index would be for this TOA in an ordered list of TOAs
        """
        TOA_at_pos = self.TOAsets[pos[0]].TOAs[pos[1]]
        all_TOAs = []
        for tset in self.TOAsets:
            for toa in tset.TOAs:
                all_TOAs.append(toa)
        if len(self.sorted_TOAs) != self.get_nTOAs():
            self.sorted_TOAs = sorted(all_TOAs)
        return self.sorted_TOAs.index(TOA_at_pos)

    def get_jump_ranges(self, chronological=False):
        jump_ranges = []
        for TOAset_index in range(self.get_nTOAsets()):
            if self.TOAsets[TOAset_index].jumped:
                jstart = self.get_nTOAs(TOAset_index)
                jend = jstart + self.TOAsets[TOAset_index].get_nTOAs()-1
                if chronological:
                    jstart_pos = self.get_position_of_TOA(jstart)
                    jend_pos = self.get_position_of_TOA(jend)
                    jstart = self.get_ordered_index_of_position(jstart_pos)
                    jend = self.get_ordered_index_of_position(jend_pos)
                jump_ranges.append((jstart,jend))
        return jump_ranges


    def rearrange_jumps(self, jump_ranges):
        phase_wrap_indices = {}
        for pw in self.phase_wraps:
            pw_pos = self.get_nTOAs(pw[0]) + pw[1]
            phase_wrap_indices[pw_pos] = self.phase_wraps[pw]

        all_TOAs = []
        for tset in self.TOAsets:
            for toa in tset.TOAs:
                all_TOAs.append(toa)

        new_TOAsets = [TOAset(jumped=False)]
        for ii,toa in enumerate(all_TOAs):
            if any([r[0] == ii for r in jump_ranges]):
                if new_TOAsets[-1].get_nTOAs():
                    new_TOAsets.append(TOAset(jumped=True))
                else:
                    new_TOAsets[-1].jumped = True
            new_TOAsets[-1].append_TOA(toa)
            if any([r[1] == ii for r in jump_ranges]):
                new_TOAsets.append(TOAset(jumped=False))
        if new_TOAsets[-1].get_nTOAs():
            self.TOAsets = new_TOAsets
        else:
            self.TOAsets = new_TOAsets[:-1]

        self.phase_wraps = {}
        for index in phase_wrap_indices:
            self.phase_wraps[self.get_position_of_TOA(index)] = \
              phase_wrap_indices[index]

    def to_tim_file(self, fname=None, toa_format=None):
        """
        If filename not provided, formatted TOAs are printed to screen.

        toa_format should be 'princeton', 'parkes', 'Tempo2', or 'ITOA'--if left
        as None, format will be self.default_format
        """
        if toa_format is None:
            toa_format = self.default_format
        lines = []
        if toa_format.lower() == 'tempo2':
            lines.append("FORMAT 1")
        if self.mode is not None:
            lines.append("MODE %d" % self.mode)
        if self.track is not None:
            lines.append("TRACK %d" % self.track)
        for ii in range(self.get_nTOAsets()):
            if self.TOAsets[ii].jumped:
                lines.append("JUMP")
            for jj in range(self.TOAsets[ii].get_nTOAs()):
                if (ii,jj) in self.phase_wraps:
                    ph_arg = "%+d" % self.phase_wraps[(ii,jj)] 
                    lines.append("PHASE %s" % ph_arg)
                if toa_format.lower() == 'parkes':
                    lines.append(self.TOAsets[ii].TOAs[jj].to_parkes_format())
                elif toa_format.lower() == 'itoa':
                    lines.append(self.TOAsets[ii].TOAs[jj].to_ITOA_format())
                elif toa_format.lower() == 'princeton':
                    lines.append(self.TOAsets[ii].TOAs[jj].to_princeton_format())
                elif toa_format.lower() == 'tempo2':
                    lines.append(self.TOAsets[ii].TOAs[jj].to_Tempo2_format())
                else:
                    print "TOA Format must be 'princeton', 'parkes', 'Tempo2', or 'ITOA'."
            if self.TOAsets[ii].jumped:
                lines.append("JUMP")

        if fname is None:
            for line in lines:
                print line
        else:
            with open(fname, 'w') as f:
                for line in lines:
                    f.write(line + "\n")

