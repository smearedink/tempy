from .. import tempy_io
from .. import un2str

class TempoHistory:
    def __init__(self, tempo_results=None):
        self.current_index = -1
        # parfiles are currently stored simply as raw strings
        self.inpars = []
        self.outpars = []
        # timfiles are TOAfile objects
        self.timfiles = []
        self.tempo_results = []
        if tempo_results is not None:
            self.append(tempo_results)

    def get_nsolutions(self):
        return len(self.tempo_results)

    def seek_next_solution(self):
        new_index = self.current_index + 1
        if new_index < self.get_nsolutions():
            self.current_index = new_index
            print "Moving ahead to solution %d of %d" % (new_index + 1,
                                                         self.get_nsolutions())
        else:
            print "Already at solution %d of %d" % (self.get_nsolutions(),
                                                    self.get_nsolutions())

    def seek_prev_solution(self):
        new_index = self.current_index - 1
        if new_index >= 0 and self.get_nsolutions():
            self.current_index = new_index
            print "Moving back to solution %d of %d" % (new_index + 1,
                                                        self.get_nsolutions())
        else:
            print "Already at solution 1 of %d" % (self.get_nsolutions())

    def seek_first_solution(self):
        if self.get_nsolutions():
            self.current_index = 0

    def seek_solution(self, n):
        if n >= 0 and n < self.get_nsolutions():
            self.current_index = n

    def clear_future_history(self):
        end = self.current_index + 1
        self.inpars = self.inpars[:end]
        self.outpars = self.outpars[:end]
        self.timfiles = self.timfiles[:end]
        self.tempo_results = self.tempo_results[:end]

    def append(self, tempo_results, increment_current=True):
        self.clear_future_history()
        with open(tempo_results.inparfn, 'r') as f:
            inpar = f.readlines()
            self.inpars.append(inpar)
        #with open(tempo_results.outparfn, 'r') as f:
        #    outpar = f.readlines()
        #    self.outpars.append(outpar)
        self.outpars.append(tempy_io.read_parfile(tempo_results.outpar.FILE))
        timfile = tempy_io.TOAfile.from_tim_file(tempo_results.intimfn)
        self.timfiles.append(timfile)
        self.tempo_results.append(tempo_results)
        if increment_current:
            self.current_index += 1

    def get_tempo_results(self, index=None):
        if index is None:
            index = self.current_index
        return self.tempo_results[index]

    def set_tempo_results(self, tempo_results, index=None):
        if index is None:
            index = self.current_index
        self.tempo_results[index] = tempo_results

    def get_parfile(self, index=None):
        if index is None:
            index = self.current_index
        return self.outpars[index]
        
    def get_timfile(self, index=None):
        if index is None:
            index = self.current_index
        return self.timfiles[index]

    def save_inpar(self, fname):
        with open(fname, 'w') as f:
            f.writelines(self.inpars[self.current_index])
        print "Wrote input parfile %s" % fname

    def save_outpar(self, fname):
        #with open(fname, 'w') as f:
        #    f.writelines(self.outpars[self.current_index])
        tempy_io.write_parfile(self.outpars[self.current_index], fname)
        print "Wrote output parfile %s" % fname

    def save_timfile(self, fname):
        self.timfiles[self.current_index].to_tim_file(fname)
        print "Wrote tim file %s" % fname

    def print_formatted_pars(self, index=None):
        if index is None:
            index = self.current_index
        no_disp_pars = list(tempy_io.no_fit_pars)
        for par in ['START', 'FINISH', 'PEPOCH']:
            if par in no_disp_pars:
                no_disp_pars.remove(par)
        formatted_par_line = "%20s: %1s %-18s"
        output_par = self.get_parfile(index)
        for par in output_par:
            if par not in no_disp_pars:
                if output_par[par].fit:
                    fit_str = '*'
                else:
                    fit_str = ''
                if output_par[par].error is None:
                    val = "%s" % output_par[par].value
                else:
                    try:
                        if par == 'RAJ' or par == 'DECJ':
                            split_str = output_par[par].value.split(':')
                            split_str[-1] = un2str(float(split_str[-1]),
                                                   output_par[par].error)
                            val = ''
                            for item in split_str:
                                val += item + ":"
                            val = val[:-1]
                        else:
                            val = un2str(output_par[par].value,
                                         output_par[par].error)
                    except:
                        err = output_par[par].error
                        val = "%s +/- %f" % (output_par[par].value, err)
                print formatted_par_line % (par, fit_str, val)
