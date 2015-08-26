#!/usr/bin/env python

# Copied from Patrick Lazarus's pyplotres 2015 Aug 18,
# then mangled by Erik Madsen (et al.) into its present form

import optparse
import sys
import re
import types
import subprocess
from shutil import copyfile

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
from matplotlib.widgets import SpanSelector
import numpy as np

import pyslalib.slalib as slalib
import binary_psr
import parfile as par
import residuals

from tempy_utils import un2str, find_freq_clusters, colors
from tempy_utils import CheckButtons
from tempy_utils import tempy_io
from tempy_utils.tempo_info import TempoHistory

# Available x-axis types
xvals = ['mjd', 'year', 'numtoa', 'orbitphase']
xind = 0
# Available y-axis types
yvals = ['phase', 'usec', 'sec']
yind = 0


class TempoResults:
    def __init__(self, freqbands):
        """Read TEMPO results (resid2.tmp, tempo.lis, timfile and parfiles)
            freqbands is a list of frequency pairs to display.
        """
        # Open tempo.lis. Parse it and find input .tim and .par files.
        # Also find output .par file.
        inputfiles_re = re.compile(r"Input data from (.*\.tim.*),  Parameters from (.*\.par.*)")
        outputfile_re = re.compile(r"Assumed parameters -- PSR (.*)$")
        tempolisfile = open("tempo.lis")
        intimfn, inparfn, outparfn = None, None, None
        for line in tempolisfile:
            match = inputfiles_re.search(line)
            if match:
                intimfn = match.group(1).strip()
                inparfn = match.group(2).strip()
            else:
                match = outputfile_re.search(line)
                if match:
                    outparfn = "%s.par" % match.group(1).strip()
            if (intimfn != None) and (inparfn != None) and (outparfn != None):
                # Found what we're looking for no need to continue parsing the file
                break
        tempolisfile.close()

        self.phase_wraps = {}
        self.jump_ranges = []

        tim = tempy_io.TOAset.from_tim_file(intimfn)
        tim_ordered_index = np.argsort(tim.TOAs)

        # if there are phase wraps in the tim file, we want to include those
        # in the intial plot
        for tim_wrap_index in tim.phase_wraps:
            wrap_index = tim_ordered_index[tim_wrap_index]
            self.phase_wraps[wrap_index] = \
              tim.phase_wraps[tim_wrap_index]

        # if there are jumps in the tim file, we want to include those in the
        # initial plot
        for tim_jstart,tim_jend in tim.jump_ranges:
            jstart = tim_ordered_index[tim_jstart]
            jend = tim_ordered_index[tim_jend]
            self.jump_ranges.append((jstart, jend))

        # Record filename
        self.inparfn = inparfn
        self.outparfn = outparfn
        self.intimfn = intimfn

        # Read parfiles
        self.inpar = par.psr_par(inparfn)
        self.outpar = par.psr_par(outparfn)

        # Read residuals
        r = residuals.read_residuals()

        self.max_TOA = r.bary_TOA.max()
        self.min_TOA = r.bary_TOA.min()

        ordered_index = np.argsort(r.bary_TOA)
        self.ordered_MJDs = r.bary_TOA[ordered_index]

        if freqbands is None:
            self.freqbands = find_freq_clusters(r.bary_freq)
        else:
            self.freqbands = freqbands
        self.residuals = {}
        for lo,hi in self.freqbands:
            indices = (r.bary_freq>=lo) & (r.bary_freq<hi)
            self.residuals[get_freq_label(lo, hi)] = \
                 Resids(r.bary_TOA[indices], r.bary_freq[indices],
                        #np.arange(r.numTOAs)[indices],
                        ordered_index[indices],
                        r.orbit_phs[indices],
                        r.postfit_phs[indices], r.postfit_sec[indices],
                        r.prefit_phs[indices], r.prefit_sec[indices],
                        r.uncertainty[indices], r.weight[indices],
                        self.inpar, self.outpar)

    def get_info(self, freq_label, index, postfit=True):
        """Given a freq_label and index return formatted text
            describing the TOA residual.
            Assume postfit period for calculating residual in phase,
            unless otherwise indicated.
        """
        r = self.residuals[freq_label]
        description = []
        description.append("TOA Selected:")
        description.append("\tNumber: %s" % r.TOA_index[index][0])
        description.append("\tEpoch (MJD): %s" % r.bary_TOA[index][0])
        if yvals[yind] == "phase":
            description.append("\tPre-fit residual (phase): %s" % \
              r.prefit_phs[index][0])
            description.append("\tPost-fit residual (phase): %s" % \
              r.postfit_phs[index][0])
            if postfit:
                description.append("\tUncertainty (phase): %s" % \
                  (r.uncertainty[index][0]/r.outpar.P0))
            else:
                description.append("\tUncertainty (phase): %s" % \
                  (r.uncertainty[index][0]/r.inpar.P0))
        elif yvals[yind] == "usec":
            description.append("\tPre-fit residual (usec): %s" % \
              (r.prefit_sec[index][0]*1e6))
            description.append("\tPost-fit residual (usec): %s" % \
              (r.postfit_sec[index][0]*1e6))
            description.append("\tUncertainty (usec): %s" % \
              (r.uncertainty[index][0]*1e6))
        elif yvals[yind] == "sec":
            description.append("\tPre-fit residual (sec): %s" % \
              r.prefit_sec[index][0])
            description.append("\tPost-fit residual (sec): %s" % \
              r.postfit_sec[index][0])
            description.append("\tUncertainty (sec): %s" % \
              r.uncertainty[index][0])
        description.append("\tFrequency (MHz): %s" % r.bary_freq[index][0])
        return description



class Resids:
    """The Resids object contains the following information
        about TEMPO residuals:
            bary_TOA
            bary_freq
            numTOAs
            orbit_phs
            postfit_phs
            postfit_sec
            prefit_phs
            prefit_sec
            uncertainty
            weight
    """
    def __init__(self, bary_TOA, bary_freq, TOA_index,
                 orbit_phs, postfit_phs, postfit_sec, prefit_phs, prefit_sec,
                 uncertainty, weight, inpar, outpar):
        self.bary_TOA = np.array(bary_TOA)
        self.bary_freq = np.array(bary_freq)
        self.TOA_index = np.array(TOA_index)
        self.orbit_phs = np.array(orbit_phs)
        self.postfit_phs = np.array(postfit_phs)
        self.postfit_sec = np.array(postfit_sec)
        self.prefit_phs = np.array(prefit_phs)
        self.prefit_sec = np.array(prefit_sec)
        self.uncertainty = np.array(uncertainty)
        self.weight = np.array(weight)
        self.inpar = inpar
        self.outpar = outpar


    def get_xdata(self, key):
        """Return label describing xaxis and the corresponding
            data given keyword 'key'.
        """
        if not isinstance(key, types.StringType):
            raise ValueError("key must be of type string.")
        xopt = key.lower()
        if xopt == 'numtoa':
            xdata = self.TOA_index.copy()
            xlabel = "TOA Number"
        elif xopt == 'mjd':
            xdata = self.bary_TOA.copy()
            xlabel = "MJD"
        elif xopt == 'orbitphase':
            xdata = self.orbit_phs.copy()
            xlabel = "Orbital Phase"
        elif xopt == 'year':
            xdata = mjd_to_year(self.bary_TOA)
            xlabel = "Year"
        else:
            raise ValueError("Unknown xaxis type (%s)." % xopt)
        return (xlabel, xdata)


    def get_ydata(self, key, postfit=True, phase_wraps={}):
        """Return label describing yaxis and the corresponding
            data/errors given keyword 'key'.
            'postfit' is a boolean argument that determines if
            postfit, or prefit data is to be returned.
        """
        if postfit==True:
            prefix='Postfit '
        if postfit==False:
            prefix='Prefit '
        if not isinstance(key, types.StringType):
            raise ValueError("key must be of type string.")
        yopt = key.lower()
        if postfit:
            if yopt == 'phase':
                ydata = self.postfit_phs.copy()
                #
                # NOTE: Should use P at TOA not at PEPOCH
                #
                yerror = self.uncertainty/self.outpar.P0
                ylabel = prefix+"Residuals (Phase)"
            elif yopt == 'usec':
                ydata = self.postfit_sec*1e6
                yerror = self.uncertainty*1e6
                ylabel = prefix+"Residuals (uSeconds)"
            elif yopt == 'sec':
                ydata = self.postfit_sec.copy()
                yerror = self.uncertainty.copy()
                ylabel = prefix+"Residuals (Seconds)"
            else:
                raise ValueError("Unknown yaxis type (%s)." % yopt)
        else:
            if yopt=='phase':
                ydata = self.prefit_phs.copy()
                #
                # NOTE: Should use P at TOA not at PEPOCH
                #
                yerror = self.uncertainty/self.inpar.P0
                ylabel = prefix+"Residuals (Phase)"
            elif yopt=='usec':
                ydata = self.prefit_sec*1e6
                yerror = self.uncertainty*1e6
                ylabel = prefix+"Residuals (uSeconds)"
            elif yopt=='sec':
                ydata = self.prefit_sec.copy()
                yerror = self.uncertainty.copy()
                ylabel = prefix+"Residuals (Seconds)"
            else:
                raise ValueError("Unknown yaxis type (%s)." % yopt)

        if postfit:
            for wrap_index in phase_wraps:
                if yopt=='phase':
                    ydata[self.TOA_index >= wrap_index] += \
                      phase_wraps[wrap_index]
                elif yopt=='usec':
                    ydata[self.TOA_index >= wrap_index] += \
                      phase_wraps[wrap_index]*self.outpar.P0*1e6
                elif yopt=='sec':
                    ydata[self.TOA_index >= wrap_index] += \
                      phase_wraps[wrap_index]*self.outpar.P0

        return (ylabel, ydata, yerror)



def plot_data(tempo_results, xkey, ykey,
              interactive=True, mark_peri=False, show_legend=True):
    subplot = 1
    numsubplots = 2
    global axes
    axes = []
    global ax_types
    ax_types = []
    global ax_phase_wraps
    ax_phase_wraps = []
    global ax_jump_ranges
    ax_jump_ranges = []
    handles = []
    labels = []

    for usepostfit in [False, True]:# Always use pre, then post
        TOAcount = 0
        # All subplots are in a single column
        if subplot == 1:
            axes.append(plt.subplot(numsubplots, 1, subplot))
        else:
            axes.append(plt.subplot(numsubplots, 1, subplot, sharex=axes[0]))

        if usepostfit:
            ax_types.append('post')
        else:
            ax_types.append('pre')

        ax_phase_wraps.append([])
        ax_jump_ranges.append([])

        # set tick formatter to not use scientific notation or an offset
        tick_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        tick_formatter.set_scientific(False)
        axes[-1].xaxis.set_major_formatter(tick_formatter)

        xmin, xmax = axes[0].get_xlim()

        for ii,(lo,hi) in enumerate(tempo_results.freqbands):
            freq_label = get_freq_label(lo, hi)
            resids = tempo_results.residuals[freq_label]

            xlabel, xdata = resids.get_xdata(xkey)
            ylabel, ydata, yerr = resids.get_ydata(ykey, usepostfit,
                                                   tempo_results.phase_wraps)
            if len(xdata):
                # Plot the residuals
                handle = plt.errorbar(xdata, ydata, yerr=yerr, fmt='.', \
                                      label=freq_label, picker=5,
                                      c=colors[len(tempo_results.freqbands)][ii])

                if subplot == 1:
                    handles.append(handle[0])
                    labels.append(freq_label)
                TOAcount += xdata.size

        # Plot phase wraps
        text_offset = offset_copy(axes[-1].transData, x=5, y=-10,
                                  units='dots')
        if usepostfit:
            pw = tempo_results.phase_wraps
        elif tempo_history.current_index > 0:
            pw = tempo_history.tempo_results[tempo_history.current_index-1].phase_wraps
        else:
            pw = []
        for wrap_index in pw:
            wrap_mjd_hi = tempo_results.ordered_MJDs[wrap_index]
            if wrap_index > 0:
                wrap_mjd_lo = tempo_results.ordered_MJDs[wrap_index-1]
            else:
                wrap_mjd_lo = wrap_mjd_hi
            if xkey == 'mjd':
                wrap_x = 0.5*(wrap_mjd_hi + wrap_mjd_lo)
            elif xkey == 'year':
                wrap_x = mjd_to_year(0.5*(wrap_mjd_hi + wrap_mjd_lo))
            elif xkey == 'numtoa':
                wrap_x = wrap_index - 0.5
            else:
                break
            wrap_color = {'pre':'pink', 'post':'red'}
            wrp = plt.axvline(wrap_x, ls=':', label='_nolegend_',
                              color=wrap_color[ax_types[-1]], lw=1.5)
            wrp_txt = plt.text(wrap_x, axes[-1].get_ylim()[1],
                               "%+d" % pw[wrap_index],
                               transform=text_offset, size='x-small',
                               color=wrap_color[ax_types[-1]])
            ax_phase_wraps[-1].append([wrp, wrp_txt])

        # set up span selector for setting new jump ranges
        options.jump_spans[ax_types[-1]] = SpanSelector(axes[-1], select_jump_range, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='orange'))
        options.jump_spans[ax_types[-1]].visible = options.jump_mode

        if subplot > 1:
            axes[0].set_xlim((xmin, xmax))

        # Finish off the plot
        plt.axhline(0, ls='--', label="_nolegend_", c='k', lw=0.5)
        axes[-1].ticklabel_format(style='plain', axis='x')

        if mark_peri and hasattr(tempo_results.outpar, 'BINARY'):
            # Be sure to check if pulsar is in a binary
            # Cannot mark passage of periastron if not a binary
            if usepostfit:
                binpsr = binary_psr.binary_psr(tempo_results.outpar.FILE)
            else:
                binpsr = binary_psr.binary_psr(tempo_results.inpar.FILE)
            xmin, xmax = axes[0].get_xlim()
            mjd_min = tempo_results.min_TOA
            mjd_max = tempo_results.max_TOA
            guess_mjds = np.arange(mjd_max + binpsr.par.PB, \
                                mjd_min - binpsr.par.PB, -binpsr.par.PB)
            for mjd in guess_mjds:
                peri_mjd = binpsr.most_recent_peri(float(mjd))
                if xkey == 'mjd':
                    plt.axvline(peri_mjd, ls=':', label='_nolegend_', c='k', lw=0.5)
                elif xkey == 'year':
                    print "plotting peri passage"
                    plt.axvline(mjd_to_year(peri_mjd), ls=':', label='_nolegend_', c='k', lw=0.5)
            axes[0].set_xlim((xmin, xmax))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        subplot += 1

    # Plot jump ranges
    for jstart,jend in tempo_results.jump_ranges:
        plot_jump_range((jstart,jend))
#        axes[-1].set_ylim((ymin, ymax))


    if numsubplots > 1:
        # Increase spacing between subplots.
        plt.subplots_adjust(hspace=0.25)

    # Write name of input files used for timing on figure
    if interactive:
        fntext = "Solution %d of %d, TOA file: %s, Parameter file: %s" % \
          (tempo_history.current_index+1, tempo_history.get_nsolutions(),
           tempo_results.intimfn, tempo_results.inparfn)
        figure_text = plt.figtext(0.01, 0.01, fntext,
                                  verticalalignment='bottom',
                                  horizontalalignment='left')

 # Make the legend and set its visibility state
    leg = plt.figlegend(handles, labels, 'upper right')
    leg.set_visible(show_legend)
    leg.legendPatch.set_alpha(0.5)
    axes[0].xaxis.tick_top()
    plt.setp(axes[0].get_yticklabels()[0], visible=False)
    plt.setp(axes[1].get_yticklabels()[-1], visible=False)
    plt.setp(axes[0].get_yticklabels()[-1], visible=False)
    plt.setp(axes[1].get_yticklabels()[0], visible=False)
    plt.subplots_adjust(wspace=0.05, hspace = 0.0, left=0.15, bottom=0.1, right=0.8, top=0.9)



    nms=[]
    fitmes=[]
    for b in tempo_history.get_parfile():
        if tempo_history.get_parfile()[b].fit==None:
            tempo_history.get_parfile()[b].fit=0
        if not any(b in s for s in tempy_io.no_fit_pars):
            nms.append(b)
            fitmes.append(tempo_history.get_parfile()[b].fit)
    rax = plt.axes([0.85, 0.1, 0.1, 0.8])
    rax.set_frame_on(False)
    options.fitcheck = CheckButtons(rax, nms, fitmes)
    options.fitcheck.on_clicked(update_fit_flag)
    redrawplot()


def update_fit_flag(label, button):
    if button=='left':
        if label:
            if tempo_history.get_parfile()[label].fit==None:
                tempo_history.get_parfile()[label].fit=0
            tempo_history.get_parfile()[label].fit=np.str((np.int(tempo_history.get_parfile()[label].fit)+1)%2)
    if button=='right':
        if label in ['RAJ', 'DECJ']:
            newvalue = raw_input("New value of %s [%s]: " % (label,  tempo_history.get_parfile()[label].value))
            if not newvalue: newvalue = "%s" % tempo_history.get_parfile()[label].value
            tempo_history.get_parfile()[label].value=newvalue
        else:
            newvalue = raw_input("New value of %s [%1.9e]: " % (label,  tempo_history.get_parfile()[label].value))
            if not newvalue: newvalue = "%1.9e" % tempo_history.get_parfile()[label].value
            tempo_history.get_parfile()[label].value=np.float(newvalue)




def create_plot():
    # Set up the plot
    fig = plt.figure(figsize=(11,8.5))
    # Force user to interact via custom inputs
    fig.canvas.toolbar.pack_forget()
    fig.set_facecolor("white")

def get_freq_label(lo, hi):
    """Return frequency label given a lo and hi
        frequency pair.
    """
    if hi is 'inf':
        return "%.0f - Inf MHz" % (lo)
    else:
        return "%.0f - %.0f MHz" % (lo, hi)


def savefigure(savefn='./resid2.tmp.ps'):
    print "Saving plot to %s" % savefn
    plt.savefig(savefn, orientation='landscape', papertype='letter')

def reloadplot(with_tempo_results=None):
    global tempo_results
    global options
    # Reload residuals and replot
    print "Plotting..."
    fig = plt.gcf()
    fig.set_visible(False)
    plt.clf() # clear figure
    if with_tempo_results is None:
        tempo_results = TempoResults(options.freqbands)
    else:
        tempo_results = with_tempo_results
    try:
        plot_data(tempo_results, options.xaxis, options.yaxis,
                  interactive=options.interactive,
                  mark_peri=options.mark_peri, show_legend=options.legend)
    except EmptyPlotValueError, msg:
        print msg
    fig.set_visible(True)
    redrawplot()


def redrawplot():
    #plt.show is keeping the plot open on nimrod, as opposed to plt.draw
    plt.draw()
    #plt.show()

def quit():
    print "Quitting..."
    sys.exit(0)

def pick(event):
    global tempo_results
    global options
    global ax_jump_ranges
#    if event.mouseevent.button == 1:
#        index = event.ind
#        axes = event.mouseevent.inaxes
#        if axes:
#            title = axes.get_title()
#            postfit = ("Postfit" in title)
#        if len(index) == 1:
#            freq_label = event.artist.get_label()
#            info = tempo_results.get_info(freq_label, index, postfit)
#            print_text(info)
#        else:
#            print "Multiple TOAs selected. Zoom in and try again."
    if event.mouseevent.button == 3:
        if options.jump_mode:
            for ax in ax_jump_ranges:
                if event.artist in ax:
                    xmin = event.artist.get_paths()[0].get_extents().xmin
                    xmax = event.artist.get_paths()[0].get_extents().xmax
                    delete_jump_range(0.5*(xmin+xmax))
                    redrawplot()
        else:
            print "Must be in jump edit mode ('j') to delete jump ranges."


def print_text(lines, *args, **kwargs):
    """Print lines of text (in a list) in the terminal."""
    print '\n'.join(lines)


def print_help():
    # Display help
    print "Helping..."
    print "-"*80
    print "Help - Hotkeys definitions:"
    print "\th - Display this help"
    print "\tq - Quit"
    print "\ts - Save current plot(s) to PostScript file"
    print "\tc - Try to determine optimal color pallete"
    print "\tz - Toggle Zoom-mode on/off"
    print "\tm - Toggle marking of periastron passages on/off"
    print "\tL - Toggle legend on/off"
    print "\tj - Toggle jump edit mode (left click to start/end jump range,\n"\
          "\t    right click to remove jump range)"
    print "\t+ - Insert positive phase wrap at cursor position"
    print "\t- - Insert negative phase wrap at cursor position"
    print "\t[Backspace] - Remove all phase wraps"
    print "\tx - Run Tempo with current postfit parameters and phase wraps"
    print "\tb - Return to previous Tempo solution"
    print "\tn - Go to next Tempo solution"
    print "\td - Dump current Tempo solution to new par/tim files"
    print "\tu - Go to original view (unzoom)"
    print "\t< - Go to previous view"
    print "\t> - Go to next view"
    print "\tr - Reload residuals"
    print "\tt - Cycle through y-axis types ('phase', 'usec', 'sec')"
    print "\t[Space] - Cycle through x-axis types ('mjd', 'year', 'numtoa',\n"\
          "\t          'orbitphase')"
    print "\t[Left mouse] - Select TOA (display info in terminal)"
    print "\t             - Select zoom region (if Zoom-mode is on)"
    print "-"*80

def run_tempo():
    global tempo_results
    global tempo_history
    par_fname = tempo_results.outpar.FILE
    tempo_history.save_outpar(par_fname)
    if par_fname.split('.')[-1] == 'tempy':
        new_par = par_fname
    else:
        new_par = par_fname + '.tempy'
    copyfile(tempo_results.outpar.FILE, new_par)
    tim = tempy_io.TOAset.from_tim_file(tempo_results.intimfn)
    tim.phase_wraps = {}
    tim.jump_ranges = []
    tempy_io.write_parfile(tempo_history.get_parfile(), new_par)

    tim_ordered_index = np.argsort(tim.TOAs)
    for wrap_index in tempo_results.phase_wraps:
        tim_wrap_index = np.where(tim_ordered_index == wrap_index)[0][0]
        tim.phase_wraps[tim_wrap_index] = tempo_results.phase_wraps[wrap_index]
    for jstart,jend in tempo_results.jump_ranges:
        tim_jstart = np.where(tim_ordered_index == jstart)[0][0]
        tim_jend = np.where(tim_ordered_index == jend)[0][0]
        tim.jump_ranges.append((tim_jstart,tim_jend))
    if tempo_results.intimfn.split('.')[-1] == 'tempy':
        new_timfn = tempo_results.intimfn
    else:
        new_timfn = tempo_results.intimfn + ".tempy"
    tim.to_tim_file(new_timfn)
    subprocess.call(["tempo", "-f", new_par, new_timfn])


def increment_phase_wrap(xdata, phase_offset):
    global tempo_results
    global options
    if options.xaxis == 'mjd':
        where_wrap = np.searchsorted(tempo_results.ordered_MJDs, xdata)
    elif options.xaxis == 'year':
        all_years = mjd_to_year(tempo_results.ordered_MJDs)
        where_wrap = np.searchsorted(all_years, xdata)
    elif options.xaxis == 'numtoa':
        where_wrap = int(np.ceil(xdata))
    else:
        return
    if where_wrap >= len(tempo_results.ordered_MJDs):
        return
    if where_wrap in tempo_results.phase_wraps:
        tempo_results.phase_wraps[where_wrap] += phase_offset
        if tempo_results.phase_wraps[where_wrap] == 0:
            del tempo_results.phase_wraps[where_wrap]
    else:
        tempo_results.phase_wraps[where_wrap] = phase_offset
    tempo_history.set_tempo_results(tempo_results)

def is_in_jump_range(index):
    """
    Returns which jump range contains given TOA index, or None if
    this index is not in a jump range
    """
    global tempo_results
    for ii,(jstart,jend) in enumerate(tempo_results.jump_ranges):
        if index >= jstart and index <= jend:
            return ii
    return None

def jump_ranges_between(index1, index2):
    """
    Returns all jump ranges that exist between index1 and index2 (even if they
    only overlap partly)
    """
    global tempo_results
    which_jump_ranges = []
    for ii,(jstart,jend) in enumerate(tempo_results.jump_ranges):
        if (jend >= index1 and jend <= index2) or \
          (jstart <= index2 and jstart >= index1):
            which_jump_ranges.append(ii)
    return which_jump_ranges

def select_jump_range(xdata_min, xdata_max):
    global tempo_results
    global options
    global ax_types
    if options.xaxis == 'mjd':
        xmin, xmax = np.searchsorted(tempo_results.ordered_MJDs,
                                     [xdata_min, xdata_max])
    elif options.xaxis == 'year':
        all_years = mjd_to_year(tempo_results.ordered_MJDs)
        xmin, xmax = np.searchsorted(all_years, [xdata_min, xdata_max])
    elif options.xaxis == 'numtoa':
        xmin = int(np.ceil(xdata_min))
        xmax = int(np.floor(xdata_max))
    if xmin >= xmax:
        return
    xmax -= 1

    xmin_in_jump_range = is_in_jump_range(xmin)
    xmax_in_jump_range = is_in_jump_range(xmax-1)
    if xmin_in_jump_range is not None and xmax_in_jump_range is not None:
        if xmin_in_jump_range == xmax_in_jump_range:
            #del tempo_results.jump_ranges[xmin_in_jump_range]
            delete_jump_range_index(xmin_in_jump_range)
    which_jump_ranges = jump_ranges_between(xmin, xmax)
    for ii in reversed(sorted(which_jump_ranges)):
        #del tempo_results.jump_ranges[ii]
        delete_jump_range_index(ii)
    tempo_results.jump_ranges.append((xmin,xmax))
    plot_jump_range((xmin,xmax))
    redrawplot()

def plot_jump_range(jump_range):
    global tempo_results
    global options
    global ax_jump_ranges
    global axes
    #ax = axes[ax_types.index(ax_type)]
    jstart,jend = jump_range
    jstart_mjd = tempo_results.ordered_MJDs[jstart]
    jend_mjd = tempo_results.ordered_MJDs[jend]
    extend_frac = 0.1
    if jstart > 0:
        jstart_mjd_before = tempo_results.ordered_MJDs[jstart-1]
    else:
        jstart_mjd_before = jstart_mjd
    if jend < len(tempo_results.ordered_MJDs)-1:
        jend_mjd_after = tempo_results.ordered_MJDs[jend+1]
    else:
        jend_mjd_after = jend_mjd
    dist_before = jstart_mjd - \
      ((1-extend_frac)*jstart_mjd + extend_frac*jstart_mjd_before)
    dist_after = (extend_frac*jend_mjd_after + \
      (1-extend_frac)*jend_mjd) - jend_mjd
    dist_mjd = min(dist_before, dist_after)
    if dist_mjd < 1e-8:
        dist_mjd = max(dist_before, dist_after)
    jstart_mjd -= dist_mjd
    jend_mjd += dist_mjd
    if options.xaxis == 'mjd':
        jstart_x = jstart_mjd
        jend_x = jend_mjd
    elif options.xaxis == 'year':
        jstart_x = mjd_to_year(jstart_mjd)[0]
        jend_x = mjd_to_year(jend_mjd)[0]
    elif options.xaxis == 'numtoa':
        jstart_x = jstart + extend_frac
        jend_x = jend - extend_frac
    else:
        return
    for ii,ax in enumerate(axes):
        ymin,ymax = ax.get_ylim()
        jmp = ax.fill_betweenx((ymin*100, ymax*100), jstart_x, jend_x,
                               edgecolor="orange", facecolor='yellow',
                               lw=0.5, alpha=0.3, picker=True)
        ax_jump_ranges[ii].append(jmp)
        ax.set_ylim(ymin,ymax)
    #redrawplot()


def delete_jump_range_index(index):
    global tempo_results
    global ax_jump_ranges
    del tempo_results.jump_ranges[index]
    for ax in ax_jump_ranges:
        ax.pop(index).remove()

def delete_jump_range(xdata):
    global tempo_results
    global options
    global ax_jump_ranges
    if options.xaxis == 'mjd':
        where_clicked = np.searchsorted(tempo_results.ordered_MJDs, xdata)
    elif options.xaxis == 'year':
        all_years = mjd_to_year(tempo_results.ordered_MJDs)
        where_clicked = np.searchsorted(all_years, xdata)
    elif options.xaxis == 'numtoa':
        where_clicked = xdata
    else:
        return
    to_delete = is_in_jump_range(where_clicked)
    if to_delete is not None:
        delete_jump_range_index(to_delete)

def keypress(event):
    global tempo_results
    global tempo_history
    global options
    global xind, xvals
    global yind, yvals
    if type(event.key) in [types.StringType, types.UnicodeType]:
        if event.key.lower() == 'q':
            quit()
        elif event.key.lower() == 's':
            savefigure()
        elif event.key.lower() == 'c':
            options.freqbands = None
            reloadplot()
        elif event.key == 'r':
            reloadplot()
        elif event.key.upper() == 'L':
            leg = plt.gcf().legends[0]
            options.legend = not options.legend
            leg.set_visible(options.legend)
            redrawplot()
        elif event.key.lower() == 'z':
            # Turn on zoom mode
            print "Toggling zoom mode..."
            options.jump_mode = False
            for k in options.jump_spans:
                options.jump_spans[k].visible = False
            event.canvas.toolbar.zoom()
        elif event.key.lower() == 'm':
            # Toggle peri markings
            print "Toggling periastron passage markings..."
            options.mark_peri = not options.mark_peri
            reloadplot()
        elif event.key == '+' or event.key == '=':
            try:
                increment_phase_wrap(event.xdata, 1)
                reloadplot(tempo_results)
            except:
                pass
        elif event.key == '-' or event.key == '_':
            try:
                increment_phase_wrap(event.xdata, -1)
                reloadplot(tempo_results)
            except: pass
        elif event.key == "backspace":
            tempo_results.phase_wraps = {}
            reloadplot(tempo_results)
        elif event.key.lower() == 'j':
            if event.canvas.toolbar._active is not None:
                if event.canvas.toolbar._active.lower() == 'zoom':
                    event.canvas.toolbar.zoom()
            options.jump_mode = not options.jump_mode
            if options.jump_mode:
                print "Jump edit mode on"
            else:
                print "Jump edit mode off"
            for k in options.jump_spans:
                options.jump_spans[k].visible = not options.jump_spans[k].visible
        elif event.key.lower() == 'x':
            run_tempo()
            tempo_results = TempoResults(options.freqbands)
            tempo_history.append(tempo_results)
            tempo_history.print_formatted_pars()
            reloadplot()
        elif event.key.lower() == 'b':
            # Previous solution
            tempo_history.seek_prev_solution()
            tempo_results = tempo_history.get_tempo_results()
            reloadplot(tempo_results)
        elif event.key.lower() == 'n':
            # Next solution
            tempo_history.seek_next_solution()
            tempo_results = tempo_history.get_tempo_results()
            reloadplot(tempo_results)
        elif event.key == 'R':
            # First solution
            tempo_history.seek_first_solution()
            tempo_results = tempo_history.get_tempo_results()
            reloadplot(tempo_results)
        elif event.key.lower() == 'd':
            basename = "%s.tpy" % tempo_results.outpar.PSR
            par_fname = raw_input("Output parfile [%s.par]: " % basename)
            if not par_fname: par_fname = "%s.par" % basename
            tim_fname = raw_input("Output timfile [%s.tim]: " % basename)
            if not tim_fname: tim_fname = "%s.tim" % basename
            tempo_history.save_outpar(par_fname)
            tempo_history.save_timfile(tim_fname)
        elif event.key.lower() == 'u':
            # Restore plot to original view
            print "Restoring plot..."
            if event.canvas.toolbar._active is not None:
                if event.canvas.toolbar._active.lower() == 'zoom':
                    event.canvas.toolbar.zoom()
            event.canvas.toolbar.home()
        elif event.key.lower() == ',' or event.key.lower() == '<':
            # Go back to previous plot view
            print "Going back..."
            event.canvas.toolbar.back()
        elif event.key.lower() == '.' or event.key.lower() == '>':
            # Go forward to next plot view
            print "Going forward..."
            event.canvas.toolbar.forward()
        elif event.key.lower() == ' ':
            xind = (xind + 1) % len(xvals)
            print "Toggling plot type...[%s]"%xvals[xind], xind
            options.xaxis = xvals[xind]
            reloadplot(tempo_results)
        elif event.key.lower() == 't':
            yind = (yind + 1) % len(yvals)
            print "Toggling plot scale...[%s]"%yvals[yind], yind
            options.yaxis = yvals[yind]
            reloadplot(tempo_results)
        elif event.key.lower() == 'h':
            print_help()


def mjd_to_year(mjds):
    mjds = np.asarray(mjds)
    if mjds.size < 1:
        return mjds
    old_shape = mjds.shape # Remember original shape
    mjds.shape = (mjds.size, 1)
    years, months, days, fracs, stats = np.apply_along_axis(slalib.sla_djcl, 1, mjds).transpose()
    # Take into account leap years
    daysperyear = (((years % 4) == 0) & (((years % 100) != 0) | ((years % 400) == 0))) * 1 + 365.0
    years, days, stats = np.array([slalib.sla_clyd(*ymd) for ymd in np.vstack((years, months, days)).transpose()]).transpose()
    mjds.shape = old_shape # Change back to original shape
    return (years + (days + fracs) / daysperyear)


def parse_options():
    (options, other_args) = parser.parse_args()
    if not options.freqs:
        # Default frequency bands
        freqbands = [['0', '400'],
                     ['400', '600'],
                     ['600', '1000'],
                     ['1000', '1600'],
                     ['1600', '2400'],
                     ['2400', 'inf']]
    else:
        freqbands = []
        for fopt in options.freqs:
            f = fopt.split(':')
            if f[0]=='':
                f[0] = '0'
            if f[-1]=='':
                f[-1] = 'inf'
            if len(f) > 2:
                for i in range(0, len(f)-1):
                    freqbands.append(f[i:i+2])
            else:
                freqbands.append(f)
    freqbands = np.array(freqbands).astype(float)
    freqbands[freqbands.argsort(axis=0).transpose()[0]]
    if np.any(freqbands.flat != sorted(freqbands.flat)):
        raise ValueError("Frequency bands have overlaps or are inverted.")
    options.freqbands = freqbands

    if not options.prefit and not options.postfit:
        # If neither prefit or postfit are selected
        # show both
        options.postfit = True
        options.prefit = True

    if options.xaxis.lower() not in xvals:
        raise BadOptionValueError("Option to -x/--x-axis (%s) is not "\
          "permitted." % options.xaxis)
    if options.yaxis.lower() not in yvals:
        raise BadOptionValueError("Option to -y/--y-axis (%s) is not "\
          "permitted." % options.yaxis)

    options.jump_spans = {}
    options.jump_mode = False

    if options.initial_parfile and len(other_args):
        options.initial_timfile = other_args[-1]
        options.run_initial_fit = True
    else:
        options.run_initial_fit = False

    return options


def main():
    global tempo_results
    global tempo_history
    global options
    options = parse_options()

    if options.run_initial_fit:
        print "Running TEMPO with parfile %s and tim file %s" % \
          (options.initial_parfile, options.initial_timfile)
        subprocess.call(["tempo", "-f", options.initial_parfile,
                        options.initial_timfile])
    else:
        print "Initial par/tim files not provided, attempting to load " \
              "existing TEMPO results."

    tempo_results = TempoResults(options.freqbands)
    tempo_history = TempoHistory(tempo_results)

    tempo_history.print_formatted_pars()

    create_plot()
    reloadplot()

    if options.interactive:
        fig = plt.gcf() # current figure

        # Before setting up our own event handlers delete matplotlib's
        # default 'key_press_event' handler.
        defcids = fig.canvas.callbacks.callbacks['key_press_event'].keys()
        for cid in defcids:
            fig.canvas.callbacks.disconnect(cid)

        # Now, register our event callback functions
        cid_keypress = fig.canvas.mpl_connect('key_press_event', keypress)
        #cid_keypress = fig.canvas.mpl_connect('pick_event', click)
        cid_pick = fig.canvas.mpl_connect('pick_event', pick)

        # Finally, let the show begin!
        #plt.ion()
        plt.show()
    else:
        # Save figure and quit
        savefigure()
        quit()


class BadOptionValueError(ValueError):
    """Bad value passed to option parser.
    """
    pass


class EmptyPlotValueError(ValueError):
    """Empty plot.
    """
    pass


if __name__=='__main__':
    parser = optparse.OptionParser(prog="tempy.py")
    parser.add_option('-f', dest='initial_parfile', type='string', \
                        help="A TEMPO parfile for initial fit. If provided," \
                             " a tim file must appear as the final command"\
                             " line argument. If not provided, the last run of"\
                             " TEMPO in the current directory is used.",\
                        default="")
    parser.add_option('--freq', dest='freqs', action='append', \
                        help="Band of frequencies, in MHz, to be plotted " \
                             "(format xxx:yyy). Each band will have a " \
                             " different colour. Multiple --freq options " \
                             " are allowed. (Default: Plot all frequencies " \
                             "in single colour.)", \
                        default=[])
    parser.add_option('-x', '--x-axis', dest='xaxis', type='string', \
                        help="Values to plot on x-axis. Must be one of " \
                             "%s. (Default: '%s')" % (str(xvals), xvals[xind]),
                        default=xvals[xind])
    parser.add_option('-y', '--y-axis', dest='yaxis', type='string', \
                        help="Values to plot on y-axis. Must be one of "
                             "%s. (Default: '%s')" % (str(yvals), yvals[yind]), \
                        default=yvals[yind])
    parser.add_option('--post', dest='postfit', action='store_true', \
                        help="Show postfit residuals. (Default: Don't show " \
                             "postfit.)", \
                        default=False)
    parser.add_option('--pre', dest='prefit', action='store_true', \
                        help="Show prefit residuals. (Default: Don't show " \
                             "prefit.)", \
                        default=False)
    parser.add_option('-l', '--legend', dest='legend', action='store_true', \
                        help="Show legend of frequencies. (Default: Do not " \
                             "show legend.)", \
                        default=False)
    parser.add_option('--mark-peri', dest='mark_peri', action='store_true', \
                        help="Mark passage of periastron. (Default: don't " \
                             "mark periastron.)", \
                        default=False)
    parser.add_option('--non-interactive', dest='interactive', \
                        action='store_false', default=True, \
                        help="Save figure and exit. (Default: Show plot, " \
                             "only save if requested.)")
    main()
