#!/usr/bin/env python

# Copied from Patrick Lazarus's pyplotres 2015 Aug 18,
# then mangled by Erik Madsen into its present form

import optparse
import sys
import re
import os
import types
import warnings
import subprocess
from shutil import copyfile

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
import numpy as np

import pyslalib.slalib as slalib
import binary_psr
import parfile as par
import residuals

from scipy.cluster.vq import kmeans2

import toa

# Available x-axis types
xvals = ['mjd', 'year', 'numtoa', 'orbitphase']
xind = 0
# Available y-axis types
yvals = ['phase', 'usec', 'sec']
yind = 0

colors = {1: ['#000000'], # black
          2: ['#ff0000', '#0000ff'], # red blue
          3: ['#ff0000', '#008000', '#0000ff'], # red green blue
          4: ['#ff0000', '#FFA500', '#008000', '#0000ff'], # red orange green blue
          # red orange green blue violet
          5: ['#ff0000', '#FFA500', '#008000', '#0000ff', '#EE82EE'], 
          # red orange green blue indigo violet
          6: ['#ff0000', '#FFA500', '#008000', '#0000ff', '#4B0082', '#EE82EE'], 
          # red orange yellow green blue indigo violet
          7: ['#ff0000', '#FFA500', '#FFFF00', '#008000', '#0000ff', '#4B0082', '#EE82EE'], 
          # red orange yellow green blue indigo violet black
          8: ['#ff0000', '#FFA500', '#FFFF00', '#008000', '#0000ff', '#4B0082', '#EE82EE', '#000000']} 


def find_freq_clusters(freqs):
    # first make a histogram
    minf, maxf = freqs.min(), freqs.max()
    maxbins = 8  # related to the max colors defined...
    df = 4.0 # MHz
    if ((maxf - minf) < df):  # Only a single freq to our resolution
        return [[0.0, 'inf']]
    numbins = int((maxf - minf) / df) + 2
    lobound = minf - 0.5 * df
    hibound = lobound + numbins * df
    hist, edges = np.histogram(freqs, numbins, [lobound, hibound])
    # Now choose the maxbins biggest bins where there are TOAs
    hibins = hist.argsort()[::-1]
    hibins = hibins[hist[hibins] > 0]
    if len(hibins) > maxbins:
        hibins = hibins[:maxbins]
    ctrs = edges[hibins] + 0.5 * df
    ctrs.sort()
    # and use these as starting points for kmeans
    kmeans, indices = kmeans2(freqs, ctrs)
    if len(kmeans)==1:
        return [[0.0, 'inf']]
    elif len(kmeans)==2:
        return [[0.0, kmeans.mean()], [kmeans.mean(), 'inf']]
    else:
        freqbands = [[0.0, kmeans[0:2].mean()]]
        for ii in range(len(kmeans)-2):
            freqbands.append([kmeans[ii:ii+2].mean(), kmeans[ii+1:ii+3].mean()])
        freqbands.append([kmeans[-2:].mean(), 'inf'])
        return freqbands

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
                ylabel = "Residuals (Phase)"
            elif yopt == 'usec':
                ydata = self.postfit_sec*1e6
                yerror = self.uncertainty*1e6
                ylabel = "Residuals (uSeconds)"
            elif yopt == 'sec':
                ydata = self.postfit_sec.copy()
                yerror = self.uncertainty.copy()
                ylabel = "Residuals (Seconds)"
            else:
                raise ValueError("Unknown yaxis type (%s)." % yopt)
        else:
            if yopt=='phase':
                ydata = self.prefit_phs.copy()
                #
                # NOTE: Should use P at TOA not at PEPOCH
                #
                yerror = self.uncertainty/self.inpar.P0
                ylabel = "Residuals (Phase)"
            elif yopt=='usec':
                ydata = self.prefit_sec*1e6
                yerror = self.uncertainty*1e6
                ylabel = "Residuals (uSeconds)"
            elif yopt=='sec':
                ydata = self.prefit_sec.copy()
                yerror = self.uncertainty.copy()
                ylabel = "Residuals (Seconds)"
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


def plot_data(tempo_results, xkey, ykey, postfit=True, prefit=False,
              interactive=True, mark_peri=False, show_legend=True):
    # figure out what should be plotted
    # True means to plot postfit
    # False means to plot prefit
    if postfit and prefit:
        to_plot_postfit = [False, True]
    elif postfit and not prefit:
        to_plot_postfit = [True]
    elif not postfit and prefit:
        to_plot_postfit = [False]
    else:
        raise EmptyPlotValueError("At least one of prefit and postfit must be True.")
    subplot = 1
    numsubplots = len(to_plot_postfit)
    global axes
    axes = []
    handles = []
    labels = []

    for usepostfit in to_plot_postfit:
        TOAcount = 0
        # All subplots are in a single column
        if subplot == 1:
            axes.append(plt.subplot(numsubplots, 1, subplot))
        else:
            axes.append(plt.subplot(numsubplots, 1, subplot, sharex=axes[0]))

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
        for wrap_index in tempo_results.phase_wraps:
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
            wrap_color = ['pink', 'red'] # [prefit, postfit]
            plt.axvline(wrap_x, ls=':', label='_nolegend_',
                        color=wrap_color[usepostfit], lw=1.5)
            plt.text(wrap_x, axes[-1].get_ylim()[1],
                     "%+d" % tempo_results.phase_wraps[wrap_index],
                     transform=text_offset, size='x-small',
                     color=wrap_color[usepostfit])
        
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
        if usepostfit:
            plt.title("Postfit Residuals (Number of TOAs: %d)" % TOAcount)
        else:
            plt.title("Prefit Residuals (Number of TOAs: %d)" % TOAcount)
        subplot += 1

    if numsubplots > 1:
        # Increase spacing between subplots.
        plt.subplots_adjust(hspace=0.25)

    # Write name of input files used for timing on figure
    if interactive:
        fntext = "Solution %d of %d, TOA file: %s, Parameter file: %s" % \
          (tempo_history.current_index+1, tempo_history.get_nsolutions(),
           tempo_results.intimfn, tempo_results.inparfn)
        figure_text = plt.figtext(0.01, 0.01, fntext, verticalalignment='bottom', \
                            horizontalalignment='left')

    # Make the legend and set its visibility state
    leg = plt.figlegend(handles, labels, 'upper right')
    leg.set_visible(show_legend)
    leg.legendPatch.set_alpha(0.5)

def create_plot():
    # Set up the plot
    fig = plt.figure(figsize=(11,8.5))


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

def reloadplot(tempo_results=None):
    global options
    #global phase_wraps
    # Reload residuals and replot
    print "Plotting..."
    fig = plt.gcf()
    fig.set_visible(False)
    plt.clf() # clear figure
    if tempo_results is None:
        tempo_results = TempoResults(options.freqbands)
    try:
        plot_data(tempo_results, options.xaxis, options.yaxis,
                  postfit=options.postfit, prefit=options.prefit,
                  interactive=options.interactive,
                  mark_peri=options.mark_peri, show_legend=options.legend)
    except EmptyPlotValueError, msg:
        print msg
        print "Press 'p'/'P' to add prefit/postfit plot."
        plt.figtext(0.5, 0.5, (str(msg) + "\n" + \
                        "Press 'p'/'P' to add prefit/postfit plot."), \
                    horizontalalignment='center', \
                    verticalalignment='center', \
                    bbox=dict(facecolor='white', alpha=0.75))
    fig.set_visible(True)
    redrawplot()


def redrawplot():
    plt.draw()         #plt.show is keeping the plot open on nimrod, as opposed to plt.draw
    #plt.show()

def quit():
    print "Quitting..."
    sys.exit(0)


def pick(event):
    global tempo_results
    index = event.ind
    axes = event.mouseevent.inaxes
    if axes:
        title = axes.get_title()
        postfit = ("Postfit" in title)
    if len(index) == 1:
        freq_label = event.artist.get_label()
        info = tempo_results.get_info(freq_label, index, postfit)
        print_text(info)
    else:
        print "Multiple TOAs selected. Zoom in and try again."


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
    print "\tp - Toggle prefit display on/off"
    print "\tP - Toggle postfit display on/off"
    print "\tz - Toggle Zoom-mode on/off"
    print "\tm - Toggle marking of periastron passages on/off"
    print "\tL - Toggle legend on/off"
    print "\t+ - Insert positive phase wrap at cursor position"
    print "\t- - Insert negative phase wrap at cursor position"
    print "\t[Backspace] - Remove all phase wraps"
    print "\tT - Run Tempo with current postfit parameters and phase wraps"
    print "\tb - Return to previous Tempo solution"
    print "\tn - Go to next Tempo solution"
    print "\to - Go to original view"
    print "\t< - Go to previous view"
    print "\t> - Go to next view"
    print "\tx - Set x-axis limits (terminal input required)"
    print "\ty - Sey y-axis limits (terminal input required)"
    print "\tr - Reload residuals"
    print "\tt - Cycle through y-axis types ('phase', 'usec', 'sec')"
    print "\t[Space] - Cycle through x-axis types ('mjd', 'year', 'numtoa', 'orbitphase')"
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
    tim = toa.TOAset.from_princeton_file(tempo_results.intimfn)
    tim.phase_wraps = {}
    tim_ordered_index = np.argsort(tim.TOAs)
    for wrap_index in tempo_results.phase_wraps:
        tim_wrap_index = np.where(tim_ordered_index == wrap_index)[0][0]
        tim.phase_wraps[tim_wrap_index] = tempo_results.phase_wraps[wrap_index]
    if tempo_results.intimfn.split('.')[-1] == 'tempy':
        new_timfn = tempo_results.intimfn
    else:
        new_timfn = tempo_results.intimfn + ".tempy"
    tim.to_princeton_file(new_timfn)
    subprocess.call(["tempo", "-f", new_par, new_timfn])

class TempoHistory:
    def __init__(self, tempo_results=None):
        self.current_index = -1
        # parfiles are currently stored simply as raw strings
        self.inpars = []
        self.outpars = []
        # timfiles are TOAset objects
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
        with open(tempo_results.outparfn, 'r') as f:
            outpar = f.readlines()
            self.outpars.append(outpar)
        timfile = toa.TOAset.from_princeton_file(tempo_results.intimfn)
        self.timfiles.append(timfile)
        self.tempo_results.append(tempo_results)
        if increment_current:
            self.current_index += 1

    def get_current_tempo_results(self):
        return self.tempo_results[self.current_index]

    def save_inpar(self, fname):
        with open(fname, 'w') as f:
            f.writelines(self.inpars[self.current_index])
        print "Wrote input parfile %s" % fname

    def save_outpar(self, fname):
        with open(fname, 'w') as f:
            f.writelines(self.outpars[self.current_index])
        print "Wrote output parfile %s" % fname

    def save_timfile(self, fname):
        self.timfiles[self.current_index].to_princeton_file(fname)
        print "Wrote tim file %s" % fname

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
        elif event.key == 'T':
            run_tempo()
            tempo_results = TempoResults(options.freqbands)
            tempo_history.append(tempo_results)
            reloadplot()
        elif event.key.lower() == 'b':
            # Previous solution
            tempo_history.seek_prev_solution()
            tempo_results = tempo_history.get_current_tempo_results()
            reloadplot(tempo_results)
        elif event.key.lower() == 'n':
            # Next solution
            tempo_history.seek_next_solution()
            tempo_results = tempo_history.get_current_tempo_results()
            reloadplot(tempo_results)
        elif event.key == 'R':
            # First solution
            tempo_history.seek_first_solution()
            tempo_results = tempo_history.get_current_tempo_results()
            reloadplot(tempo_results)
        elif event.key.lower() == 'o':
            # Restore plot to original view
            print "Restoring plot..."
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
        elif event.key == 't':
            yind = (yind + 1) % len(yvals)
            print "Toggling plot scale...[%s]"%yvals[yind], yind
            options.yaxis = yvals[yind]
            reloadplot(tempo_results)
        elif event.key == 'p':
            options.prefit = not options.prefit
            print "Toggling prefit-residuals display to: %s" % \
                    ((options.prefit and "ON") or "OFF")
            reloadplot(tempo_results)
        elif event.key == 'P':
            options.postfit = not options.postfit
            print "Toggling postfit-residuals display to: %s" % \
                    ((options.postfit and "ON") or "OFF")
            reloadplot(tempo_results)
        elif event.key.lower() == 'x':
            # Set x-axis limits
            print "Setting x-axis limits. User input required..."
            xmin = raw_input("X-axis minimum: ")
            xmax = raw_input("X-axis maximum: ")
            try:
                xmin = float(xmin)
                xmax = float(xmax)
                if xmax <= xmin:
                    raise ValueError
            except ValueError:
                print "Bad values provided!"
                return
            plt.xlim(xmin, xmax)
        elif event.key.lower() == 'y':
            global axes
            # Set y-axis limits
            print "Setting y-axis limits. User input required..."
            if len(axes) == 2:
                axes_to_adjust = raw_input("Axes to adjust (pre/post): ")
                if axes_to_adjust.lower().startswith('pre'):
                    plt.axes(axes[0])
                elif axes_to_adjust.lower().startswith('post'):
                    plt.axes(axes[1])
                else:
                    raise ValueError
            ymin = raw_input("Y-axis minimum: ")
            ymax = raw_input("Y-axis maximum: ")
            try:
                ymin = float(ymin)
                ymax = float(ymax)
                if ymax <= ymin:
                    raise ValueError
            except ValueError:
                print "Bad values provided!"
                return
            plt.ylim(ymin, ymax)
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
    (options, sys.argv) = parser.parse_args()
    if sys.argv==[]:
        sys.argv = ['tempy.py']
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
    return options


def main():
    global tempo_results
    global tempo_history
    global options
    options = parse_options()
    tempo_results = TempoResults(options.freqbands)
    tempo_history = TempoHistory(tempo_results)

    tim = toa.TOAset.from_princeton_file(tempo_results.intimfn)
    tim_ordered_index = np.argsort(tim.TOAs)
    for tim_wrap_index in tim.phase_wraps:
        wrap_index = tim_ordered_index[tim_wrap_index]
        tempo_results.phase_wraps[wrap_index] = tim.phase_wraps[tim_wrap_index]
    
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
    parser.add_option('-f', '--freq', dest='freqs', action='append', \
                        help="Band of frequencies, in MHz, to be plotted " \
                             "(format xxx:yyy). Each band will have a " \
                             " different colour. Multiple -f/--freq options " \
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
