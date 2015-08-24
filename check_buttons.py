import numpy as np

from matplotlib.mlab import dist
from matplotlib.patches import Circle, Rectangle
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory
import pylab as plt
import toa
class Widget(object):
    """
    Abstract base class for GUI neutral widgets
    """
    drawon = True
    eventson = True

class CheckButtons(Widget):
    """
    A GUI neutral radio button

    The following attributes are exposed

     *ax*
        The :class:`matplotlib.axes.Axes` instance the buttons are
        located in

     *labels*
        List of :class:`matplotlib.text.Text` instances

     *lines*
        List of (line1, line2) tuples for the x's in the check boxes.
        These lines exist for each box, but have ``set_visible(False)``
        when its box is not checked.

     *rectangles*
        List of :class:`matplotlib.patches.Rectangle` instances

    Connect to the CheckButtons with the :meth:`on_clicked` method
    """
    def __init__(self, ax, labels, actives):
        """
        Add check buttons to :class:`matplotlib.axes.Axes` instance *ax*

        *labels*
            A len(buttons) list of labels as strings

        *actives*
            A len(buttons) list of booleans indicating whether
             the button is active
        """

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_navigate(False)

        if len(labels)>1:
            dy = 1./(len(labels)+1)
            ys = np.linspace(1-dy, dy, len(labels))
        else:
            dy = 0.25
            ys = [0.5]

        cnt = 0
        axcolor = ax.get_axis_bgcolor()
        self.labels = []
        self.lines = []
        self.rectangles = []

        lineparams = {'color':'k', 'linewidth':1.25, 'transform':ax.transAxes,
                      'solid_capstyle':'butt'}
        for y, label in zip(ys, labels):
            t = ax.text(0.25, y, label, transform=ax.transAxes,
                        horizontalalignment='left',
                        verticalalignment='center')

            w, h = dy/2., dy/2.
            x, y = 0.05, y-h/2.

            #p = Circle(xy=(x,y), radius=dy/2., facecolor=axcolor)#,transform=ax.transAxes)
            p = Rectangle(xy=(x,y), width=0.9, height=h, facecolor=axcolor,transform=ax.transAxes)


            l1 = Rectangle(xy=(x,y), width=0.9, height=h, facecolor='red',alpha=0.5, transform=ax.transAxes)
            #l2 = Line2D([x, x+w], [y, y+h], **lineparams)

            l1.set_visible(actives[cnt])
            #l2.set_visible(actives[cnt])
            self.labels.append(t)
            self.rectangles.append(p)
            self.lines.append((l1))#,l2))
            ax.add_patch(p)
            ax.add_line(l1)
            #ax.add_line(l2)
            cnt += 1

        ax.figure.canvas.mpl_connect('button_press_event', self._clicked)
        self.ax = ax


        self.cnt = 0
        self.observers = {}
    def _clicked(self, event):
        if event.button !=1 : return
        if event.inaxes != self.ax: return

        for p,t,lines in zip(self.rectangles, self.labels, self.lines):
            if (t.get_window_extent().contains(event.x, event.y) or
                p.get_window_extent().contains(event.x, event.y) ):
                l1=lines #, l2 = lines
                l1.set_visible(not l1.get_visible())
                #l2.set_visible(not l2.get_visible())
                thist = t
                break
        else:
            return


        if self.drawon: self.ax.figure.canvas.draw()

        if not self.eventson: return
        for cid, func in self.observers.items():
            func(thist.get_text())

    def on_clicked(self, func):
        """
        When the button is clicked, call *func* with button label

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        'remove the observer with connection id *cid*'
        try: del self.observers[cid]
        except KeyError: pass
