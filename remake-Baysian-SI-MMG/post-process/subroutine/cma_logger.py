import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

class Logger:
    """Logger for dd-CMA"""
    def __init__(self, prefix='opt_result/', variable_list=['xmean', 'D', 'S', 'sigma', 'beta']):
        """
        Parameters
        ----------
        cma : DdCma instance
        prefix : string
            prefix for the log file path
        variable_list : list of string
            list of names of attributes of `cma` to be monitored
        """
        self.neval_offset = 0
        self.t_offset = 0
        self.prefix = prefix
        self.variable_list = variable_list
        self.logger = dict()
        # self.fmin_logger = self.prefix + 'fmin.dat'
        # with open(self.fmin_logger, 'w') as f:
        #     f.write('#' + type(self).__name__ + "\n")
        # for key in self.variable_list:
        #     self.logger[key] = self.prefix + key + '.dat'
        #     with open(self.logger[key], 'w') as f:
        #         f.write('#' + type(self).__name__ + "\n")
                
    def my_formatter(self, x, pos):
        """Float Number Format for Axes"""
        float_str = "{0:2.1e}".format(x)
        if self.latexmode:
            if "e" in float_str:
                base, exponent = float_str.split("e")
                return r"{0}e{1}".format(base, int(exponent))
            else:
                return r"" + float_str + ""
        else:
            if "e" in float_str:
                base, exponent = float_str.split("e")
                return "{0}e{1}".format(base, int(exponent))
            else:
                return "" + float_str + ""
    def plot(self,
             xaxis=0,
             ncols=None,
             figsize=None,
             cmap_='Spectral',
             latexmode=False):
        
        """Plot the result
        Parameters
        ----------
        xaxis : int, optional (default = 0)
            0. vs iterations
            1. vs function evaluations
        ncols : int, optional (default = None)
            number of columns
        figsize : tuple, optional (default = None)
            figure size
        cmap_ : string, optional (default = 'spectral')
            cmap
        latexmode : bool, optional (default = False)
            LaTeX is enabled if this is True
        
        Returns
        -------
        fig : figure object.
            figure object
        axdict : dictionary of axes
            the keys are the names of variables given in `variable_list`
        """
        self.latexmode = latexmode 

        mpl.rc('lines', linewidth=2, markersize=8)
        mpl.rc('font', size=12)
        mpl.rc('grid', color='0.75', linestyle=':')
        # if self.latexmode:
        #     mpl.rc('text', usetex=True)  # for a paper submision
        #     mpl.rc('ps', useafm=True)  # Force to use
        #     mpl.rc('pdf', use14corefonts=True)  # only Type 1 fonts
        prefix = self.prefix
        variable_list = self.variable_list

        # Default settings
        nfigs = 1 + len(variable_list)
        if ncols is None:
            ncols = int(np.ceil(np.sqrt(nfigs)))
        nrows = int(np.ceil(nfigs / ncols))
        if figsize is None:
            figsize = (4 * ncols, 3 * nrows)
        axdict = dict()
        
        # Figure
        fig = plt.figure(figsize=figsize)
        # The first figure
        x = np.loadtxt('/Users/asakanashi/Documents/SI_MMG/1113_Baysian_MMG/post-process/opt_result/fmin.dat')
        x = x[~np.isnan(x[:, xaxis]), :]  # remove columns where xaxis is nan
        # Axis
        ax = plt.subplot(nrows, ncols, 1)
        ax.set_title('fmin')
        ax.grid(True)
        ax.grid(which='major', linewidth=0.50)
        ax.grid(which='minor', linewidth=0.25)
        plt.plot(x[:, xaxis], x[:, 2:])
        ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(self.my_formatter))
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(self.my_formatter))
        axdict['fmin'] = ax

        # The other figures
        idx = 1
        for key in variable_list:
            idx += 1
            x = np.loadtxt('log/' +  key + '.dat')
            # x = np.genfromtxt(prefix + '_' + key + '.dat',dtype=None)

            x = x[~np.isnan(x[:, xaxis]), :]  # remove columns where xaxis is nan
            ax = plt.subplot(nrows, ncols, idx)
            if self.latexmode:
                ax.set_title(r'\detokenize{' + key + '}')
            else:
                ax.set_title(key)
            ax.grid(True)
            ax.grid(which='major', linewidth=0.50)
            ax.grid(which='minor', linewidth=0.25)
            cmap = plt.get_cmap(cmap_)
            cNorm = mpl.colors.Normalize(vmin=0, vmax=x.shape[1] - 2)
            scalarMap = mpl.cm.ScalarMappable(norm=cNorm, cmap=cmap)
            for i in range(x.shape[1] - 2):
                plt.plot(
                    x[:, xaxis], x[:, 2 + i], color=scalarMap.to_rgba(i))
            ax.xaxis.set_major_formatter(
                mpl.ticker.FuncFormatter(self.my_formatter))
            ax.yaxis.set_major_formatter(
                mpl.ticker.FuncFormatter(self.my_formatter))
            axdict[key] = ax

        plt.tight_layout() # NOTE: not sure if it works fine
        return fig, axdict