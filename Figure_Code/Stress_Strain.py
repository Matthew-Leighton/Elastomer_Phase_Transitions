import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib import rc
import matplotlib.gridspec as gridspec
from pyswarm import pso


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)