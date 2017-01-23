import matplotlib.pyplot as plt
import matplotlib.cm as cm
import platform
import pickle as pl

class plotting(object):
	"""Plotting class"""
	def __init__(self, num_bins=500, fig_size=(16,9), output='.', **kwargs):
		self.num_bins = num_bins
		self.fig_size = fig_size 
		self.output = output

		self.histtype = kwargs.get('histtype', u'stepfilled')
		self.facecolor = kwargs.get('facecolor', 'g')
		self.alpha = kwargs.get('alpha', 0.45)

		if platform.system() == 'Linux':
			self.my_cmap = cm.get_cmap('Spectral')
		else:
			self.my_cmap = cm.get_cmap('viridis')

		self.my_cmap.set_over('w')
		plt.ioff()
		plt.rc('text', usetex=True)
		
	def hist2D(self, x_val, y_val, **kwargs):
		self.output = kwargs.get('output', self.output)
		self.num_bins = kwargs.get('num_bins', self.num_bins)
		self.fig_size = kwargs.get('fig_size', self.fig_size)
		self.x_lab = kwargs.get('x_lab', r'x (units)')
		self.y_lab = kwargs.get('y_lab', r'y (units)')
		self.plot_name = kwargs.get('plot_name', 'hist2d.pdf')

		fig = plt.figure(num=None, figsize=self.fig_size, dpi=200, facecolor='w', edgecolor='k')
		ax = fig.gca()
		if 'range' in kwargs:
			plt.hist2d(x_val,y_val,bins=self.num_bins,range=kwargs['range'],cmap=self.my_cmap)
		else:
			plt.hist2d(x_val,y_val,bins=self.num_bins, cmap=self.my_cmap)

		plt.ylabel(self.y_lab, fontsize=18)
		plt.xlabel(self.x_lab, fontsize=20)
		plt.colorbar()

		if self.output[-1] != '/':
			self.output = self.output + '/'

		pl.dump(ax, file(self.output + self.plot_name + '.pkl','wb'))
		plt.savefig(self.output + self.plot_name)


	def hist1D(self, x_val, **kwargs):
		self.output = kwargs.get('output', self.output)
		self.num_bins = kwargs.get('num_bins', self.num_bins)
		self.fig_size = kwargs.get('fig_size', self.fig_size)
		self.x_lab = kwargs.get('x_lab', r'$x (units)$')
		self.y_lab = kwargs.get('y_lab', r'$counts (\#)$')
		self.plot_name = kwargs.get('plot_name', 'hist1d.pdf')

		fig = plt.figure(num=None, figsize=self.fig_size, dpi=200, facecolor='w', edgecolor='k')
		if 'range' in kwargs:
			plt.hist(x_val,	bins=self.num_bins,
				range=kwargs['range'],
				histtype=self.histtype,
				facecolor = self.facecolor,
				alpha = self.alpha)
		else:
			plt.hist(x_val,	bins=self.num_bins,
				histtype=self.histtype,
				facecolor = self.facecolor,
				alpha = self.alpha)

		plt.ylabel(self.y_lab, fontsize=18)
		plt.xlabel(self.x_lab, fontsize=20)

		if self.output[-1] != '/':
			self.output = self.output + '/'

		pl.dump(fig, file(self.output + self.plot_name + '.pkl','wb'))
		plt.savefig(self.output + self.plot_name)


	def WvsQ2(self,W,Q2,**kwargs):
		self.hist2D(W, Q2, 
			range=[[0,4],[0,10]], 
			x_lab = r'$W (GeV)$', 
			y_lab = r'$Q^{2}$ $(GeV^{2})$', 
			plot_name = 'WvsQ2.pdf',
			**kwargs)

		self.hist1D(W,
			range=[0,4],
			x_lab = r'$W (GeV)$',
			plot_name = 'W_Hist.pdf',
			**kwargs)

		self.hist1D(Q2,
			range=[0,10],
			x_lab = r'$Q^{2}$ $(GeV^{2})$',
			plot_name = 'Q2_Hist.pdf',
			**kwargs)




