import sys, os
import numpy as np
#import scipy
#import skimage

from scipy import ndimage as ndi

# from skimage.transforms; standalone version of the skimage resize methods
# REF https://github.com/scikit-image/scikit-image/blob/v0.22.0/skimage/transform/_warps.py
def _preprocess_resize_output_shape(image, output_shape):
    output_shape = tuple(output_shape)
    output_ndim = len(output_shape)
    input_shape = image.shape
    if output_ndim > image.ndim:
        # append dimensions to input_shape
        input_shape += (1, ) * (output_ndim - image.ndim)
        image = np.reshape(image, input_shape)
    elif output_ndim == image.ndim - 1:
        # multichannel case: append shape of last axis
        output_shape = output_shape + (image.shape[-1], )
    elif output_ndim < image.ndim:
        raise ValueError("output_shape length cannot be smaller than the "
                         "image number of dimensions")

    return image, output_shape

def resize(image, output_shape, 
           order=None, mode='reflect', 
           cval=0, 
           clip=True,
           preserve_range=False, 
           anti_aliasing=None, 
           anti_aliasing_sigma=None):

    image, output_shape = _preprocess_resize_output_shape(image, output_shape)
    input_shape = image.shape
    input_type = image.dtype

    if input_type == np.float16:
        image = image.astype(np.float32)

    if anti_aliasing is None:
        anti_aliasing = (
            not input_type == bool and
            not (np.issubdtype(input_type, np.integer) and order == 0) and
            any(x < y for x, y in zip(output_shape, input_shape)))

    if input_type == bool and anti_aliasing:
        raise ValueError("anti_aliasing must be False for boolean images")

    factors = np.divide(input_shape, output_shape)
    #order = _validate_interpolation_order(input_type, order)

    #if order > 0:
    #    image = convert_to_float(image, preserve_range)

    # Translate modes used by np.pad to those used by scipy.ndimage
    #ndi_mode = _to_ndimage_mode(mode)
    ndi_mode = mode
    if anti_aliasing:
        if anti_aliasing_sigma is None:
            anti_aliasing_sigma = np.maximum(0, (factors - 1) / 2)
        else:
            anti_aliasing_sigma = \
                np.atleast_1d(anti_aliasing_sigma) * np.ones_like(factors)
            if np.any(anti_aliasing_sigma < 0):
                raise ValueError("Anti-aliasing standard deviation must be "
                                 "greater than or equal to zero")
            elif np.any((anti_aliasing_sigma > 0) & (factors <= 1)):
                warn("Anti-aliasing standard deviation greater than zero but "
                     "not down-sampling along all axes")
        filtered = ndi.gaussian_filter(image, anti_aliasing_sigma,
                                       cval=cval, mode=ndi_mode)
    else:
        filtered = image

    zoom_factors = [1 / f for f in factors]
    out = ndi.zoom(filtered, zoom_factors, order=order, mode=ndi_mode,
                   cval=cval, grid_mode=True)

    #_clip_warp_output(image, out, mode, cval, clip)

    return out




# Render images to terminal
class TerminalPlot:

    def __init__(self, nx, ny):
        self.nx = nx
        self.ny = ny
        self.screen = np.zeros((self.nx, self.ny))

    def norm(self, x):
        if x < 0.05:
            return "  "
        elif x < 0.25:
            return ". "
        elif x < 0.4:
            return "o "
        elif x < 0.7:
            return "x "
        elif x < 0.9:
            return "X "
        else:
            return "W "

    def rescale(self, im):
        im2 = resize(im, (self.nx, self.ny), order=1)
        return im2

    
    #--------------------------------------------------
    # plot single panel
    def plot(self, 
             data, # array to plot
             name='',
             ):
        lines = self.gen_panel(data, name=name)

        for line in lines: print(line)


    #--------------------------------------------------
    def plot_panels(self, 
                    dims, # tuple of (rows, cols)
                   *all_data, # dictionary of panels):
                   ): 

        #--------------------------------------------------
        # prepare fig
        lines = []
        nrows, ncols = dims

        for i in range(nrows):
            for j in range(self.ny+2):
                lines.append('')
                #lines.append(str(j))

        #--------------------------------------------------
        # collect panels

        for ir in range(nrows):
            for jr in range(ncols):

                for data in all_data:
                    axs = data['axs']
                    if axs[0] == ir and axs[1] == jr:
                        l1 = self.gen_panel(data['data'], name=data['name']) # get panel

                        for i, line in enumerate(l1):
                            irr = ir*(self.ny+1) + i
                            
                            if ir < nrows-1 and i == self.ny+1: continue # skip last horizontal frame

                            if jr == 0: # first col
                                lines[irr] += line
                            else:
                                lines[irr] += line[1:] # skip first vertical frame


        #--------------------------------------------------
        for line in lines: print(line)


    #--------------------------------------------------
    # generate panel to be plotted
    def gen_panel(self, 
             data, # array to plot
             name='',
             ):

        # transpose and flip to orient the img properly
        data = np.fliplr(data.T) 

        #for i in range(self.ny):
        #    self.screen[i,i] = 1

        self.screen = self.rescale(data)

        lines = [] # screen is collected here

        #--------------------------------------------------
        # first line
        if name == '':
            line = '+'
            for i in range(self.nx):
                line += "--"
            line += '+'

            lines.append(line)
        else: # add name
            line = '+'
            line += '{:-{align}{width}}'.format(name, align='^', width=str(2*self.nx))
            line += '+'

            lines.append(line)

        #--------------------------------------------------
        # print content

        for j in range(self.ny):
            line = ""
            stripe = self.screen[:,j]

            line += "|"

            for i in range(self.nx):
                line += self.norm(stripe[i])

            line += "|"

            lines.append(line)

        #--------------------------------------------------
        # last line
        line = '+'
        for i in range(self.nx):
            line += "--"
        line += '+'
        lines.append(line)
        #--------------------------------------------------

        return lines
            

                
def print_format_table():
    """
    prints table of formatted text format options
    """

    for style in range(8):
        for fg in range(30,38):
            s1 = ''
            for bg in range(40,48):
                format = ';'.join([str(style), str(fg), str(bg)])
                s1 += '\x1b[%sm %s \x1b[0m' % (format, format)
            print(s1)
        print('\n')



if __name__ == "__main__":

    plt = TerminalPlot(12, 12)

    data = np.ones((64, 64))
    x = np.linspace(-3, 3, 64)
    y = np.linspace(-3, 3, 64)
    X,Y = np.meshgrid(x,y, indexing='ij')
    r = X**2 + Y**2
    #print(r)

    data[:,:] = np.exp(-r**2)
    #print(data)

    #plt.plot(data, name='ex')

    #print_format_table()

    plt.plot_panels(
            (3,3),
            dict(axs=(0,0), data=data, name='ex'),
            dict(axs=(0,1), data=data, name='ey'),
            dict(axs=(0,2), data=data, name='ez'),
            dict(axs=(1,0), data=data, name='bx'),
            dict(axs=(1,1), data=data, name='by'),
            dict(axs=(1,2), data=data, name='bz'),
            dict(axs=(2,0), data=data, name='e-'),
            dict(axs=(2,1), data=data, name='e+'),
            dict(axs=(2,2), data=data, name='ph'),
            )









