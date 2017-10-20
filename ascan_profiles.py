import numpy as np
import matplotlib.pyplot as plt
from App import readCirrusOct as cirrusOct


def smooth(x,window_len=11,window='hanning', align=True):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    if align:
        y = y[window_len:]
        y = np.append(y, np.zeros(window_len))
    return y

if __name__ == '__main__':
    fname = 'Data/Sample/Cirrus/Sample1/Macular Cube 512x128_01-01-2001_01-01-01_OS_sn41343_cube_z.img'
    x_coord = 64
    y_coord = 256

    window_size = np.arange(2, 25, 3)
    ncols = 3

    nrows = int(np.ceil(len(window_size) / ncols))

    oct = cirrusOct.readCirrusOct(fname)
    line = cirrusOct.getAScan(oct, x_coord, y_coord)

    xvals = range(len(line))

    window = 'blackman'


    fig, axarr = plt.subplots(nrows, ncols, sharex=True, sharey=True)

    lines =[smooth(line, ws, window) for ws in window_size]

    for idx,smoothed_line  in enumerate(lines):
        row = int(np.floor(idx / ncols))
        col = int(idx % ncols)
        axarr[row, col].set_title('Size:{}'.format(window_size[idx]))

        axarr[row, col].imshow(oct['image_data'][x_coord,: , :])
        axarr[row, col].plot(smoothed_line + y_coord, range(len(smoothed_line)))
        axarr[row, col].set_ylim([400, 850])
        axarr[row, col].set_xlim([0, 512])

    plt.show()
