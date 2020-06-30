from PyQt5.QtCore import pyqtSlot

from ui.QtTest import Ui_Form
import scipy.signal as signal
from numpy import linspace, logspace, cos, sin, heaviside, log10, floor, zeros
from math import pi

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from PyQt5.QtWidgets import QWidget
from PyQt5 import QtGui, QtWidgets


class myWidget(QWidget, Ui_Form):
    def __init__(self):
        super(myWidget, self).__init__()
        self.setupUi(self)
        self.setWindowTitle("Electrotecnia - TP Final")

        self.filter_flag = 'none'
        self.first_pole_input.hide()
        self.char_value_input.hide()
        self.pole_label.hide()

        # Input Validators
        self.first_pole_real.setValidator(QtGui.QDoubleValidator(-10e12, 0, 8, self))
        self.first_pole_im.setValidator(QtGui.QDoubleValidator(-10e12, 10e12, 8, self))

        self.w0_input.setValidator(QtGui.QDoubleValidator(0, 10e12, 8, self))
        self.xi_input.setValidator(QtGui.QDoubleValidator(0, 10e12, 8, self))

        self.amplitude_input.setValidator(QtGui.QDoubleValidator(0, 10e12, 8, self))
        self.freq_input.setValidator(QtGui.QDoubleValidator(0, 10e12, 8, self))

        self.square_amplitude_input.setValidator(QtGui.QDoubleValidator(0, 10e12, 8, self))
        self.square_duty_input.setValidator(QtGui.QDoubleValidator(1e-12, 1, 8, self))

        self.ut_amplitude_input.setValidator(QtGui.QDoubleValidator(0, 10e12, 8, self))

        self.gain_input.setValidator(QtGui.QDoubleValidator(0, 10e12, 8, self))

        # Canvas
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.axes = self.figure.add_subplot()
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.canvas)

        canvas_index = self.plotter_container.addWidget(self.canvas)
        self.plotter_container.setCurrentIndex(canvas_index)

        # Filter buttons
        # First Order
        self.low_pass_1.clicked.connect(self.set_1st_lowpass)
        self.high_pass_1.clicked.connect(self.set_1st_highpass)
        self.high_all_pass_1.clicked.connect(self.set_1st_high_allpass)
        self.low_all_pass_1.clicked.connect(self.set_1st_low_allpass)
        # Second Order
        self.high_pass_2.clicked.connect(self.set_2nd_highpass)
        self.low_pass_2.clicked.connect(self.set_2nd_lowpass)
        self.band_pass.clicked.connect(self.set_bandpass)
        self.notch.clicked.connect(self.set_notch)
        self.high_all_pass_2.clicked.connect(self.set_2nd_high_allpass)
        self.low_all_pass_2.clicked.connect(self.set_2nd_low_allpass)

        # Graph buttons
        self.bode_plot_button.clicked.connect(lambda: self.plot('bode'))
        self.phase_plot_button.clicked.connect(lambda: self.plot('phase'))
        self.pole_plot_button.clicked.connect(lambda: self.plot('pole'))
        self.sine_plot_button.clicked.connect(lambda: self.plot('sine'))
        self.cosine_plot_button.clicked.connect(lambda: self.plot('cosine'))
        self.u_plot_button.clicked.connect(lambda: self.plot('ut'))
        self.square_plot_button.clicked.connect(lambda: self.plot('square'))
        self.sawtooth_plot_button.clicked.connect(lambda: self.plot('sawtooth'))
        self.delta_plot_button.clicked.connect(lambda: self.plot('delta'))

    @pyqtSlot()
    def set_1st_lowpass(self):
        self.char_value_input.hide()
        self.first_pole_input.show()
        self.first_pole_im.hide()
        self.first_pole_im_label.hide()
        self.filter_flag = '1st low pass'
        self.filter_label.setText('1st Order Low Pass Filter')

    def set_1st_highpass(self):
        self.char_value_input.hide()
        self.first_pole_input.show()
        self.first_pole_im.hide()
        self.first_pole_im_label.hide()
        self.filter_flag = '1st high pass'
        self.filter_label.setText('1st Order High Pass Filter')

    def set_1st_high_allpass(self):
        self.char_value_input.hide()
        self.first_pole_input.show()
        self.first_pole_im.hide()
        self.first_pole_im_label.hide()
        self.filter_flag = '1st high all pass'
        self.filter_label.setText('1st Order High All Pass Filter')

    def set_1st_low_allpass(self):
        self.char_value_input.hide()
        self.first_pole_input.show()
        self.first_pole_im.hide()
        self.first_pole_im_label.hide()
        self.filter_flag = '1st low all pass'
        self.filter_label.setText('1st Order Low All Pass Filter')

    def set_2nd_lowpass(self):
        self.char_value_input.show()
        self.first_pole_input.hide()
        self.filter_flag = '2nd low pass'
        self.filter_label.setText('2nd Order Low Pass Filter')

    def set_2nd_highpass(self):
        self.char_value_input.show()
        self.first_pole_input.hide()
        self.filter_flag = '2nd high pass'
        self.filter_label.setText('2nd Order High Pass Filter')

    def set_bandpass(self):
        self.char_value_input.show()
        self.first_pole_input.hide()
        self.filter_flag = 'band pass'
        self.filter_label.setText('2nd Order Band Pass Filter')

    def set_notch(self):
        self.char_value_input.show()
        self.first_pole_input.hide()
        self.filter_flag = 'notch'
        self.filter_label.setText('2nd Order Notch Filter')

    def set_2nd_high_allpass(self):
        self.char_value_input.show()
        self.first_pole_input.hide()
        self.filter_flag = '2nd high all pass'
        self.filter_label.setText('2nd Order High All Pass Filter')

    def set_2nd_low_allpass(self):
        self.char_value_input.show()
        self.first_pole_input.hide()
        self.filter_flag = '2nd low all pass'
        self.filter_label.setText('2nd Order Low All Pass Filter')

    def plot(self, flag):

        order_of_magnitude = 1

        if self.filter_flag == '1st low pass':
            pole = (float(self.first_pole_real.text()))
            P = [-1 * pole]
            Q = [1, -1 * pole]

            order_of_magnitude = floor(log10(-1 * pole))
        elif self.filter_flag == '1st high pass':
            pole = (float(self.first_pole_real.text()))
            P = [1, 0]
            Q = [1, -1 * pole]

            order_of_magnitude = floor(log10(-1 * pole))
        elif self.filter_flag == '1st high all pass':
            pole = (float(self.first_pole_real.text()))
            P = [1, pole]
            Q = [1, -1 * pole]

            order_of_magnitude = floor(log10(-1 * pole))
        elif self.filter_flag == '1st low all pass':
            pole = (float(self.first_pole_real.text()))
            P = [1, -1 * pole]
            Q = [1, pole]

            order_of_magnitude = floor(log10(-1 * pole))
        elif self.filter_flag == '2nd low pass':
            w0 = (float(self.w0_input.text()))
            xi = (float(self.xi_input.text()))
            P = [w0 ** 2]
            Q = [1, 2 * xi * w0, w0 ** 2]

            order_of_magnitude = floor(log10(w0))
        elif self.filter_flag == '2nd high pass':
            w0 = (float(self.w0_input.text()))
            xi = (float(self.xi_input.text()))
            P = [1, 0, 0]
            Q = [1, 2 * xi * w0, w0 ** 2]

            order_of_magnitude = floor(log10(w0))
        elif self.filter_flag == 'band pass':
            w0 = (float(self.w0_input.text()))
            xi = (float(self.xi_input.text()))
            P = [2 * xi * w0, 0]
            Q = [1, 2 * xi * w0, w0 ** 2]

            order_of_magnitude = floor(log10(w0))
        elif self.filter_flag == 'notch':
            w0 = (float(self.w0_input.text()))
            xi = (float(self.xi_input.text()))
            P = [1, 0, w0 ** 2]
            Q = [1, 2 * xi * w0, w0 ** 2]

            order_of_magnitude = floor(log10(w0))
        elif self.filter_flag == '2nd high all pass':
            w0 = (float(self.w0_input.text()))
            xi = (float(self.xi_input.text()))
            P = [1, 2 * xi * w0, w0 ** 2]
            Q = [1, -2 * xi * w0, w0 ** 2]

            order_of_magnitude = floor(log10(w0))
        elif self.filter_flag == '2nd low all pass':
            w0 = (float(self.w0_input.text()))
            xi = (float(self.xi_input.text()))
            P = [1, -2 * xi * w0, w0 ** 2]
            Q = [1, 2 * xi * w0, w0 ** 2]

            order_of_magnitude = floor(log10(w0))
        elif self.filter_flag == 'none':
            return

        left_limit = order_of_magnitude - 5
        right_limit = order_of_magnitude + 5

        gain = float(self.gain_input.text())
        i = 0
        while i < len(P):
            P[i] = gain * P[i]
            i += 1

        H = signal.TransferFunction(P, Q)

        x = logspace(left_limit, right_limit, num=1000)
        Bode = signal.bode(H, x)
        freq = Bode[0] / (2 * pi)

        self.axes.clear()

        if flag == 'bode':
            self.pole_label.hide()
            self.axes.set_title('Bode plot')
            self.axes.set_xlabel('f (log) [Hz]')
            self.axes.set_xscale('log')
            self.axes.set_ylabel('|H(jw)| [dB]')
            if self.filter_flag == '1st high all pass' or self.filter_flag == '1st low all pass' or self.filter_flag == '2nd high all pass' or self.filter_flag == '2nd low all pass':
                self.axes.plot(freq, zeros(len(freq)))
            else:
                self.axes.plot(freq, Bode[1])
        elif flag == 'phase':
            self.pole_label.hide()
            self.axes.set_title('Bode Phase plot')
            self.axes.set_xlabel('f (log) [Hz]')
            self.axes.set_xscale('log')
            self.axes.set_ylabel('/_H(jw) [deg]')
            self.axes.plot(freq, Bode[2])
        elif flag == 'pole':
            self.pole_label.setText('Triangle: Zeros\nCircle: Poles')
            self.pole_label.show()
            self.axes.set_title('Zero/Pole plot')
            self.axes.set_xlabel('Re(Z))')
            self.axes.set_ylabel('Im(Z)')
            transfer_poles = H.poles
            transfer_zeros = H.zeros
            self.axes.scatter(transfer_poles.real, transfer_poles.imag, marker='o', color='r')
            self.axes.scatter(transfer_zeros.real, transfer_zeros.imag, marker='^', color='b')
        elif flag == 'sine':
            self.pole_label.hide()
            A = float(self.amplitude_input.text())
            f = float(self.freq_input.text())
            x2 = linspace(0, 10 / f, num=1000)
            sinx = A * sin(2 * pi * f * x2)
            sine = signal.lsim(H, U=sinx, T=x2)

            self.axes.set_title('sin(t) response')
            self.axes.set_xlabel('t')
            self.axes.set_ylabel('Amplitude')

            self.axes.plot(sine[0], sine[1])
        elif flag == 'cosine':
            self.pole_label.hide()
            A = float(self.amplitude_input.text())
            f = float(self.freq_input.text())
            x2 = linspace(0, 10 / f, num=1000)
            cosx = A * cos(2 * pi * f * x2)
            cosine = signal.lsim(H, U=cosx, T=x2)

            self.axes.set_title('cos(t) response')
            self.axes.set_xlabel('t')
            self.axes.set_ylabel('Amplitude')

            self.axes.plot(cosine[0], cosine[1])
        elif flag == 'ut':
            self.pole_label.hide()
            A = float(self.ut_amplitude_input.text())
            x2 = linspace(0, 10, num=1000)
            ut = A * heaviside(x2, 0.5)
            step = signal.lsim(H, U=ut, T=x2)

            self.axes.set_title('u(t) response')
            self.axes.set_xlabel('t')
            self.axes.set_ylabel('Amplitude')

            self.axes.plot(step[0], step[1])
        elif flag == 'square':
            self.pole_label.hide()
            A = float(self.square_amplitude_input.text())
            dc = float(self.square_duty_input.text())
            x2 = linspace(0, 10 / dc, num=1000)
            sqx = A * signal.square(x2, dc)
            square = signal.lsim(H, U=sqx, T=x2)

            self.axes.set_title('Π(t) response')
            self.axes.set_xlabel('t')
            self.axes.set_ylabel('Amplitude')

            self.axes.plot(square[0], square[1])
        elif flag == 'sawtooth':
            self.pole_label.hide()
            A = float(self.square_amplitude_input.text())
            dc = float(self.square_duty_input.text())
            x2 = linspace(0, 10 / dc, num=1000)
            stx = A * signal.sawtooth(x2, dc)
            sawtooth = signal.lsim(H, U=stx, T=x2)

            self.axes.set_title('Λ(t) response')
            self.axes.set_xlabel('t')
            self.axes.set_ylabel('Amplitude')

            self.axes.plot(sawtooth[0], sawtooth[1])
        elif flag == 'delta':
            self.pole_label.hide()
            x2 = linspace(0, 10, num=1000)
            delta = signal.impulse(H, T=x2)

            self.axes.set_title('δ(t) response')
            self.axes.set_xlabel('t')
            self.axes.set_ylabel('Amplitude')

            self.axes.plot(delta[0], delta[1])

        self.axes.grid(True, which='both')
        self.axes.minorticks_on()

        self.canvas.draw()