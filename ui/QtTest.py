# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'QtTest.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(1611, 980)
        self.plotter_container = QtWidgets.QStackedWidget(Form)
        self.plotter_container.setGeometry(QtCore.QRect(0, 0, 821, 701))
        self.plotter_container.setObjectName("plotter_container")
        self.page = QtWidgets.QWidget()
        self.page.setObjectName("page")
        self.plotter_container.addWidget(self.page)
        self.page_2 = QtWidgets.QWidget()
        self.page_2.setObjectName("page_2")
        self.plotter_container.addWidget(self.page_2)
        self.filter_label = QtWidgets.QLabel(Form)
        self.filter_label.setGeometry(QtCore.QRect(890, 870, 391, 81))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.filter_label.setFont(font)
        self.filter_label.setObjectName("filter_label")
        self.tabWidget = QtWidgets.QTabWidget(Form)
        self.tabWidget.setGeometry(QtCore.QRect(10, 760, 851, 211))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.tabWidget.setFont(font)
        self.tabWidget.setObjectName("tabWidget")
        self.tab_1 = QtWidgets.QWidget()
        self.tab_1.setObjectName("tab_1")
        self.low_pass_1 = QtWidgets.QPushButton(self.tab_1)
        self.low_pass_1.setGeometry(QtCore.QRect(20, 40, 151, 81))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.low_pass_1.setFont(font)
        self.low_pass_1.setObjectName("low_pass_1")
        self.high_pass_1 = QtWidgets.QPushButton(self.tab_1)
        self.high_pass_1.setGeometry(QtCore.QRect(230, 40, 151, 81))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.high_pass_1.setFont(font)
        self.high_pass_1.setObjectName("high_pass_1")
        self.high_all_pass_1 = QtWidgets.QPushButton(self.tab_1)
        self.high_all_pass_1.setGeometry(QtCore.QRect(430, 40, 151, 81))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.high_all_pass_1.setFont(font)
        self.high_all_pass_1.setObjectName("high_all_pass_1")
        self.low_all_pass_1 = QtWidgets.QPushButton(self.tab_1)
        self.low_all_pass_1.setGeometry(QtCore.QRect(630, 40, 151, 81))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.low_all_pass_1.setFont(font)
        self.low_all_pass_1.setObjectName("low_all_pass_1")
        self.tabWidget.addTab(self.tab_1, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.low_pass_2 = QtWidgets.QPushButton(self.tab_2)
        self.low_pass_2.setGeometry(QtCore.QRect(30, 10, 171, 61))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.low_pass_2.setFont(font)
        self.low_pass_2.setObjectName("low_pass_2")
        self.high_pass_2 = QtWidgets.QPushButton(self.tab_2)
        self.high_pass_2.setGeometry(QtCore.QRect(230, 10, 171, 61))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.high_pass_2.setFont(font)
        self.high_pass_2.setObjectName("high_pass_2")
        self.high_all_pass_2 = QtWidgets.QPushButton(self.tab_2)
        self.high_all_pass_2.setGeometry(QtCore.QRect(430, 10, 171, 61))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.high_all_pass_2.setFont(font)
        self.high_all_pass_2.setObjectName("high_all_pass_2")
        self.band_pass = QtWidgets.QPushButton(self.tab_2)
        self.band_pass.setGeometry(QtCore.QRect(30, 90, 171, 61))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.band_pass.setFont(font)
        self.band_pass.setObjectName("band_pass")
        self.notch = QtWidgets.QPushButton(self.tab_2)
        self.notch.setGeometry(QtCore.QRect(230, 90, 171, 61))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.notch.setFont(font)
        self.notch.setObjectName("notch")
        self.low_all_pass_2 = QtWidgets.QPushButton(self.tab_2)
        self.low_all_pass_2.setGeometry(QtCore.QRect(430, 90, 171, 61))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.low_all_pass_2.setFont(font)
        self.low_all_pass_2.setObjectName("low_all_pass_2")
        self.low_pass_notch = QtWidgets.QPushButton(self.tab_2)
        self.low_pass_notch.setGeometry(QtCore.QRect(630, 10, 171, 61))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.low_pass_notch.setFont(font)
        self.low_pass_notch.setObjectName("low_pass_notch")
        self.high_pass_notch = QtWidgets.QPushButton(self.tab_2)
        self.high_pass_notch.setGeometry(QtCore.QRect(630, 90, 171, 61))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.high_pass_notch.setFont(font)
        self.high_pass_notch.setObjectName("high_pass_notch")
        self.tabWidget.addTab(self.tab_2, "")
        self.label_2 = QtWidgets.QLabel(Form)
        self.label_2.setGeometry(QtCore.QRect(20, 730, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.first_pole_input = QtWidgets.QGroupBox(Form)
        self.first_pole_input.setGeometry(QtCore.QRect(840, 300, 361, 131))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.first_pole_input.setFont(font)
        self.first_pole_input.setObjectName("first_pole_input")
        self.first_pole_re_label = QtWidgets.QLabel(self.first_pole_input)
        self.first_pole_re_label.setGeometry(QtCore.QRect(20, 40, 101, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.first_pole_re_label.setFont(font)
        self.first_pole_re_label.setObjectName("first_pole_re_label")
        self.first_pole_im = QtWidgets.QLineEdit(self.first_pole_input)
        self.first_pole_im.setGeometry(QtCore.QRect(190, 75, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(20)
        self.first_pole_im.setFont(font)
        self.first_pole_im.setObjectName("first_pole_im")
        self.first_pole_real = QtWidgets.QLineEdit(self.first_pole_input)
        self.first_pole_real.setGeometry(QtCore.QRect(20, 75, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(20)
        self.first_pole_real.setFont(font)
        self.first_pole_real.setObjectName("first_pole_real")
        self.first_pole_im_label = QtWidgets.QLabel(self.first_pole_input)
        self.first_pole_im_label.setGeometry(QtCore.QRect(190, 40, 101, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.first_pole_im_label.setFont(font)
        self.first_pole_im_label.setObjectName("first_pole_im_label")
        self.pole_label = QtWidgets.QLabel(Form)
        self.pole_label.setGeometry(QtCore.QRect(840, 590, 211, 111))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.pole_label.setFont(font)
        self.pole_label.setObjectName("pole_label")
        self.first_pole_im_label_4 = QtWidgets.QLabel(Form)
        self.first_pole_im_label_4.setGeometry(QtCore.QRect(880, 790, 101, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.first_pole_im_label_4.setFont(font)
        self.first_pole_im_label_4.setObjectName("first_pole_im_label_4")
        self.gain_input = QtWidgets.QLineEdit(Form)
        self.gain_input.setGeometry(QtCore.QRect(880, 825, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(20)
        self.gain_input.setFont(font)
        self.gain_input.setObjectName("gain_input")
        self.char_value_input = QtWidgets.QGroupBox(Form)
        self.char_value_input.setGeometry(QtCore.QRect(840, 440, 361, 131))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.char_value_input.setFont(font)
        self.char_value_input.setObjectName("char_value_input")
        self.first_pole_re_label_2 = QtWidgets.QLabel(self.char_value_input)
        self.first_pole_re_label_2.setGeometry(QtCore.QRect(20, 40, 101, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.first_pole_re_label_2.setFont(font)
        self.first_pole_re_label_2.setObjectName("first_pole_re_label_2")
        self.xi_input = QtWidgets.QLineEdit(self.char_value_input)
        self.xi_input.setGeometry(QtCore.QRect(190, 75, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(20)
        self.xi_input.setFont(font)
        self.xi_input.setObjectName("xi_input")
        self.w0_input = QtWidgets.QLineEdit(self.char_value_input)
        self.w0_input.setGeometry(QtCore.QRect(20, 75, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(20)
        self.w0_input.setFont(font)
        self.w0_input.setObjectName("w0_input")
        self.first_pole_im_label_2 = QtWidgets.QLabel(self.char_value_input)
        self.first_pole_im_label_2.setGeometry(QtCore.QRect(190, 40, 101, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.first_pole_im_label_2.setFont(font)
        self.first_pole_im_label_2.setObjectName("first_pole_im_label_2")
        self.groupBox = QtWidgets.QGroupBox(Form)
        self.groupBox.setGeometry(QtCore.QRect(840, 20, 361, 271))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.groupBox.setFont(font)
        self.groupBox.setObjectName("groupBox")
        self.bode_plot_button = QtWidgets.QPushButton(self.groupBox)
        self.bode_plot_button.setGeometry(QtCore.QRect(10, 30, 331, 71))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.bode_plot_button.setFont(font)
        self.bode_plot_button.setObjectName("bode_plot_button")
        self.phase_plot_button = QtWidgets.QPushButton(self.groupBox)
        self.phase_plot_button.setGeometry(QtCore.QRect(10, 110, 331, 71))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.phase_plot_button.setFont(font)
        self.phase_plot_button.setObjectName("phase_plot_button")
        self.pole_plot_button = QtWidgets.QPushButton(self.groupBox)
        self.pole_plot_button.setGeometry(QtCore.QRect(10, 190, 331, 71))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.pole_plot_button.setFont(font)
        self.pole_plot_button.setObjectName("pole_plot_button")
        self.groupBox_2 = QtWidgets.QGroupBox(Form)
        self.groupBox_2.setGeometry(QtCore.QRect(1210, 20, 391, 271))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.groupBox_2.setFont(font)
        self.groupBox_2.setObjectName("groupBox_2")
        self.sine_plot_button = QtWidgets.QPushButton(self.groupBox_2)
        self.sine_plot_button.setGeometry(QtCore.QRect(10, 50, 161, 51))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.sine_plot_button.setFont(font)
        self.sine_plot_button.setObjectName("sine_plot_button")
        self.cosine_plot_button = QtWidgets.QPushButton(self.groupBox_2)
        self.cosine_plot_button.setGeometry(QtCore.QRect(210, 50, 161, 51))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.cosine_plot_button.setFont(font)
        self.cosine_plot_button.setObjectName("cosine_plot_button")
        self.signal_input = QtWidgets.QGroupBox(self.groupBox_2)
        self.signal_input.setGeometry(QtCore.QRect(20, 110, 361, 131))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.signal_input.setFont(font)
        self.signal_input.setObjectName("signal_input")
        self.first_pole_re_label_3 = QtWidgets.QLabel(self.signal_input)
        self.first_pole_re_label_3.setGeometry(QtCore.QRect(20, 40, 101, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.first_pole_re_label_3.setFont(font)
        self.first_pole_re_label_3.setObjectName("first_pole_re_label_3")
        self.freq_input = QtWidgets.QLineEdit(self.signal_input)
        self.freq_input.setGeometry(QtCore.QRect(190, 75, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(20)
        self.freq_input.setFont(font)
        self.freq_input.setObjectName("freq_input")
        self.amplitude_input = QtWidgets.QLineEdit(self.signal_input)
        self.amplitude_input.setGeometry(QtCore.QRect(20, 75, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(20)
        self.amplitude_input.setFont(font)
        self.amplitude_input.setObjectName("amplitude_input")
        self.first_pole_im_label_3 = QtWidgets.QLabel(self.signal_input)
        self.first_pole_im_label_3.setGeometry(QtCore.QRect(190, 40, 121, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.first_pole_im_label_3.setFont(font)
        self.first_pole_im_label_3.setObjectName("first_pole_im_label_3")
        self.groupBox_3 = QtWidgets.QGroupBox(Form)
        self.groupBox_3.setGeometry(QtCore.QRect(1210, 310, 391, 121))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.groupBox_3.setFont(font)
        self.groupBox_3.setObjectName("groupBox_3")
        self.signal_input_2 = QtWidgets.QGroupBox(self.groupBox_3)
        self.signal_input_2.setGeometry(QtCore.QRect(190, 10, 191, 101))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.signal_input_2.setFont(font)
        self.signal_input_2.setObjectName("signal_input_2")
        self.first_pole_re_label_4 = QtWidgets.QLabel(self.signal_input_2)
        self.first_pole_re_label_4.setGeometry(QtCore.QRect(10, 20, 101, 31))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.first_pole_re_label_4.setFont(font)
        self.first_pole_re_label_4.setObjectName("first_pole_re_label_4")
        self.ut_amplitude_input = QtWidgets.QLineEdit(self.signal_input_2)
        self.ut_amplitude_input.setGeometry(QtCore.QRect(10, 50, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(20)
        self.ut_amplitude_input.setFont(font)
        self.ut_amplitude_input.setObjectName("ut_amplitude_input")
        self.u_plot_button = QtWidgets.QPushButton(self.groupBox_3)
        self.u_plot_button.setGeometry(QtCore.QRect(20, 50, 151, 51))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.u_plot_button.setFont(font)
        self.u_plot_button.setObjectName("u_plot_button")
        self.groupBox_5 = QtWidgets.QGroupBox(Form)
        self.groupBox_5.setGeometry(QtCore.QRect(1210, 450, 391, 271))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.groupBox_5.setFont(font)
        self.groupBox_5.setObjectName("groupBox_5")
        self.square_plot_button = QtWidgets.QPushButton(self.groupBox_5)
        self.square_plot_button.setGeometry(QtCore.QRect(10, 50, 161, 51))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.square_plot_button.setFont(font)
        self.square_plot_button.setObjectName("square_plot_button")
        self.signal_input_4 = QtWidgets.QGroupBox(self.groupBox_5)
        self.signal_input_4.setGeometry(QtCore.QRect(20, 110, 361, 131))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.signal_input_4.setFont(font)
        self.signal_input_4.setObjectName("signal_input_4")
        self.first_pole_re_label_6 = QtWidgets.QLabel(self.signal_input_4)
        self.first_pole_re_label_6.setGeometry(QtCore.QRect(20, 40, 101, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.first_pole_re_label_6.setFont(font)
        self.first_pole_re_label_6.setObjectName("first_pole_re_label_6")
        self.square_duty_input = QtWidgets.QLineEdit(self.signal_input_4)
        self.square_duty_input.setGeometry(QtCore.QRect(190, 75, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(20)
        self.square_duty_input.setFont(font)
        self.square_duty_input.setObjectName("square_duty_input")
        self.square_amplitude_input = QtWidgets.QLineEdit(self.signal_input_4)
        self.square_amplitude_input.setGeometry(QtCore.QRect(20, 75, 151, 41))
        font = QtGui.QFont()
        font.setPointSize(20)
        self.square_amplitude_input.setFont(font)
        self.square_amplitude_input.setObjectName("square_amplitude_input")
        self.first_pole_im_label_5 = QtWidgets.QLabel(self.signal_input_4)
        self.first_pole_im_label_5.setGeometry(QtCore.QRect(190, 40, 161, 31))
        font = QtGui.QFont()
        font.setPointSize(9)
        self.first_pole_im_label_5.setFont(font)
        self.first_pole_im_label_5.setObjectName("first_pole_im_label_5")
        self.sawtooth_plot_button = QtWidgets.QPushButton(self.groupBox_5)
        self.sawtooth_plot_button.setGeometry(QtCore.QRect(200, 50, 161, 51))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.sawtooth_plot_button.setFont(font)
        self.sawtooth_plot_button.setObjectName("sawtooth_plot_button")
        self.groupBox_4 = QtWidgets.QGroupBox(Form)
        self.groupBox_4.setGeometry(QtCore.QRect(1220, 730, 391, 121))
        font = QtGui.QFont()
        font.setPointSize(12)
        self.groupBox_4.setFont(font)
        self.groupBox_4.setObjectName("groupBox_4")
        self.delta_plot_button = QtWidgets.QPushButton(self.groupBox_4)
        self.delta_plot_button.setGeometry(QtCore.QRect(20, 50, 151, 51))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.delta_plot_button.setFont(font)
        self.delta_plot_button.setObjectName("delta_plot_button")

        self.retranslateUi(Form)
        self.tabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.filter_label.setText(_translate("Form", "Type of Filter"))
        self.low_pass_1.setText(_translate("Form", "\n"
"Low Pass Filter\n"
""))
        self.high_pass_1.setText(_translate("Form", "\n"
"High Pass Filter\n"
""))
        self.high_all_pass_1.setText(_translate("Form", "\n"
"High All Pass Filter\n"
""))
        self.low_all_pass_1.setText(_translate("Form", "\n"
"Low All Pass Filter\n"
""))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_1), _translate("Form", "First Order"))
        self.low_pass_2.setText(_translate("Form", "\n"
"Low Pass Filter\n"
""))
        self.high_pass_2.setText(_translate("Form", "\n"
"High Pass Filter\n"
""))
        self.high_all_pass_2.setText(_translate("Form", "\n"
"High All Pass Filter\n"
""))
        self.band_pass.setText(_translate("Form", "\n"
"Band Pass Filter\n"
""))
        self.notch.setText(_translate("Form", "\n"
"Notch Filter\n"
""))
        self.low_all_pass_2.setText(_translate("Form", "\n"
"Low All Pass Filter\n"
""))
        self.low_pass_notch.setText(_translate("Form", "\n"
"Low-Pass Notch Filter\n"
""))
        self.high_pass_notch.setText(_translate("Form", "\n"
"High-Pass Notch Filter\n"
""))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("Form", "Second Order"))
        self.label_2.setText(_translate("Form", "Filters:"))
        self.first_pole_input.setTitle(_translate("Form", "First Pole"))
        self.first_pole_re_label.setText(_translate("Form", "Re(pole):"))
        self.first_pole_im.setText(_translate("Form", "0"))
        self.first_pole_real.setText(_translate("Form", "-1"))
        self.first_pole_im_label.setText(_translate("Form", "Im(pole):"))
        self.pole_label.setText(_translate("Form", "Triangles: Zeros\n"
"Circles: Poles"))
        self.first_pole_im_label_4.setText(_translate("Form", " Filter Gain:"))
        self.gain_input.setText(_translate("Form", "1"))
        self.char_value_input.setTitle(_translate("Form", "Characteristic Values"))
        self.first_pole_re_label_2.setText(_translate("Form", "ω0 [rad/s]:"))
        self.xi_input.setText(_translate("Form", "0.5"))
        self.w0_input.setText(_translate("Form", "1"))
        self.first_pole_im_label_2.setText(_translate("Form", " ξ:"))
        self.groupBox.setTitle(_translate("Form", "Bode Graphs"))
        self.bode_plot_button.setText(_translate("Form", "Bode Plot"))
        self.phase_plot_button.setText(_translate("Form", "Phase Plot"))
        self.pole_plot_button.setText(_translate("Form", "Zero/Pole Plot"))
        self.groupBox_2.setTitle(_translate("Form", "Sinusoid Signal Input"))
        self.sine_plot_button.setText(_translate("Form", "sin(t) Response"))
        self.cosine_plot_button.setText(_translate("Form", "cos(t) Response"))
        self.signal_input.setTitle(_translate("Form", "Signal Properties"))
        self.first_pole_re_label_3.setText(_translate("Form", "Amplitude:"))
        self.freq_input.setText(_translate("Form", "1"))
        self.amplitude_input.setText(_translate("Form", "1"))
        self.first_pole_im_label_3.setText(_translate("Form", "Frequency [Hz]:"))
        self.groupBox_3.setTitle(_translate("Form", "Step Signal Input"))
        self.signal_input_2.setTitle(_translate("Form", "Signal Properties"))
        self.first_pole_re_label_4.setText(_translate("Form", "Amplitude:"))
        self.ut_amplitude_input.setText(_translate("Form", "1"))
        self.u_plot_button.setText(_translate("Form", "u(t) Response"))
        self.groupBox_5.setTitle(_translate("Form", "Square/Sawtooth Signal Input"))
        self.square_plot_button.setText(_translate("Form", "Π(t) Response"))
        self.signal_input_4.setTitle(_translate("Form", "Signal Properties"))
        self.first_pole_re_label_6.setText(_translate("Form", "Amplitude:"))
        self.square_duty_input.setText(_translate("Form", "0.5"))
        self.square_amplitude_input.setText(_translate("Form", "1"))
        self.first_pole_im_label_5.setText(_translate("Form", "Duty Cycle/Width (0,1]:"))
        self.sawtooth_plot_button.setText(_translate("Form", "Λ(t) Response"))
        self.groupBox_4.setTitle(_translate("Form", "Impulse Signal Input"))
        self.delta_plot_button.setText(_translate("Form", "𝛿(t) Response"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())
