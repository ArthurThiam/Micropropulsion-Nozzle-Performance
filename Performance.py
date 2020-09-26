from numpy import *
import plotly.express as px
import plotly.graph_objects as go
from scipy import signal
from PyQt5 import QtCore, QtGui, QtWidgets
import pandas

# =====================================================================================================================
# Area ratio calculator
def pressure_ratio(epsilon_true, gamma, Gamma):

    pressure_ratio_calculated = 0.001
    pressure_stepsize = 0.00005
    epsilon_calculated = 0

    error = abs(epsilon_true - epsilon_calculated)
    threshold = 1

    while error > threshold:
        epsilon_calculated = Gamma / sqrt((2 * gamma / (gamma - 1)) * pressure_ratio_calculated ** (2 / gamma) * (
                    1 - pressure_ratio_calculated ** ((gamma - 1) / gamma)))
        error = abs(epsilon_true - epsilon_calculated)

        pressure_ratio_calculated += pressure_stepsize

    return pressure_ratio_calculated


# Performance Calculator
def determine_performance(alpha, w_t, p_c):

    # Import Constant Parameters
    gamma = ui.gamma.value()
    p_e = ui.p_e.value()
    R = ui.R_s.value()
    T_c = ui.T_c.value()
    mu = ui.mu.value()*10**(-6)
    rho_c = ui.rho_c.value()
    T_w = ui.T_w.value()
    d_n = ui.D_n.value()
    epsilon = ui.epsilon.value()

    # Calculate Flow Parameters
    Gamma = sqrt(gamma) * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
    rho_e = rho_c * (p_e / p_c) ** (1 / gamma)
    u_e = sqrt(2 * gamma / (gamma - 1) * (R) * T_c * (1 - (p_e / p_c) ** ((gamma - 1) / gamma)))
    M_e = u_e / sqrt(gamma * R * T_c)
    T_e = T_c * (p_e / p_c) ** ((gamma - 1) / gamma)
    T_total = T_e * (1 + ((gamma - 1) / 2) * M_e ** 2)

    # Calculate Nozzle Parameters
    A_t = w_t * d_n
    d_h = 4 * A_t / (2 * w_t + 2 * d_n)
    r_exit = epsilon * w_t / 2
    l_divergent = r_exit / sin(alpha * pi / 180)

    variable_parameters = [w_t, alpha, p_c]  # For data logging

    # ======================================================================================================================
    # IDEAL PERFORMANCE

    m_ideal = (Gamma * p_c * A_t) / sqrt(R * T_c)
    F_ideal = m_ideal * sqrt((2 * gamma / (gamma - 1)) * R * T_c * (1 - (p_e / p_c) ** ((gamma - 1) / gamma)))
    C_F_ideal = F_ideal / (p_c * A_t)
    Re_t = (m_ideal * d_h) / (mu * A_t)

    # ======================================================================================================================
    # LOSS CALCULATOR

    divergence_loss = 1 - 0.5 * (1 - cos((pi * alpha) / 180))

    c_f = 0.664 / sqrt(Re_t)
    compressible_c_f = c_f * (((T_w / T_total) + 1) / 2 + 0.22 * (gamma - 1) / 2 * M_e ** 2) ** (-0.6)
    momentum_thickness = compressible_c_f * l_divergent  # at exit
    displacement_thickness = 2.59036 * momentum_thickness

    delta_F_momentum = rho_e * u_e * d_h * momentum_thickness * u_e
    epsilon_true = 2 * (r_exit - displacement_thickness) / w_t
    true_pressure_ratio = pressure_ratio(epsilon_true, gamma, Gamma)
    p_e_true = true_pressure_ratio * p_c

    # ==================================================================================================================
    # TRUE PERFORMANCE

    F_ideal_new = m_ideal * sqrt((2 * gamma / (gamma - 1)) * R * T_c * (1 - (p_e_true / p_c) ** ((gamma - 1) / gamma)))

    F_analytical = F_ideal_new * divergence_loss - delta_F_momentum
    n_F = F_analytical / F_ideal
    I_sp = F_analytical / (m_ideal * 9.81)

    performance_parameters = [I_sp, n_F, Re_t, divergence_loss, delta_F_momentum, epsilon_true]  # For data logging

    # ==================================================================================================================
    # DATA LOGGING

    return [variable_parameters, performance_parameters]


# Plot function
def plot_microprop_performance():
    print(ui.filter_status)

    # Import iteration settings
    alpha_initial = ui.alpha_initial.value()
    w_t_initial = ui.w_t_initial.value()*10**(-6)
    p_c_initial = ui.p_c_initial.value()

    n_alpha = ui.n_alpha.value()
    n_width = ui.n_width.value()
    n_pressure = ui.n_pressure.value()

    delta_alpha = ui.delta_alpha.value()
    delta_width = ui.delta_width.value()*10**(-6)
    delta_pressure = ui.delta_pressure.value()

    alpha_range = arange(alpha_initial, alpha_initial + n_alpha * delta_alpha, delta_alpha)
    # print('alpha values: ', alpha_range)

    width_range = arange(w_t_initial, round(w_t_initial + n_width * delta_width, 6), delta_width)
    # print('throat width values: ', width_range)

    pressure_range = arange(p_c_initial, p_c_initial + n_pressure * delta_pressure, delta_pressure)
    # print('chamber pressure values: ', pressure_range)
    # print('')

    loop = [alpha_range, width_range, pressure_range]  # Matrix containing all variable parameter values to calculate
    # performance for.

    VP_1 = ui.VP_1.value()
    VP_2 = ui.VP_2.value()
    VP_3 = ui.VP_3.value()
    PP = ui.PP.value()

    # Determine performance for all variable parameter configurations
    cut_off = len(loop[VP_1])  # Length of data for each curve
    x_val = []
    curve_val = []
    plot_val = []
    y_val = []
    total_output = []

    # Exit Condition
    # TBD

    print('Starting computations...')
    # print(loop[VP_1])
    # Performance Calculation
    for a in (loop[VP_3]):  # Loop for different plots
        for b in (loop[VP_2]):  # Loop for different curves
            for c in (loop[VP_1]):  # Loop for different points

                if (VP_1 == 0) & (VP_2 == 1):
                    output = determine_performance(c, b, a)

                elif (VP_1 == 0) & (VP_2 == 2):
                    output = determine_performance(c, a, b)

                elif (VP_1 == 1) & (VP_2 == 0):
                    output = determine_performance(b, c, a)

                elif (VP_1 == 1) & (VP_2 == 2):
                    output = determine_performance(a, c, b)

                elif (VP_1 == 2) & (VP_2 == 0):
                    output = determine_performance(b, a, c)

                elif (VP_1 == 2) & (VP_2 == 1):
                    output = determine_performance(a, b, c)

                else:
                    print('Error in performance calculation.')
                    sys.exit(1)

                curve_val.append(output[0][VP_2])
                plot_val.append(output[0][VP_3])

                total_output.append(output)

                y_val.append(output[1][PP])

    print('Computations Finished.')
    print('')
    print('Plotting data...')
    print('')
    # ======================================================================================================================
    # DATA PLOTTING

    n_curve = 0  # Curve counter
    n_plot = 0  # plot counter

    # Define output parameters
    x_labels = ['exit angle [°]', 'nozzle width [um]', 'chamber pressure [Pa]']
    y_labels = ['Specific Impulse [s]', 'Thrust Efficiency', 'Reynolds Number', 'Divergence Loss', 'Momentum loss',
                'Area ratio efficiency']
    curve_labels = {}

    # Loop running through each plot
    while n_plot < len(loop[VP_3]):

        # Create figure
        fig = go.Figure()
        graph_title = str(y_labels[PP]) + ' as function of ' + str(x_labels[VP_1]) + ' for ' + str(
            x_labels[VP_3]) + " = " + str(loop[VP_3][n_plot])

        # Loop running through each curve

        while n_curve < len(loop[VP_2]):
            # print('plotting curve for n_plot = ', n_plot, ', n_curve = ', n_curve)

            # extract data for current curve
            # print('x values: ', loop[VP_1])
            # print('y_values: ', y_val)
            # print('')

            y_out = y_val[0:cut_off]

            # resize data lists
            x_val = x_val[cut_off:]
            y_val = y_val[cut_off:]

            if ui.filter_bool.isChecked():
                fig.add_trace(go.Scatter(x=loop[VP_1], y=signal.savgol_filter(y_out, 53, 3),
                                     name=str(x_labels[VP_2]) + " = " + str(loop[VP_2][n_curve])))

            elif not ui.filter_bool.isChecked():
                fig.add_trace(go.Scatter(x=loop[VP_1], y=y_out,
                                     name=str(x_labels[VP_2]) + " = " + str(loop[VP_2][n_curve])))

            n_curve += 1

        # Format figure
        fig.update_layout(title=graph_title, xaxis_title=x_labels[VP_1], yaxis_title=y_labels[PP])
        fig.show()

        n_plot += 1
        n_curve = 0

    print('')
    print('Data plotting complete.')

# UI class
class Ui_MainWindow(object):

    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1073, 843)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.groupBox = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox.setGeometry(QtCore.QRect(10, 350, 361, 251))
        self.groupBox.setObjectName("groupBox")
        self.label = QtWidgets.QLabel(self.groupBox)
        self.label.setGeometry(QtCore.QRect(20, 30, 47, 14))
        self.label.setObjectName("label")
        self.gamma = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.gamma.setGeometry(QtCore.QRect(210, 30, 62, 22))
        self.gamma.setMaximum(2.0)
        self.gamma.setProperty("value", 1.4)
        self.gamma.setObjectName("gamma")
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.label_2.setGeometry(QtCore.QRect(20, 60, 141, 16))
        self.label_2.setObjectName("label_2")
        self.p_e = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.p_e.setGeometry(QtCore.QRect(210, 60, 62, 22))
        self.p_e.setProperty("value", 30.0)
        self.p_e.setObjectName("p_e")
        self.label_3 = QtWidgets.QLabel(self.groupBox)
        self.label_3.setGeometry(QtCore.QRect(20, 90, 171, 16))
        self.label_3.setObjectName("label_3")
        self.R_s = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.R_s.setGeometry(QtCore.QRect(210, 90, 71, 22))
        self.R_s.setMaximum(500.0)
        self.R_s.setProperty("value", 297.0)
        self.R_s.setObjectName("R_s")
        self.label_4 = QtWidgets.QLabel(self.groupBox)
        self.label_4.setGeometry(QtCore.QRect(20, 120, 171, 16))
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(self.groupBox)
        self.label_5.setGeometry(QtCore.QRect(20, 150, 161, 16))
        self.label_5.setObjectName("label_5")
        self.rho_c = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.rho_c.setGeometry(QtCore.QRect(210, 150, 101, 22))
        self.rho_c.setDecimals(4)
        self.rho_c.setProperty("value", 3.3648)
        self.rho_c.setObjectName("rho_c")
        self.label_6 = QtWidgets.QLabel(self.groupBox)
        self.label_6.setGeometry(QtCore.QRect(20, 180, 161, 16))
        self.label_6.setObjectName("label_6")
        self.T_c = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.T_c.setGeometry(QtCore.QRect(210, 180, 71, 22))
        self.T_c.setMaximum(500.0)
        self.T_c.setProperty("value", 300.55)
        self.T_c.setObjectName("T_c")
        self.mu = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.mu.setGeometry(QtCore.QRect(210, 120, 101, 22))
        self.mu.setDecimals(2)
        self.mu.setProperty("value", 17.86)
        self.mu.setObjectName("mu")
        self.label_7 = QtWidgets.QLabel(self.groupBox)
        self.label_7.setGeometry(QtCore.QRect(20, 210, 131, 16))
        self.label_7.setObjectName("label_7")
        self.T_w = QtWidgets.QDoubleSpinBox(self.groupBox)
        self.T_w.setGeometry(QtCore.QRect(210, 210, 81, 22))
        self.T_w.setMaximum(500.0)
        self.T_w.setProperty("value", 160.0)
        self.T_w.setObjectName("T_w")
        self.label_8 = QtWidgets.QLabel(self.groupBox)
        self.label_8.setGeometry(QtCore.QRect(320, 120, 31, 16))
        self.label_8.setObjectName("label_8")
        self.groupBox_2 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_2.setGeometry(QtCore.QRect(10, 620, 291, 101))
        self.groupBox_2.setObjectName("groupBox_2")
        self.label_13 = QtWidgets.QLabel(self.groupBox_2)
        self.label_13.setGeometry(QtCore.QRect(20, 30, 101, 16))
        self.label_13.setObjectName("label_13")
        self.epsilon = QtWidgets.QDoubleSpinBox(self.groupBox_2)
        self.epsilon.setGeometry(QtCore.QRect(190, 30, 62, 22))
        self.epsilon.setDecimals(3)
        self.epsilon.setProperty("value", 16.971)
        self.epsilon.setObjectName("epsilon")
        self.label_14 = QtWidgets.QLabel(self.groupBox_2)
        self.label_14.setGeometry(QtCore.QRect(20, 60, 111, 16))
        self.label_14.setObjectName("label_14")
        self.D_n = QtWidgets.QDoubleSpinBox(self.groupBox_2)
        self.D_n.setGeometry(QtCore.QRect(190, 60, 81, 22))
        self.D_n.setDecimals(4)
        self.D_n.setProperty("value", 0.0001)
        self.D_n.setObjectName("D_n")
        self.groupBox_3 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_3.setGeometry(QtCore.QRect(410, 350, 311, 371))
        self.groupBox_3.setObjectName("groupBox_3")
        self.label_15 = QtWidgets.QLabel(self.groupBox_3)
        self.label_15.setGeometry(QtCore.QRect(20, 30, 151, 16))
        self.label_15.setObjectName("label_15")
        self.w_t_initial = QtWidgets.QDoubleSpinBox(self.groupBox_3)
        self.w_t_initial.setGeometry(QtCore.QRect(220, 30, 81, 22))
        self.w_t_initial.setDecimals(2)
        self.w_t_initial.setProperty("value", 45.0)
        self.w_t_initial.setObjectName("w_t_initial")
        self.label_16 = QtWidgets.QLabel(self.groupBox_3)
        self.label_16.setGeometry(QtCore.QRect(20, 60, 141, 16))
        self.label_16.setObjectName("label_16")
        self.delta_width = QtWidgets.QDoubleSpinBox(self.groupBox_3)
        self.delta_width.setGeometry(QtCore.QRect(220, 60, 81, 22))
        self.delta_width.setDecimals(0)
        self.delta_width.setProperty("value", 5.0)
        self.delta_width.setObjectName("delta_width")
        self.label_17 = QtWidgets.QLabel(self.groupBox_3)
        self.label_17.setGeometry(QtCore.QRect(20, 90, 121, 16))
        self.label_17.setObjectName("label_17")
        self.label_18 = QtWidgets.QLabel(self.groupBox_3)
        self.label_18.setGeometry(QtCore.QRect(20, 180, 141, 16))
        self.label_18.setObjectName("label_18")
        self.delta_alpha = QtWidgets.QDoubleSpinBox(self.groupBox_3)
        self.delta_alpha.setGeometry(QtCore.QRect(220, 180, 62, 22))
        self.delta_alpha.setProperty("value", 0.1)
        self.delta_alpha.setObjectName("delta_alpha")
        self.label_19 = QtWidgets.QLabel(self.groupBox_3)
        self.label_19.setGeometry(QtCore.QRect(20, 210, 141, 16))
        self.label_19.setObjectName("label_19")
        self.label_20 = QtWidgets.QLabel(self.groupBox_3)
        self.label_20.setGeometry(QtCore.QRect(20, 150, 131, 16))
        self.label_20.setObjectName("label_20")
        self.alpha_initial = QtWidgets.QDoubleSpinBox(self.groupBox_3)
        self.alpha_initial.setGeometry(QtCore.QRect(220, 150, 62, 22))
        self.alpha_initial.setProperty("value", 10.0)
        self.alpha_initial.setObjectName("alpha_initial")
        self.label_21 = QtWidgets.QLabel(self.groupBox_3)
        self.label_21.setGeometry(QtCore.QRect(20, 300, 201, 16))
        self.label_21.setObjectName("label_21")
        self.label_22 = QtWidgets.QLabel(self.groupBox_3)
        self.label_22.setGeometry(QtCore.QRect(20, 330, 141, 16))
        self.label_22.setObjectName("label_22")
        self.p_c_initial = QtWidgets.QDoubleSpinBox(self.groupBox_3)
        self.p_c_initial.setGeometry(QtCore.QRect(220, 270, 81, 22))
        self.p_c_initial.setMaximum(100000.0)
        self.p_c_initial.setProperty("value", 80000.0)
        self.p_c_initial.setObjectName("p_c_initial")
        self.delta_pressure = QtWidgets.QDoubleSpinBox(self.groupBox_3)
        self.delta_pressure.setGeometry(QtCore.QRect(220, 300, 81, 22))
        self.delta_pressure.setMaximum(10000000.0)
        self.delta_pressure.setProperty("value", 10000.0)
        self.delta_pressure.setObjectName("delta_pressure")
        self.label_23 = QtWidgets.QLabel(self.groupBox_3)
        self.label_23.setGeometry(QtCore.QRect(20, 270, 171, 16))
        self.label_23.setObjectName("label_23")
        self.n_width = QtWidgets.QSpinBox(self.groupBox_3)
        self.n_width.setGeometry(QtCore.QRect(220, 90, 43, 24))
        self.n_width.setProperty("value", 2)
        self.n_width.setDisplayIntegerBase(5)
        self.n_width.setObjectName("n_width")
        self.n_pressure = QtWidgets.QSpinBox(self.groupBox_3)
        self.n_pressure.setGeometry(QtCore.QRect(220, 330, 43, 24))
        self.n_pressure.setProperty("value", 5)
        self.n_pressure.setObjectName("n_pressure")
        self.n_alpha = QtWidgets.QSpinBox(self.groupBox_3)
        self.n_alpha.setGeometry(QtCore.QRect(220, 210, 51, 24))
        self.n_alpha.setMaximum(200)
        self.n_alpha.setProperty("value", 100)
        self.n_alpha.setDisplayIntegerBase(10)
        self.n_alpha.setObjectName("n_alpha")

        self.groupBox_4 = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_4.setGeometry(QtCore.QRect(760, 520, 291, 201))
        self.groupBox_4.setObjectName("groupBox_4")
        self.label_24 = QtWidgets.QLabel(self.groupBox_4)
        self.label_24.setGeometry(QtCore.QRect(10, 30, 111, 16))
        self.label_24.setObjectName("label_24")
        self.label_26 = QtWidgets.QLabel(self.groupBox_4)
        self.label_26.setGeometry(QtCore.QRect(10, 110, 111, 16))
        self.label_26.setObjectName("label_26")
        self.label_27 = QtWidgets.QLabel(self.groupBox_4)
        self.label_27.setGeometry(QtCore.QRect(10, 70, 111, 16))
        self.label_27.setObjectName("label_27")
        self.VP_1 = QtWidgets.QSpinBox(self.groupBox_4)
        self.VP_1.setGeometry(QtCore.QRect(120, 30, 43, 24))
        self.VP_1.setMaximum(3)
        self.VP_1.setObjectName("VP_1")
        self.VP_2 = QtWidgets.QSpinBox(self.groupBox_4)
        self.VP_2.setGeometry(QtCore.QRect(120, 70, 43, 24))
        self.VP_2.setMaximum(3)
        self.VP_2.setProperty("value", 2)
        self.VP_2.setObjectName("VP_2")
        self.VP_3 = QtWidgets.QSpinBox(self.groupBox_4)
        self.VP_3.setGeometry(QtCore.QRect(120, 110, 43, 24))
        self.VP_3.setMaximum(3)
        self.VP_3.setProperty("value", 1)
        self.VP_3.setObjectName("VP_3")
        self.label_25 = QtWidgets.QLabel(self.groupBox_4)
        self.label_25.setGeometry(QtCore.QRect(10, 160, 111, 16))
        self.label_25.setObjectName("label_25")
        self.PP = QtWidgets.QSpinBox(self.groupBox_4)
        self.PP.setGeometry(QtCore.QRect(120, 160, 43, 24))
        self.PP.setProperty("value", 1)
        self.PP.setObjectName("PP")
        self.line = QtWidgets.QFrame(self.groupBox_4)
        self.line.setGeometry(QtCore.QRect(10, 140, 251, 16))
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")

        # Filter Boolean
        self.filter_bool = QtWidgets.QCheckBox(self.groupBox_4)
        self.filter_bool.setGeometry(QtCore.QRect(180, 160, 101, 21))
        self.filter_bool.setObjectName("filter_bool")
        if self.filter_bool.isChecked:
            self.filter_status = True
        elif not self.filter_bool.isChecked:
            pass



        self.textBrowser = QtWidgets.QTextBrowser(self.centralwidget)
        self.textBrowser.setGeometry(QtCore.QRect(10, 20, 1031, 311))
        self.textBrowser.setObjectName("textBrowser")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1073, 20))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        # Plot Button
        self.plot_performance = QtWidgets.QPushButton(self.centralwidget)
        self.plot_performance.setGeometry(QtCore.QRect(870, 740, 171, 51))
        self.plot_performance.setObjectName("plot_performance")
        self.plot_performance.clicked.connect(plot_microprop_performance)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.groupBox.setTitle(_translate("MainWindow", "Flow Properties"))
        self.label.setText(_translate("MainWindow", "gamma"))
        self.label_2.setText(_translate("MainWindow", "exit pressure [Pa]"))
        self.label_3.setText(_translate("MainWindow", "specific gas constant [J/kgK]"))
        self.label_4.setText(_translate("MainWindow", "dynamic viscosity [Pa*s]"))
        self.label_5.setText(_translate("MainWindow", "chamber density [kg/m3]"))
        self.label_6.setText(_translate("MainWindow", "chamber temperature [K]"))
        self.label_7.setText(_translate("MainWindow", "wall temperature [K]"))
        self.label_8.setText(_translate("MainWindow", "e-6"))
        self.groupBox_2.setTitle(_translate("MainWindow", "Constant nozzle properties"))
        self.label_13.setText(_translate("MainWindow", "expansion ratio"))
        self.label_14.setText(_translate("MainWindow", "throat depth [m]"))
        self.groupBox_3.setTitle(_translate("MainWindow", "Variable nozzle properties"))
        self.label_15.setText(_translate("MainWindow", "initial throat width [um]"))
        self.label_16.setText(_translate("MainWindow", "throat width step [um]"))
        self.label_17.setText(_translate("MainWindow", "number of widths"))
        self.label_18.setText(_translate("MainWindow", "exit angle step size [°]"))
        self.label_19.setText(_translate("MainWindow", "number of exit angles"))
        self.label_20.setText(_translate("MainWindow", "initial exit angle [°]"))
        self.label_21.setText(_translate("MainWindow", "chamber pressure step size [Pa]"))
        self.label_22.setText(_translate("MainWindow", "number of pressures"))
        self.label_23.setText(_translate("MainWindow", "initial chamber pressure [Pa]"))
        self.plot_performance.setText(_translate("MainWindow", "Plot Performance"))
        self.groupBox_4.setTitle(_translate("MainWindow", "Plot Options"))
        self.label_24.setText(_translate("MainWindow", "x-axis variable"))
        self.label_26.setText(_translate("MainWindow", "plot variable"))
        self.label_27.setText(_translate("MainWindow", "curve variable"))
        self.label_25.setText(_translate("MainWindow", "y-axis variable"))
        self.filter_bool.setText(_translate("MainWindow", "Filter output"))
        self.textBrowser.setHtml(_translate("MainWindow",
                                            "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
                                            "<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
                                            "p, li { white-space: pre-wrap; }\n"
                                            "</style></head><body style=\" font-family:\'Sans Serif\'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600; text-decoration: underline;\">INSTRUCTIONS</span></p>\n"
                                            "<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">&quot;Flow properties&quot; and &quot;Constant nozzle properties&quot; will be kept constant for all plots. The variable nozzle properties determine what ranges of values are applied to the variable nozzle parameters. The plot parameters determine how values are plotted across different plots. the curve-variable will be changed for each durve in a plot, while the plot variable will be altered only between plots. </p>\n"
                                            "<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">x-axis, curve and plot variables are chosen from the following list:</p>\n"
                                            "<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">0 = exit angle</p>\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">1 = throat width</p>\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">2 = chamber pressure</p>\n"
                                            "<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">The y-axis variable can be chosen from the following options:</p>\n"
                                            "<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">0 = Specific Impulse</p>\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">1 = Thrust efficiency</p>\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">2 = Reynolds number</p>\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">3 = Divergence Loss</p>\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">4 = Momentum Loss</p>\n"
                                            "<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">5 = True/ideal exit area ratio</p></body></html>"))

# Main
if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
