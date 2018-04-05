# Misa, Bernadette B.
# CS 131 MP1

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import linalg, optimize
import pprint


def least_squares_cubic_poly(data):
    #n= no. of datapoints/length of col

    # x_arr = depth
    x_arr = data[0,:]
    # y_arr = temp
    y_arr = data[1,:]
    n_x = len(x_arr)
    n_y = len(y_arr)

    # print "x_arr", x_arr
    # print "y_arr", y_arr
    # print "n", n

    # Gets summation of the x, see derivation in the documentation
    least_square_matrix = []
    for i in range(0,n_x-1):
        if(i==0):
            least_square_matrix.append(float(n_x))
        else:
            x = sum(np.power(x_arr,i))
            least_square_matrix.append(x)

    # Makes it a jagged array
    lsm_final = []
    for i in range(0,4):
        lsm_final.append([])
        for j in range(i,4+i):

            lsm_final[i].append(least_square_matrix[j])
    # print "least_square_matrix final", lsm_final

    # Gets the proper y-vector
    y_vector = []
    for i in range(0, 4):
        if(i==0):
            y_vector.append(sum(y_arr))
        else:
            y_vector.append(np.matmul(y_arr,(np.power(x_arr,i))))
    # print "y_vector", y_vector

    # Makes it a jagged array
    y_vector_final = []
    for i in range(0,4):
            y_vector_final.append([y_vector[i]])
    # print "\ny_vector_final", y_vector_final

    return lsm_final

def LDLt(A_matrix):
    lsm_final = A_matrix
    data = np.array([[0.0, 2.3, 4.9, 9.1, 13.7, 18.3, 22.9, 27.2], [22.8, 22.8, 22.8, 20.6, 13.9, 11.7, 11.1, 11.1]], dtype=np.float64)
    x_arr = data[0,:]
    y_arr = data[1,:]
    n_x = len(x_arr)
    n_y = len(y_arr)

    # Gets the proper y-vector
    y_vector = []
    for i in range(0, 4):
        if(i==0):
            y_vector.append(sum(y_arr))
        else:
            y_vector.append(np.matmul(y_arr,(np.power(x_arr,i))))

    print "y_vector", y_vector

    # Makes it a jagged array
    y_vector_final = []
    for i in range(0,4):
            y_vector_final.append([y_vector[i]])
        
    print "\ny_vector_final", y_vector_final

    A = np.matrix(lsm_final)
    b = np.matrix(y_vector_final)
    print "\nA", A
    # print "\nb", b

    # Call LU function here next time
    # Ax=b
    # x = A/b
    x = np.linalg.solve(A,b)

    # print "\nx", x
    return x


def is_square(matrix):    
    matrix_row_len = len(matrix[0,:])
    matrix_col_len = len(matrix[:,0])
    # print "matrix_row_len: ", matrix_row_len
    # print "matrix_col_len: ", matrix_col_len
    if(matrix_row_len==matrix_col_len):
        return True
    else:
        return False

def generate_identity_matrix(matrix_size):
    identity_matrix = np.identity(matrix_size)
    identity_matrix = np.matrix(identity_matrix, dtype=np.float64)
    return identity_matrix

def LDLT(matrix):
    A = matrix
    matrix_size = len(A[0,:])
    L = generate_identity_matrix(matrix_size)
    D = generate_identity_matrix(matrix_size)

    i=0
    k=0
    while(i<matrix_size) :
        
        # # Generating the D matrix of D

        if(k==0):
            D[0,0] = A[0,0]
        elif(k==1):
            D[1,1] = A[1,1] - ((L[1,0]**2)*D[0,0])
        elif(k==2):
            D[2,2] = A[2,2] - ((L[2,0]**2)*D[0,0]) - ((L[2,1]**2)*D[1,1])
        elif(k==3):
            D[3,3] = A[3,3] - ((L[3,0]**2)*D[0,0]) - ((L[3,1]**2)*D[1,1]) -  ((L[4,2]**2)*D[2,2])
        else:
            print ""

        # print "D", D 
        # Generating the L matrix

        j=0
        while (j<matrix_size):
            print "\ni, j, matrix_size", i,j, matrix_size
            # Generating the first column of the L matrix
            if((i==0)):
                if(j==0):
                    L[j,i]=1.0

                L[j,i]= A[j,i]/D[i,i]
                print "A[" + str(j) + "," + str(j)+"]: ", A[j,i]
                print "D[i,i]", D[i,i]
                print "L[" + str(j) + "," + str(j)+"]: ", L[j,i]
                print "\n"


            # nth column 
            if(i == j):
                L[j,i]=1.0

            elif(i!=j):
                # L[j,i]= A[j,i]-L[j,i-1]*D[i-1,i-1]*L[i,i-1]/D[i,i]
                e=0
                val = 0
                while (e<i):
                    val= val + L[j,i-1]*D[i-1,i-1]*L[i,i-1]

                    L[j,i]= A[j,i]-val/D[i,i]
                    e=e+1

            else:
                print "else"
            j=j+1

        # ctr=ctr+1
        i=i+1
        k=k+1

        print "\nL", L
    LT = np.transpose(L).copy()

    return None  

def plot_graph(graph, number):

    item_no = number
    fcn = graph
    #     velocity_fcn_np = graph2
    #     acceleration_fcn_np = graph3

    plt.plot(fcn[0][:], fcn[1][:])
    
    if(item_no == 1):
        plt.title('Curve Fitting: Water Temperature (C) vs Depth (m)')
    
    plt.ylabel("Water Temperature (Celsius)")
    plt.xlabel("Depth (meter)")
    plt.show()


    #     # plt.plot(acceleration_fcn_np[0][0:21], acceleration_fcn_np[1][0:21])
    #     # plt.title('The acceleration function')
    #     # plt.show()

    #     fig, ax = plt.subplots()
    #     ax.plot(position_fcn_np[0][0:21], position_fcn_np[1][0:21], 'k--', color='b', label='Position (cm) vs Time (secs) Graph')
    #     ax.plot(velocity_fcn_np[0][0:21], velocity_fcn_np[1][0:21], 'k:', color='r', label='Velocity (cm/s) vs Time (secs) Graph')
    #     ax.plot(acceleration_fcn_np[0][0:21], acceleration_fcn_np[1][0:21], 'k', color='g', label='Acceleration (cm/s^2) vs Time (secs) Graph')

    #     legend = ax.legend(loc='best', shadow=True, fontsize='small')

    #     legend.get_frame().set_facecolor('#00FFCC')

    #     plt.show()


    # return x

def generate_Vandermonde_matrix(data):

    x_arr = data

    n = len(x_arr)

    V = generate_identity_matrix(n)
    # print "\nV before:", V
    
    for i in range(0, n):
        for j in range(0, n):
            V[i,j]=np.power(x_arr[i],j)
    # print "\nV after:",V
    # print "\nV after:"

    # pprint.pprint(V)
    return V

# Solves Ax=b using LU decomposition
def LU_decomp(A_matrix, b):
    # Ax=b
    # A=LU
    # LUx=b
    # Ux=y
    # Ly=b

    A = A_matrix
    P, L, U = scipy.linalg.lu(A)

    print "A:"
    pprint.pprint(A)

    # print "P:"
    # pprint.pprint(P)

    print "\L:"
    pprint.pprint(L)

    print "\U:"
    pprint.pprint(U)

    # Ly=b
    y = np.linalg.solve(L,b)
    # print "\ny", y

    # Ux=y
    x = np.linalg.solve(U,y)
    # print "\nx", x

    return x

def get_heat_flux(point):
    # heatflux = -alpha*density*specific_heat*point
    
    alpha = 0.1
    density = 1.0
    specific_heat = 1

    heatflux = -alpha*density*specific_heat*point
    # returns heatflux with unit cal/m^2*s

    # converting heatflux to unit cal/cm^2*day
    heatflux = heatflux * 24 * 60 * 60 *10000 

    return heatflux

def cubic_splines_wiki(C):
    # Inspired by the pseudocode in https://en.wikipedia.org/wiki/Spline_(mathematics)
    # C is set of coordinates

    # x_arr = depth
    x_arr = C[0,:]
    # print "x_arr", x_arr
    # y_arr = temp
    y_arr = C[1,:]
    # print "y_arr", y_arr

    n = len(x_arr)
    print "n", n

    a = y_arr.copy()
    b = np.zeros(n-1)
    d = np.zeros(n-1)
    h = np.zeros(n-1)

    print "h", h

    # Populates the step size vector
    for i in range(0, n-1):
        h[i] = x_arr[i+1] - x_arr[i]

    print "h", h

    alpha = np.zeros(n-1)

    for i in range(0, n-1):
        alpha[i] = (3/h[i])*(a[i+1]-a[i]) - (3/h[i-1])*(a[i]-a[i-1])

    print "alpha", alpha

    c = np.zeros(n)
    l = np.zeros(n)
    mu = np.zeros(n)
    z = np.zeros(n)

    print "c", c
    print "l", l
    print "mu", mu
    print "z", z

    # Setting up initial variables
    l[0] = 1.0
    mu[0] = 0.0
    z[0] = 0.0

    for i in range(1, n-1):
        l[i]=2*(x_arr[i+1] - x_arr[i-1]) - h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]
        z[i] = (alpha[i]-(h[i-1]*z[i-1]))/l[i]

    print "l", l
    print "mu", mu
    print "z", z

    l[-1] = 1.0
    z[-1] = 0.0
    c[-1] = 0.0

    print "\nl", l
    print "\nz", z
    print "\nc", c

    for j in range(n-2, 0, -1):
        c[j] = z[j] - mu[j]*c[j+1]
        b[j] = ((a[j+1]-a[j])/h[j]) - h[j]*((c[j+1]+2*c[j])/3)
        d[j] = (c[j+1]-c[j])/(3*h[j])

    print "\nc:", c
    print "\nb:", b
    print "\nd:", d

    output_set = []


    for i in range(0, n-1):
        S = [a[i], b[i], c[i], d[i], x_arr[i]]
        output_set.append(S)

    print "\noutput_set", output_set

    return output_set      
def main():

    # Given
    # # temp vs depth for Platte Lake
    temp_vs_depth = np.array([[0.0, 2.3, 4.9, 9.1, 13.7, 18.3, 22.9, 27.2], [22.8, 22.8, 22.8, 20.6, 13.9, 11.7, 11.1, 11.1]], dtype=np.float64)
    print temp_vs_depth

    # x_arr = depth
    x_arr = temp_vs_depth[0, :]
    # y_arr = temp
    y_arr = temp_vs_depth[1, :]
    n = len(x_arr)    

    # ****************************************************************************************
    # ****************************************************************************************
    # 1. Curve Fitting
    # ****************************************************************************************
    # ****************************************************************************************
    
    # 1.1 Performing least squares regression to a cubic polynomial
    # ========================================================================================
    print "\n 1.1 Performing least squares regression to a cubic polynomial"
    poly_fcn = least_squares_cubic_poly(temp_vs_depth)
    print "\nLeast squares regression matrix:", poly_fcn
    # print "\nnew_y", new_y
    # ========================================================================================


    # 1.2 LDLT
    # ========================================================================================
    print "\n 1.2 Determining a0 - a3 using LDLT decomposition\n"

    curve_fitting_model = LDLt(poly_fcn)
    print "\na0 - a3 coefficients: \n", curve_fitting_model
    # curve_fitting = yung result ng ldlt
    # ========================================================================================

    # 1.4 Plotting the Water Temperature vs. Depth graph using the Curve Fitting Model
    # ========================================================================================
    print "\nPlotting the Water Temperature vs. Depth graph using the Curve Fitting Model"
    
    # plugging in datapoints(depth) in the curve fitting model to generate new_y
    new_y = []
    for i in range(0, n):
        j = np.polyval(curve_fitting_model, x_arr[i])
        new_y.append(float(j))
    # print "\nnew_y", new_y
    curve_fitting = [x_arr,new_y]
    # print "\ncurve_fitting", curve_fitting
    plot_graph(curve_fitting, 1)
    # ========================================================================================

    
    # 1.5 Solving for the gradient
    # ========================================================================================
    print "\n 1.5 Solving for the gradient"

    # poly_fcn_grad = [poly_fcn[0, 0], poly_fcn[1, 0],poly_fcn[2, 0], poly_fcn[3, 0]]
    curve_fitting_gradient = [curve_fitting_model[0, 0], curve_fitting_model[1, 0],curve_fitting_model[2, 0], curve_fitting_model[3, 0]]
    
    # curve_fitting_gradient = [curve_fitting_model[0,0]]

    curve_fitting_gradient = np.polyder(curve_fitting_gradient)
    print "\ncurve_fitting_gradient:", curve_fitting_gradient
    
    curve_fitting_gradient_abs = abs(curve_fitting_gradient)
    # print "\nabs value of curve_fitting_gradient", curve_fitting_gradient_abs 
    # print "abs value of curve_fitting_gradient, (max elem:)", max(curve_fitting_gradient_abs)
    # print "abs value of curve_fitting_gradient, (max elem index:)", np.argmax(curve_fitting_gradient_abs)
    # ========================================================================================

    # 1.3 Solve for the Inflection Point (thermocline)
    # ========================================================================================
    # Second derivative of f
    print "\n 1.3. Solving for the infleciton point"

    curve_fitting_thermocline = np.polyder(curve_fitting_gradient)
    print "\ncurve_fitting_thermocline", curve_fitting_thermocline
    curve_fitting_thermocline_roots = np.roots(curve_fitting_thermocline)
    print "curve_fitting_thermocline", curve_fitting_thermocline_roots
    # ========================================================================================


    # 1.6 Plotting Gradient vs Lake Depth
    # ========================================================================================
    print "\n 1.6 Plotting Gradient vs Lake Depth"

    # plugin datapoints(depth) sa gradient polynomial to generate new_y
    # new_y = np.polyval(curve_fitting, iterate datapoints here)
    # plot_graph(curve_fitting)
    # ========================================================================================


    # 1.7 Determininng the heat flux across the thermocline using the gradient and using Fourier's Law
    # ========================================================================================
    print "\n 1.7 Computing heatflux at the thermocline"
    curve_fitting_heatflux = get_heat_flux(curve_fitting_thermocline)
    # print "curve_fitting_heatflux", curve_fitting_heatflux
    # ========================================================================================

    
    # ****************************************************************************************
    # ****************************************************************************************
    # 2. Monomial Basis Interpolation
    # ****************************************************************************************
    # ****************************************************************************************
    
    # 2.1 Constructing a Vandermode Matrixx
    # ========================================================================================
    
    print "\n 2.1 Generating the Vandermonde matrix...\n"
    V = generate_Vandermonde_matrix(x_arr)
    V = np.fliplr(V)
    V = np.flipud(V)
    print "V", V
    # ========================================================================================

    # 2.2 Determining a0-a7 using LU Decomposition
    # ========================================================================================
    # using built-in LU decomp
    print "\n 2.2 Determining a0 - a7 using LU decomposition\n"
    monomial_poly_fcn = LU_decomp(V, y_arr)
    # Ly=b
    w = np.linalg.solve(V, y_arr)
    print "\ntry", w
    
    # print "\ny", y
    # print "monomial_poly_fcn:", monomial_poly_fcn

    # ========================================================================================

    # 2.5 Solving for the gradient
    # ========================================================================================
    print "\n 2.5 Solving for the gradient"

    monomial_poly_fcn = [monomial_poly_fcn[0], monomial_poly_fcn[1],monomial_poly_fcn[2], monomial_poly_fcn[3], monomial_poly_fcn[4], monomial_poly_fcn[5], monomial_poly_fcn[6], monomial_poly_fcn[7]]
    # monomial_poly_fcn = np.poly1d(monomial_poly_fcn)
    print "\nmonomial poly_fcn:", monomial_poly_fcn
    monomial_gradient = np.polyder(monomial_poly_fcn)
    print "\nmonomial gradient:", monomial_gradient

    monomial_gradient_abs = abs(monomial_gradient)
    print "\nabs value of monomial gradient", monomial_gradient_abs 
    print "abs value of monomial gradient, (max elem:)", max(monomial_gradient_abs )
    max_elem_index = np.argmax(monomial_gradient_abs)
    print "abs value of monomial gradient, (max elem index:)", max_elem_index
    # ========================================================================================

    # 2.3 Solve for the inflection point
    # ========================================================================================
    print "\n 2.3 Solving for the inflection point"

    monomial_thermocline = np.polyder(monomial_poly_fcn)
    print "\nmonomial thermocline", monomial_thermocline
    monomial_thermocline = np.roots(monomial_thermocline)
    print "\nmonomial thermocline roots", monomial_thermocline
    # monomial_thermocline = scipy.optimize.fsolve(monomial_thermocline, 0)

    # monomial_thermocline = monomial_thermocline[max_elem_index]
    # print "\nmonomial thermocline max", monomial_thermocline
    # ========================================================================================


    # 2.4 Plotting the Water Temperature vs. Depth graph using Monomial Basis Interpolation
    # ========================================================================================
    print "\n 2.4 Plotting the Water Temperature vs Depth graph"
    # ========================================================================================


    # 2.6 Determininng the heat flux across the thermocline using the gradient and using Fourier's Law
    # ========================================================================================
    print "\n 2.6 Computing for the heatflux"
    monomial_heatflux = get_heat_flux(monomial_thermocline)
    print "\nmonomial_heatflux", monomial_heatflux
    # ========================================================================================

    # 2.7 Plotting the Gradient vs Lake Depth
    # ========================================================================================


    # ========================================================================================

    # ****************************************************************************************
    # ****************************************************************************************
    # 3. Piecewise Cubic Polynomial Interpolation
    # ****************************************************************************************
    # ****************************************************************************************

    # 3.1 Polynomial Cubic Splines Modelling  
    # ========================================================================================
    print "\n3.1 Cubic Splines Setup"
    print "\nSolving for the cubic splines:"
    cubic_splines_wiki(temp_vs_depth)
    # ========================================================================================

    # 3.2 Determining the depth of the thermocline
    # ========================================================================================
    print "\n 3.2 Determining the depth of the thermocline"

    # ========================================================================================
    
    # 3.3 Solving for the gradient
    # ========================================================================================
    print "\n 3.3 Solving for the gradient"

    # ========================================================================================
    
    # 3.5 Plotting the water temperature vs lake depth
    # ========================================================================================
    print "\n 3.5 Plotting the water temperature vs lake depth"

    # ========================================================================================

    # 3.6 Plotting the gradient vs lake depth
    # ========================================================================================
    print "\n 3.6 Plotting the gradient vs lake depth"

    # ================================================
if __name__ == "__main__":
    main()



