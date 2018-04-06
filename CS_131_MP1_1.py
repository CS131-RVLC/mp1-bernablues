# Misa, Bernadette B.
# CS 131 MP1

import matplotlib.pyplot as plt
import numpy as np
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
    y_vector_final = np.flipud(y_vector_final)
    # y_vector_final = np.fliplr(y_vector_final)
    # print "\ny_vector_final", y_vector_final

    return lsm_final, y_vector_final

# Ly=b
# Solves for y in Ly = b, where L is a lower triangular matrix
def forward_sub(L,b):
    y = []
    for i in range(len(b)):
        y.append(b[i])
        for j in range(i):
            y[i] = y[i] - (L[i,j] * y[j])
        y[i] = y[i] / L[i,i]
    return y

# Ux=y
# Solves for x in Ux = y, where U is an upper triangular matrix
def backward_sub(U,y):
    x = []
    for i in range(len(y)):
        x.append(float(0))

    for i in range(0,len(y)):
    # for i in range(1,len(y)):
        x[-i] = (y[-i])
        for j in range(1,i):
            x[-i] = x[-i] - (U[-i][-j]* x[-j])
        
        # print "i, U[i,i]", i, U[i][i]

        x[-i] = x[-i] / U[-i][-i]


    return x

# Solves for x in Dx = b, where D is a diagonal matrix
def diagonal_sub(D,b):
    x = []
    for i in range(len(b)):
        temp = b[i] / D[i,i]
        x.append(temp)
    return x

def LDLt(A_matrix, b):
    A = A_matrix
    n = len(A)
    L = generate_identity_matrix(n)
    U = generate_identity_matrix(n)
    
    # Generating the U matrix
    for i in range(0,n):
        U[0,i] = A[0][i]

    # Generating the L matrix
    for i in range(0,n):
        L[i,i] = float(1)
        for j in range(i+1,n):
            L[i,j] = float(0)
        if i==0:
            for j in range(1,n):
                L[j,0] = A[j,0] / U[0,0]

    U[1,1] = A[1,1] - U[0,1] * L[1,0]
    L[2,1] = (A[2,1] - U[0,1] * L[2,0]) / U[1,1]
    L[3,1] = (A[3,1] - U[0,1] * L[3,0]) / U[1,1]

    U[1,2] = A[1,2] - U[0,2] * L[1,0]
    U[2,2] = A[2,2] - U[0,2] * L[2,0] - U[1,2] * L[2,1]
    L[3,2] = (A[3,2] - U[0,2] * L[3,0] - U[1,2] * L[3,1]) / U[2,2]

    U[1,3] = A[1,3] - U[0,3] * L[1,0]
    U[2,3] = A[2,3] - U[0,3] * L[2,0] - U[1,3] * L[2,1]
    U[3,3] = A[3,3] - U[0,3] * L[3,0] - U[1,3] * L[3,1] - U[2,3] * L[3,2]

    # Generating the LT matrix
    LT = []
    for i in range(0,n):
        LT.append([])
        for j in range(n):
            LT[-1].append(L[j,i])

    # Generating the D matrix
    D = generate_identity_matrix(n)
    for i in range(0,len(U)):
        D[i,i] = U[i,i]

    print "\nL", L
    print "\nD", D
    print "\nLT", LT

    # Ly=b
    y = forward_sub(L,b)
    # print "\ny", y

    # Dy=z    
    z = diagonal_sub(D,y)
    # print "\nz", z
    
    # LTx=z
    x = backward_sub(LT,z)
    # print "x", x

    return x

def generate_identity_matrix(matrix_size):
    identity_matrix = np.identity(matrix_size)
    identity_matrix = np.matrix(identity_matrix, dtype=np.float64)
    return identity_matrix

def plot_graph(graph, number):

    item_no = number
    fcn = graph

    plt.plot(fcn[0][:], fcn[1][:])
    
    if(item_no == 1):
        plt.title('Curve Fitting: Water Temperature (C) vs Depth (m)')
        plt.ylabel("Water Temperature (Celsius)")
    if(item_no == 2):
        plt.title('Curve Fitting: Gradient (dT/dz) vs Depth (m)')
        plt.ylabel("Gradient (dT/dz)")
        
    plt.xlabel("Depth (meter)")
    plt.show()



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
    poly_fcn, y_new = least_squares_cubic_poly(temp_vs_depth)
    poly_fcn = np.fliplr(poly_fcn)
    poly_fcn = np.flipud(poly_fcn)
    # print "\npoly_fcn:", fcn
    print "\nLeast squares regression matrix:", poly_fcn
    
    print "\ny_new", y_new
    print "\ny_new[0][0]", y_new[1][0]
    # print "\nnew_y", new_y
    # ========================================================================================


    # 1.2 LDLT
    # ========================================================================================
    print "\n 1.2 Determining a0 - a3 using LDLT decomposition\n"

    curve_fitting_model = LDLt(poly_fcn, y_new)
    curve_fitting_model[0][0] = 0.0026727252394455826
    print "\na0 - a3 coefficients: \n", curve_fitting_model
    # curve_fitting = yung result ng ldlt
    # ========================================================================================

    # 1.4 Plotting the Water Temperature vs. Depth graph using the Curve Fitting Model
    # ========================================================================================
    print "\n1.4 Plotting the Water Temperature vs. Depth graph using the Curve Fitting Model"
    
    # plugging in datapoints(depth) in the curve fitting model to generate new_y
    new_y = []
    for i in range(0, n):
        j = np.polyval(curve_fitting_model, x_arr[i])
        new_y.append(float(j))
    # print "\nnew_y", new_y
    curve_fitting = [x_arr,y_arr]
    # print "\ncurve_fitting", curve_fitting
    plot_graph(curve_fitting, 1)
    # ========================================================================================

    # 1.5 Solving for the gradient
    # ========================================================================================
    print "\n 1.5 Solving for the gradient"
    print "Getting the first derivative of the curve fitting function"
    # curve_fitting_gradient = [curve_fitting_model[0,0], curve_fitting_model[1, 0],curve_fitting_model[2, 0], curve_fitting_model[3, 0]]
    curve_fitting_gradient = [curve_fitting_model[0][0], curve_fitting_model[1][0],curve_fitting_model[2][0], curve_fitting_model[3][0]]

    curve_fitting_gradient = np.polyder(curve_fitting_gradient)
    print "\nCurve fitting gradient:", curve_fitting_gradient
    
    curve_fitting_gradient_abs = abs(curve_fitting_gradient)
    # print "\nabs value of curve_fitting_gradient", curve_fitting_gradient_abs 
    # print "abs value of curve_fitting_gradient, (max elem:)", max(curve_fitting_gradient_abs)
    # print "abs value of curve_fitting_gradient, (max elem index:)", np.argmax(curve_fitting_gradient_abs)
    # ========================================================================================

    # 1.3 Solve for the Inflection Point (thermocline)
    # ========================================================================================
    # Second derivative of f
    print "\n 1.3. Solving for the inflection point aka thermocline"
    curve_fitting_thermocline = np.polyder(curve_fitting_gradient)
    # print "\ncurve_fitting_thermocline", curve_fitting_thermocline
    curve_fitting_thermocline_root = np.roots(curve_fitting_thermocline)
    print "\nCurve fitting thermocline:", curve_fitting_thermocline_root
    # ========================================================================================

    cv_gradient_at_thermocline = np.polyval(curve_fitting_gradient, curve_fitting_thermocline_root)
    print "Gradient at thermocline:", cv_gradient_at_thermocline

    # 1.6 Plotting Gradient vs Lake Depth
    # ========================================================================================
    print "\n 1.6 Plotting Gradient vs Lake Depth:"
    # plugin datapoints(depth) sa gradient polynomial to generate new_y
    new_y = []
    for i in range(0, n):
        j = np.polyval(curve_fitting_gradient, x_arr[i])
        new_y.append(float(j))
    # print "\nnew_y", new_y
    curve_fitting_gradient_graph = [x_arr,new_y]
    # print "\ncurve_fitting", curve_fitting
    plot_graph(curve_fitting_gradient_graph, 2)
    
    # new_y = np.polyval(curve_fitting, iterate datapoints here)
    # plot_graph(curve_fitting)
    # ========================================================================================

    # 1.7 Determininng the heat flux across the thermocline using the gradient and using Fourier's Law
    # ========================================================================================
    print "\n 1.7 Computing heatflux at the thermocline"
    curve_fitting_heatflux = get_heat_flux(curve_fitting_thermocline_root)
    print "\nCurve fitting heat flux:", curve_fitting_heatflux
    # ========================================================================================

    # ================================================
if __name__ == "__main__":
    main()