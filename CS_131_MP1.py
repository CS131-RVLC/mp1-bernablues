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
    print "\ny_vector_final", y_vector_final

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

def plot_graph(graph, number, category):

    item_no = number
    graph_type = category
    fcn = graph

    plt.plot(fcn[0][:], fcn[1][:])
    
    if(item_no == 1):
        label = "Curve Fitting"
    if(item_no == 2):
        label = "Monomial"        
    if(item_no == 3):
        label = "Cubic Splines"


    if(graph_type == 1):
        plt.title(label + 'Water Temperature (C) vs Depth (m)')
        plt.ylabel("Water Temperature (Celsius)")
    if(graph_type == 2):
        plt.title(label + 'Gradient (dT/dz) vs Depth (m)')
        plt.ylabel("Gradient (dT/dz)")
        
    plt.xlabel("Depth (meter)")
    plt.show()


def generate_Vandermonde_matrix(data):

    x_arr = data

    n = len(x_arr)

    V = generate_identity_matrix(n)
    # print "\nV before:", V
    
    for i in range(0, n):
        for j in range(0, n):
            V[i,j]=np.power(x_arr[i],j)
    # print "\nV after:", V
    # print "\nV after:"

    # pprint.pprint(V)
    return V

# Solves Ax=b using LU decomposition
def LU_decomp(A_matrix, b):
    A = A_matrix
    n = len(A)
    L = generate_identity_matrix(n)
    U = generate_identity_matrix(n)


    for i in range(0,n):
        L[i,0] = A[i,0]
        U[i,i] = float(1)
        for j in range(0,i):
            U[i,j] = float(0)
        if i==0:
            for j in range(1,n):
                U[0,j] = A[0,j] / L[0,0]
    
    for i in range(1,n):
        for j in range(i,n):
            L[j,i] = A[j,i]
            U[i,j] = A[i,j]
            # print "i, j", i, j
            for k in range(0,i):
                L[j,i] -= (L[j,k] * U[k,i])
            for k in range(0,i):
                U[i,j] -= (L[i,k] * U[k,j])
            U[i,j] = U[i,j] / L[i,i]

    print "\nL", L
    print "\nU", U

    print "\nb", b
    # Ly=b
    y = forward_sub(L,b)
    print "\ny", y
    # Ux=y
    # x = backward_sub(U,y)
    x = np.linalg.solve(U,y)

    print "\nx", x

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

def poly_der(x):
    output_arr = []
    for i in range(1,len(x)):
        temp = float(i) * x[i]
        output_arr.append(temp)
    return output_arr

def poly_val(f,x):
    result = f[0]
    for i in range(1,len(f)):
        temp = f[i] * ((x)**(i))
        result += temp
    return result

def get_inflection_pt(poly_fcn):
    poly_fcn_1 = poly_der(poly_fcn)
    poly_fcn_2 = poly_der(poly_fcn_1)

    x = -poly_fcn_2[1] / poly_fcn_2[0]
    y = 0
    for i in range(0,len(poly_fcn)):
        y += poly_fcn[i] * (x**(float(len(poly_fcn)-i-1)))

    return (x,y)


def cubic_splines_wiki(C):
    # Inspired by the pseudocode in https://en.wikipedia.org/wiki/Spline_(mathematics)
    # C is set of coordinates

    # x_arr = depth
    x_arr = C[0,:]
    # print "x_arr", x_arr
    # y_arr = temp
    y_arr = C[1,:]
    # print "y_arr", y_arr

    n = len(x_arr)-1
    # print "\n", n
    a=[]
    b=[]
    c=[]
    d=[]
    h=[]
    alpha=[]
    l=[]
    mu=[]
    z=[]

    print "\nh", h

    # Step 1
    for i in range(0,n):
        h.append(x_arr[i+1]-x_arr[i])
        alpha.append(0.0)
        a.append(y_arr[i])
        b.append(0.0)
        c.append(0.0)
        d.append(0.0)
    a.append(y_arr[n])
    b.append(0.0)
    d.append(0.0)

    # Step 2
    alpha[0] = (3.0*(a[1]-a[0])/h[0])
    alpha.append(-3.0*(a[n]-a[n-1])/h[n-1])

    # Step 3
    for i in range(1, n):
        alpha[i] = (3.0/h[i])*(a[i+1]-a[i]) - (3.0/h[i-1])*(a[i]-a[i-1])

    print "\nalpha", alpha

    # Step 4
    # Setting up initial variables
    l.append(2.0*h[0])
    mu.append(0.5)
    z.append(alpha[0]/l[0])

    print "\nc", c
    print "\nl", l
    print "\nmu", mu
    print "\nz", z

    # Step 5
    for i in range(1, n):
        l.append((2.0*(x_arr[i+1] - x_arr[i-1])) - (h[i-1]*mu[i-1]))
        mu.append(h[i]/l[i])
        z.append(((alpha[i]-(h[i-1]*z[i-1]))/l[i]))

    print "\nl", l
    print "\nmu", mu
    print "\nz", z

    # Step 6
    l.append( h[n-1]*(2.0-mu[n-1]) )
    z.append( (alpha[n]-(h[n-1]*z[n-1])) / l[n] )
    c.append(z[n])
    
    print "\nl", l
    print "\nz", z
    print "\nc", c

    # Step 7
    for j in range(n-1, 0, -1):
        c[j] = z[j] - (mu[j]*c[j+1])
        b[j] = ((a[j+1]-a[j])/h[j]) - ((h[j]*(c[j+1]+2.0*c[j]))/3.0)
        d[j] = (c[j+1]-c[j])/(3.0*h[j])

    print "\nfinal"
    print "\nh", h
    print "\na", a
    print "\nb", b
    print "\nc", c
    print "\nd", d
    print "\nalpha", alpha
    print "\nl", l
    print "\nmu", mu
    print "\nz", z

    output_set = []


    for i in range(0, n):
        S = [a[i], b[i], c[i], d[i]]
        # S = [a[i], b[i], c[i], d[i], x_arr[i]]
        # S = [d[i], c[i], b[i], a[i]]
        output_set.append(S)

    # print "\noutput_set", output_set
    print "\nSplines coefficients:"
    pprint.pprint(output_set)

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
    plot_graph(curve_fitting, 1, 1)
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
    plot_graph(curve_fitting_gradient_graph, 1, 2)
    
    # new_y = np.polyval(curve_fitting, iterate datapoints here)
    # plot_graph(curve_fitting)
    # ========================================================================================


    # 1.7 Determininng the heat flux across the thermocline using the gradient and using Fourier's Law
    # ========================================================================================
    print "\n 1.7 Computing heatflux at the thermocline"
    curve_fitting_heatflux = get_heat_flux(curve_fitting_thermocline_root)
    print "\nCurve fitting heat flux:", curve_fitting_heatflux
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
    # V = np.fliplr(V)
    # V = np.flipud(V)
    print "\nV:", V
    # ========================================================================================

    # 2.2 Determining a0-a7 using LU Decomposition
    # ========================================================================================
    print "\n 2.2 Determining a0 - a7 using LU decomposition\n"

    monomial_poly_fcn = LU_decomp(V, y_arr)
    # LU(V)
    print "\na0 - a7 coefficients:", monomial_poly_fcn
    
    # print "\ny", y
    # print "monomial_poly_fcn:", monomial_poly_fcn

    # ========================================================================================

  # 2.5 Solving for the gradient
    # ========================================================================================
    print "\n 2.5 Solving for the gradient"
    print "Getting the first derivative of the monomial function"

    monomial_poly_gradient = [monomial_poly_fcn[0], monomial_poly_fcn[1],monomial_poly_fcn[2], monomial_poly_fcn[3], monomial_poly_fcn[4], monomial_poly_fcn[5], monomial_poly_fcn[6], monomial_poly_fcn[7]]
    monomial_poly_gradient = poly_der(monomial_poly_gradient)
    print "\nmonomial gradient:", monomial_poly_gradient

    # # ========================================================================================

    # 2.3 Solve for the inflection point
    # # ========================================================================================
    print "\n 2.3 Solving for the inflection point"

    monomial_thermocline = get_inflection_pt(monomial_poly_fcn)
    print "\nmonomial thermocline", monomial_thermocline

    # 2.4 Plotting the Water Temperature vs. Depth graph using Monomial Basis Interpolation
    # ========================================================================================
    print "\n 2.4 Plotting the Water Temperature vs Depth graph using the Monomial Model"
    # ========================================================================================
    new_x = list(np.linspace(0, 27.2, num=100))
    new_y = []
    for i in range(0, len(new_x)):
        j = poly_val(monomial_poly_fcn, new_x[i])
        new_y.append(float(j))
    monomial_temp_v_depth = [new_x,new_y]
    plot_graph(monomial_temp_v_depth, 2, 1)
    # ============================


    # 2.6 Determining the heat flux across the thermocline using the gradient and using Fourier's Law
    # ========================================================================================
    print "\n 2.6 Computing for the heatflux"
    monomial_heatflux = get_heat_flux(monomial_thermocline[0])
    monomial_real_thermocline = poly_val(monomial_poly_gradient, monomial_thermocline[0])
    print "\nmonomial_heatflux", monomial_real_thermocline
    # # ========================================================================================

    # 2.7 Plotting the Gradient vs Lake Depth
    # ========================================================================================
    print "\n 2.7 Plotting Gradient vs Lake Depth:"
    # plugin datapoints(depth) sa gradient polynomial to generate new_y
    new_y = []
    # for i in range(n-1, -1, -1):
    for i in range(len(new_x)-1, -1, -1):
        # j = np.polyval(monomial_poly_gradient, x_arr[i])
        # j = poly_val(monomial_poly_gradient, x_arr[i])
        j = poly_val(monomial_poly_gradient, new_x[i])
        new_y.append(float(j))
    # monomial_gradient_graph = [x_arr, new_y]
    monomial_gradient_graph = [new_x, new_y]
    plot_graph(monomial_gradient_graph, 2, 2)
    
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