# Misa, Bernadette B.
# CS 131 MP1

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import linalg, optimize
import pprint

# Ly=b
# Solves for y in Ly = b, where L is a lower triangular matrix
def forward_sub(L,b):
    y = []
    for i in range(len(b)):
        y.append(b[i])
        for j in range(i):
            y[i] = y[i] - (L[i,j] * y[j])
        y[i] = y[i] / L[i,i]

        # print "y", y
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
        
        print "i, U[i,i]", i, U[-i][-i]

        x[-i] = x[-i] / U[-i][-i]

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
        plt.title('Monomial Basis: Water Temperature (C) vs Depth (m)')
        plt.ylabel("Water Temperature (Celsius)")
    if(item_no == 2):
        plt.title('Monomial Basis: Gradient (dT/dz) vs Depth (m)')
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
    # fillMatrix(L,len(A))
    # fillMatrix(U,len(A))

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
    specific_heat = 1.0

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

    # ========================================================================================

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
    plot_graph(monomial_temp_v_depth, 1)
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
    plot_graph(monomial_gradient_graph, 2)
    

    # ========================================================================================

if __name__ == "__main__":
    main()