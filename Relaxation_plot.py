import matplotlib.pyplot as plt
import numpy as np


if __name__ == "__main__":
    fname_Epidmicdata = 'problem6_1_data.txt'
    data = np.loadtxt(fname_Epidmicdata,delimiter = ',', skiprows = 1)

    r0_values = data[:,0]
    solutions = data[:,1]

    plt.plot(r0_values,solutions,color = 'green')
    plt.ylabel('Probability of Epidemic (P)')
    plt.xlabel('Those infected by sick individual (R0)')
    plt.title('Epidemic probability as r0 increases')
    plt.show()