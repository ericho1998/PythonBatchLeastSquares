import numpy as np
import matplotlib.pyplot as plt
import math as m

def batch_LSQ():
    #Batch Least Squares
    #Load in data
    #Data is an epoch, and a range to [0,0], [100,0], [100,100], and [0,100] in that order
    data = np.loadtxt(open("Lab2data.txt", "r"), delimiter = " ")

    #Define the locations of the 4 targets which provide the ranges
    t = np.array([[0,0],[100,0],[100,100],[100,0]])

    #Create an array to store all adjusted predicted x values
    x_final = np.zeros((np.size(data,0), 2))

    #Create a residual vector for all ranges in each epoch
    v = np.zeros((np.size(data,0), 4))

    for i in range(0,np.size(data,0)):
        #Set the measurement for the epoch as the 4 ranges
        l = data[i,1:]
        #Initial coordinates
        x = np.array([50,50])
        x = x.astype(float)

        #set a convergence criterion
        delta = m.inf

        #Initialize Design Matrix (A) and Misclosure Vector (w)
        A = np.zeros((np.size(t,0), np.size(x,0)))
        w = np.zeros((np.size(t,0), 1))
        while np.sum(np.abs(delta)) > 0.0001:

            #Calculate A and w
            for j in range(0,np.size(t,0)):
                dist = m.sqrt(pow((t[j,0] - x[0]),2) + pow((t[j,1] - x[1]),2))
                A[j,0] = (t[j,0] - x[0])/dist
                A[j,1] = (t[j,1] - x[1])/dist
                w[j] = l[j] - dist

            #Update delta
            delta = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.transpose(A),A)),np.transpose(A)),w)

            #Update x for iteration
            x[0] = x[0] - delta[0,0]
            x[1] = x[1] - delta[1,0]
        #When convergence is met, add x to the list of finalized x coordinates
        x_final[i,:] = x

        #Calculate residuals
        for j in range(0, np.size(t,0)):
            dist = m.sqrt(pow((t[j,0] - x[0]),2) + pow((t[j,1] - x[1]),2))
            v[i,j] = l[j] - dist

        #Calculate aposteriori variance factor
        apost = np.matmul(np.transpose(v[i]), v[i]) / (np.size(l,0) - np.size(x))

        if i == 0:
            A1 = A
            apost1 = apost
            print(A1)
        cxhat = apost1 * np.linalg.inv(np.matmul(np.transpose(A1), A1))

        #Error ellipse stuff
        a2 = 0.5 * ((pow(cxhat[0,0],2) + pow(cxhat[1,1],2) + m.sqrt(pow((cxhat[0,0] - cxhat[1,1]),2) + (4 * cxhat[0,1]))))
        b2 = 0.5 * ((pow(cxhat[0,0],2) + pow(cxhat[1,1],2) - m.sqrt(pow((cxhat[0,0] - cxhat[1,1]),2) + (4 * cxhat[0,1]))))
        theta = 90 - (0.5 * m.atan2((2*cxhat[0,1]), (pow(cxhat[0,0],2) - pow(cxhat[1,1],2)))) * 180 / m.pi




    #Plot the positions obtained
    fig, ax = plt.subplots()
    ax.plot(x_final[:,0], x_final[:,1])
    ax.set_xlabel('X Coordinate (m)')
    ax.set_ylabel('Y Coordinate (m)')
    ax.set_title('X and Y Coordinates with Bundle Least Squares Adjustment')
    plt.show()





def main():
    batch_LSQ()

if __name__ == '__main__':
    main()
print('complete')
