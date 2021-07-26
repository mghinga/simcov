# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main(dimension):
    airways = pd.read_csv('airway.csv', header=None)
    airways = airways.values
    airways = airways[~(np.isnan(airways))]
    X, Y, Z = [],[],[]
    for cell in airways:
        Z.append(cell / (dimension * dimension))
        cell = cell % (dimension * dimension)
        Y.append(cell / dimension)
        X.append(cell % dimension)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(X, Y, Z,  color='red', marker='o')
    ax.set_xlim([0, dimension])
    ax.set_ylim([0, dimension])
    ax.set_zlim([0, dimension])
    ax.set_xlabel('X', fontsize=24)
    ax.set_ylabel('Y', fontsize=24)
    ax.set_zlabel('Z', fontsize=24)
    plt.show()
    #plt.savefig('lung.png', dpi=300)
    #plt.clf()


if __name__ == "__main__":
    main(dimension=300)
