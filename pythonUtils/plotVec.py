import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys
import csv
import os


delimeter = ' ';
dataset = []
samples = []
originalLabels = []
mapOld2NewLabels = dict()
newLabels = []
colorTable = []
nSamples = 0
nFeatures = 0
nClasses = 0;


def readDataset():
    messageInputDataSet = "\n" \
                          "Digite o filepath do dataset\n"
    try:
        filepath = raw_input(messageInputDataSet)
        print filepath
        if(len(filepath) == 0 ):
            return
        global dataset,samples,originalLabels,colorTable,newLabels
        global nSamples,nFeatures,nClasses
        global mapOld2NewLabels
        dataset = np.loadtxt(filepath, delimiter=delimeter)
        print dataset
        [nSamples, nFeatures] = dataset.shape
        samples = dataset[:,:]
        newLabels = range(0,nSamples);

        print "Dataset carregado com sucesso\n"
    except ValueError:
        print "Houve um erro ao tentar carregar o dataset: " + str(ValueError)

    return

def projection():
    if (len(dataset) == 0):
        print "nenhum dataset carregado\n"
        return


    colors = ['b', 'c', 'y', 'm', 'r','k','g']
    plots = list();
    legends = list();
    print samples
    for i in range(0,nSamples):
        plt.plot(range(0,nFeatures), samples[i,0:nFeatures], color=colors[i],label=str(i));

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()



def clear():
    if (os.name == 'nt'):
        os.system('cls')
    else:
        os.system('clear')

def changeDelimeter():
    global delimeter;
    messageInputDelimeter = "\n" \
                          "Digite o novo delimitador:\n"
    temp = raw_input(messageInputDelimeter)
    if (len(temp) == 0):
        return
    delimeter = temp
    print "novo delimitador definido com sucesso\n"


inFile = sys.argv[1]
dataset = np.loadtxt(inFile, delimiter=' ')
[nSamples, nFeatures] = dataset.shape
colors = ['b', 'c', 'y', 'm', 'r','k','g']
for i in range(0, nSamples):
    plt.plot(range(0, nFeatures), dataset[i, 0:nFeatures], color=colors[i], label=str(i));

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()

