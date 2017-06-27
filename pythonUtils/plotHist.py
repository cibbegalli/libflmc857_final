import csv
import numpy as np
import matplotlib.pyplot as plt
import sys

#1 - input File
#2 - row or col based
#3 - split figures (color image only)

def plotGrayHistogram(x,y):
  plt.bar(x,y)
  plt.title("Histogram")
  plt.xlabel("bins")
  plt.ylabel("counts")
  plt.show();

def plotColorHistogram(x,r,g,b,split):
  if(split):
    plt.figure(0)
    plt.bar(x,r)
    plt.title("Histogram Red Channel")
    plt.xlabel("bins")
    plt.ylabel("counts")

    plt.figure(1)
    plt.bar(x,g)
    plt.title("Histogram Green Channel")
    plt.xlabel("bins")
    plt.ylabel("counts")

    plt.figure(2)
    plt.bar(x,b)
    plt.title("Histogram Blue Channel")
    plt.xlabel("bins")
    plt.ylabel("counts")

    plt.show();

  else:
    plt.figure(0)
    plt.subplot(221)
    plt.bar(x,r)
    plt.title("Histogram Red Channel")
    plt.xlabel("bins")
    plt.ylabel("counts")

    plt.subplot(222)
    plt.bar(x,g)
    plt.title("Histogram Green Channel")
    plt.xlabel("bins")
    plt.ylabel("counts")

    plt.subplot(223)
    plt.bar(x,b)
    plt.title("Histogram Blue Channel")
    plt.xlabel("bins")
    plt.ylabel("counts")

    plt.show();


inFile = sys.argv[1]
colBased = True;
split = False;
if(len(sys.argv) >= 3):
  if (sys.argv[2] == 'False'):
    colBased = False;
  else:
    colBased = True;

if(len(sys.argv) >= 4):
  if (sys.argv[3] == 'False'):
    split = False;
  else:
    split = True;

inputData = np.loadtxt(inFile,delimiter=' ')
[row,col] = inputData.shape

if(colBased):
  if(col == 2):
    x = inputData[:, 0]
    y = inputData[:, 1]
    plotGrayHistogram(x,y)
  else:
    x = inputData[:, 0]
    r = inputData[:, 1]
    g = inputData[:, 2]
    b = inputData[:, 3]
    plotColorHistogram(x, r, g, b,split)

else:
  if(row == 2):
    x = inputData[0, :]
    y = inputData[1, :]
    plotGrayHistogram(x, y)
  else:
    x = inputData[0, :]
    r = inputData[1, :]
    g = inputData[2, :]
    b = inputData[3, :]
    plotColorHistogram(x,r,g,b,split)


