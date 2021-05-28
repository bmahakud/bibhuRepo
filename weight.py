import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import ROOT




def fillhist(h,x1,weight):
    itr=0
    for var1 in x1:
        h.Fill(var1,weight[0][itr])
        #print (var1,"--",weight[0][itr])
        itr=itr+1
    h.Draw()


# an nuphy array with two variables and 12 events
x = np.array([[1, 2, 2,2,3,3,4,5,5,1,1,1],[1, 4, 4,2,3,5,4,5,5,1,1,1]])

eventWeight= np.array([[0.1,0.2,1,1,1,1,0.5,0.1,0.3,0.6,0.4,0.7]])



print ("shape: ",x.shape)
print ("weight: ", eventWeight.shape)


h = ROOT.TH1F("h", "h", 5, 1, 6)

fillhist(h,x[0],eventWeight)

































