import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy
from PIL import Image
from scipy import ndimage


nCol = 6


#trainingDataFile = "TrainDataStLineDecimal.txt"
trainingDataFile = "../input/training-data/TrainData.txt"
#trainingDataFile = "TrainDataStLine.txt" 
#testDataFile = "TestDataStLine.txt"
#trainingDataFile = "TrainData.txt"
#testDataFile = "TestData.txt"





def LoadDataFromTxtFile(InputFile, nCol):
    X = np.zeros((nCol,1))
    Y = np.zeros((1,1))
    with open(InputFile) as file:
       for line in file:
          linedata=line.rstrip()
          linedataArr=linedata.split(" ");
          npArr = np.array([[1]])
          performance = None
          itr=0;
          for var in linedataArr:
              itr=itr+1
              if itr != nCol:
                  element = np.array([[var]])
                  npArr=np.append(npArr, element, axis=1)
              if itr == nCol:
                  element = np.array([[var]])
                  Y=np.append(Y, element, axis=1);
          X=np.append(X, npArr.T, axis=1)

    X = X[:,1:]
    Y = Y[:,1:]
    trainX = np.array(X,dtype=float)
    trainY = np.array(Y,dtype=float)

    return trainX,trainY;





Xtrain,Ytrain = LoadDataFromTxtFile(trainingDataFile, nCol)
#Xtest,Ytest = LoadDataFromTxtFile(testDataFile, nCol)

print ("Xtrain.shape: ", Xtrain.shape)
print ("Xtrain: ", Xtrain)
print ("Ytrain.shape: ", Ytrain.shape)
print ("Ytrain: ", Ytrain)

#print ("Xtest.shape: ", Xtest.shape)
#print ("Ytest.shape: ", Ytest.shape)


mTrain = Xtrain.shape[1]

print ("Number of training examples: ",mTrain)
Theta = np.random.rand(nCol,1);

#Theta=np.array([[0.53635237],[0.34235451]])

Theta=np.array(Theta,dtype=float);

print ("Initial value of Theta: ", Theta)



def AvgCost(Theta, Xtrain, Ytrain):
    #print ("------------Theta.shape: ", Theta.shape)
    mTrain = Xtrain.shape[1];
    mTrainFY = Ytrain.shape[1];
    if (mTrain != mTrainFY):
        raise ValueError('There different number of training examples in X and Y')
    hx = np.dot(Theta.T,Xtrain);
    dy = (hx - Ytrain);
    cost = (1./(2*mTrain))*dy*dy;
    SumCost = np.sum(cost);
    return SumCost



def train(Xtrain, Ytrain, Theta, n_iterations,alpha):
    mTrain = Xtrain.shape[1];
    mTrainFY = Ytrain.shape[1];
    if (mTrain != mTrainFY):
        raise ValueError('There different number of training examples in X and Y')

   
    for iteration in range(1,n_iterations):
        hx = np.dot(Theta.T,Xtrain);
        dy = (hx - Ytrain)
        cost = AvgCost(Theta, Xtrain, Ytrain)
        if iteration % 10 == 0:
            binNum=round(iteration/10);
            print ("iteration: ", iteration, " : cost",cost)
            
        Theta = Theta - (alpha/mTrain) * (np.dot(Xtrain , dy.T))
    return Theta, histCost



Theta, histCost = train(Xtrain, Ytrain, Theta, n_iterations=100000, alpha=0.05);


print ("Fitted Theta: ", Theta)

#costTest = AvgCost(Theta, Xtest, Ytest);

#print ("Cost Test set: ", costTest)

def getCostVsThetaHist(nbin_theta0,theta0_start,theta0_end,nbin_theta1,theta1_start,theta1_end):
    #h_cost_vs_theta = ROOT.TH2F("h_cost_vs_theta","h_cost_vs_theta",nbin_theta0,-5,5,nbin_theta1,-5,5);
    dtheta0 = (theta0_end - theta0_start)/nbin_theta0;
    dtheta1 = (theta1_end - theta1_start )/nbin_theta1;
    nBins = nbin_theta0*nbin_theta1;
    Tg = ROOT.TGraph2D(nBins);
    Tg.SetTitle("Cost vs (#theta_{0}, #theta_{1}); #theta_{0}; #theta_{1}; cost");
    binNo=0
    for ibin0 in range(1,nbin_theta0):
        for ibin1 in range(1, nbin_theta1):
            binNo=binNo+1
            Theta0 = theta0_start + ibin0*(dtheta0)
            Theta1 = theta1_start + ibin1*(dtheta1)
            Theta = np.array([[Theta0],[ Theta1]]) 
            print ("Theta.shape:-----Jai ho ", Theta.shape)
            cost =  AvgCost(Theta, Xtrain, Ytrain)
            Tg.SetPoint(binNo, Theta0, Theta1, cost);
    return Tg



#tGraph = getCostVsThetaHist(nbin_theta0=500,theta0_start=-5,theta0_end=5,nbin_theta1=500,theta1_start=-5,theta1_end=5)


#tGraph.Draw("surf1");




y_pred = np.dot(Theta.T,Xtrain);

XtrainNew = np.delete(Xtrain, 0, 0);

XtrainAttendance=np.array([Xtrain[1]]);
XtrainAlertness=np.array([Xtrain[2]]);
XtrainHomework=np.array([Xtrain[3]]);
XtrainUnderstanding=np.array([Xtrain[4]]);
XtrainPrevPerformance=np.array([Xtrain[5]]);

print (np.array([Xtrain[1]]))
#print ("Xtrain.shape: ", Xtrain.shape);
#print ("XtrainAttendance: ", XtrainAttendance);
#print ("XtrainAttendance.shape: ", XtrainAttendance.T.shape);


#print ("y_pred.T.shape: ", y_pred.T.shape);


#plt.scatter( XtrainAttendance.T, Ytrain.T, color='green');
#plt.xlabel("Attendance")
#plt.ylabel("Performance")
#plt.show();
#plt.savefig('AttendanceVsPerformanceProjection1D.png')




#plt.scatter( XtrainHomework.T, Ytrain.T, color='green');
#plt.xlabel("Homework")
#plt.ylabel("Performance")
#plt.show();
#plt.savefig('HomeworkVsPerformanceProjection1D.png')




#plt.plot(XtrainAttendance.T, y_pred.T, color='red')

#plt.show();
#plt.savefig('StraightLinFitScratchUncertAdded.png')






