dictr={"Short":'T',"Buy":'B',"Sell":'S'}



def solution(hrt_trades, market_trades):
    missing_trades = []
    HRTArr=[]
    MRTArr=[]
    for a in hrt_trades:
        hrtlist= a.split(",")
        HRTArr.append(hrtlist)
    for b in market_trades:
        mrttrd = b.split(",")
        MRTArr.append(mrttrd)
    #print (len(HRTArr))
    #print (len(MRTArr))

    for itr1 in HRTArr:
        #print ( "-------------searching for id in MRTTr",itr1[-1])
        idtotest=itr1[-1]
        for itr2 in MRTArr:
             if itr2[1]==idtotest:
                 #print("found id")
                 #check these following
                 appnd="no"
                 print (dictr.get(str(itr1[3])))
                 if itr1[0] !=itr2[4]:
                     appnd="yes"
                 if  itr1[1] !=itr2[5]:
                     appnd="yes"
                 if itr1[2] != itr2[3]:
                     appnd="yes"
                 if  dictr.get(str(itr1[3])) != itr2[2]:
                     appnd="yes"
                 if itr1[5] != itr1[0]: 
                     appnd="yes"

                 missing_trades.append(itr2)




    return missing_trades

print (solution(hrt_trades, market_trades))

print (len(solution(hrt_trades, market_trades)))
