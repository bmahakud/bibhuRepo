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
                 if itr1[0] !=itr2[4] or itr1[1] !=itr2[5] or itr1[2] != itr2[3] or itr1[3][0] != itr2[2]:

                    #print("found problem")
                    missing_trades.append(itr1)
                    missing_trades.append(itr2)




    return missing_trades

print (solution(hrt_trades, market_trades))
