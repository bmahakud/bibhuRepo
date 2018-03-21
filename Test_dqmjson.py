from x509auth import * 
from ROOT import TBufferFile, TH1F, TProfile, TH1F, TH2F
import re
from math import *
from dqmjson import *
from ROOT import TFile, gStyle, TCanvas, TH1F, TLegend
from optparse import OptionParser
from xml.dom.minidom import parseString
from rrapi import RRApi, RRApiError
import xmlrpclib
#import elementtree.ElementTree as ET
import sys, os, os.path, time, re, subprocess
import urllib
import json



X509CertAuth.ssl_key_file, X509CertAuth.ssl_cert_file = x509_params()



def DQMget(server, run, dataset, folder, rootContent=False):
      postfix = "?rootcontent=1" if rootContent else ""
      datareq = urllib2.Request(('%s/data/json/archive/%s/%s/%s%s') % (server, run, dataset, folder, postfix))
      datareq.add_header('User-agent', ident)
      data = eval(re.sub(r"\bnan\b", "0", urllib2.build_opener(X509CertOpen()).open(datareq).read()),
	                     { "__builtins__": None }, {})
      if rootContent:
          print "printing len(data): ",len(data['contents'])
	  print "dict: ", data['contents']
          for idx,item in enumerate(data['contents']):
               if 'obj' in item.keys():
		   print item.keys()
                   if 'rootobj' in item.keys():
		       a = array('B')
                       a.fromstring(item['rootobj'].decode('hex'))
		       print "lena:  ",len(a)
	               t = TBufferFile(TBufferFile.kRead, len(a), a, False)
		       rootType = item['properties']['type']
		       print rootType
		       if rootType=='TH1F': rootType ='TH1F'
                       if rootType == 'TPROF': rootType = 'TProfile'
                       if rootType == 'TPROF2D': rootType = 'TProfile'
                       data['contents'][idx]['rootobj'] = t.ReadObject(eval(rootType+'.Class()'))

      return dict( [ (x['obj'], x) for x in data['contents'][1:] if 'obj' in x] )




lastweek=[311532,  311533,  311534,  311535,  311536,  311537,311651,  311614,  311649,  311650,  311662,  311663,  311665,  311668,  311669,  311677,  311678,  311683,  311688,311696,  311702,311727,  311728,  311730,  311733,  311785,  311787,  311789,  311790,  311791,  311792 , 311793,  311795,  311796,  311797 , 311798,  311799,  311800,  311802,  311808,311811,  311814,  311822,311879,  311880,  311881 , 311882,  311883,  311884,311974,  311983,  311989,  311993,  311994,  311995,  312001,  312002,  312003,  312005,  312018,  312023,  312073,  312081, 311949,  311957,  312050,  312059,  312060,  312061,  312091,  312092,  312094,  312096,  312097,  312098,  312099,  312100,  312101]

thisweek=[ 312002, 312059, 312060, 312061, 312073, 312081, 312091, 312092, 312094, 312001, 321098, 312099, 312100, 312142, 312144, 312219, 312135, 312136, 312174, 312176, 312202, 312217, 312218,312141, 312105]

runlist=[]

runlist=lastweek
totalAlcaTracks=0
for run in range(0,len(runlist)):
   dataset="/StreamExpressCosmics/Commissioning2018-Express-v1/DQMIO"
   tmp=DQMget(serverurl, runlist[run], dataset, "PixelPhase1/Tracks/")
   ntrackpixel = tmp['ntracks']['nentries']
   print "alltracks",ntrackpixel 
print "totalNum of alca tracks:  ",totalAlcaTracks    



















