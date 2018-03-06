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

##Run classification
#groupName = "Collisions17" #TO_BE_CHANGED_EACH_YEAR
groupName = "Commissioning2018" #TO_BE_CHANGED_EACH_YEAR
#groupName = "Cosmics18CRUZET" #TO_BE_CHANGED_EACH_YEAR
##Dataset for GUI query
express = ['/StreamExpress/', '/StreamExpressCosmics/']
expresshi = ['/StreamHIExpress/', '/StreamExpressCosmics/']
prompt  = ['/ZeroBias/',   '/Cosmics/']
prompt1  = ['/ZeroBias1/',   '/Cosmics/']
prompthi  = ['/HIMinimumBias1/',   '/Cosmics/']
express0t = ['/StreamExpress0T/', '/StreamExpressCosmics/']
prompt0t  = ['/ZeroBias_0T/',   '/Cosmics/']
#yearPattern = ".*17" # select anything with a 17 in the name #TO_BE_CHANGED_EACH_YEAR
yearPattern = ".*18" # select anything with a 17 in the name #TO_BE_CHANGED_EACH_YEAR
##Workspace and dset types
Wkspace = ["GLOBAL", "TRACKER"]
####This is under construction...
##Recotype= ["Online", "Prompt"]
##Selection of GUI query array element
Dtype = 0

### List of people who are not shifters, and whose open runs should be considered "TODO"
NonShifters = [ "DQMGUI Trigger" ] 

####Cosmics settings are set after loading config options#####

os.environ['X509_USER_CERT']='/data/users/11.11a/auth/proxy/proxy.cert'
parser = OptionParser()
parser.add_option("-c", "--cosmics", dest="cosmics", action="store_true",  default=False, help="Check cosmic instead of collision")
parser.add_option("-m", "--min", dest="min", type="int", default=0,      help="Minimum run")
parser.add_option("-M", "--max", dest="max", type="int", default=999999, help="Maximum run")
parser.add_option("--min-ls",    dest="minls",  type="int", default="10",   help="Ignore runs with less than X lumis (default 10)")
parser.add_option("-v", "--verbose", dest="verbose", action="store_true",  default=False, help="Print more info")
parser.add_option("-p", "--pretend", dest="pretend", action="store_true",  default=False, help="Use cached RR result")
parser.add_option("-f", "--force", dest="force", action="store_true",  default=False, help="Never cached RR result")
parser.add_option("-n", "--notes", dest="notes", type="string", default="notes.txt", help="Text file with notes")
(options, args) = parser.parse_args()

eras = []
erafile = open("eras.txt","r")
for line in erafile:
    if "from" in line: continue
    cols = line.split(); 
    if len(cols) != 5: continue
    eras.append( ( int(cols[0]), int(cols[1]), cols[2], cols[3], cols[4] ) )

def eraForRun(run):
    for min,max,era,pr,er in eras:
        if run >= min and run <= max: return era
    return "Unknown"
def getPrForRun(run):
    for min,max,era,pr,er in eras:
        if run >= min and run <= max: return pr
    return "Unknown"
def getErForRun(run):
    for min,max,era,pr,er in eras:
        if run >= min and run <= max: return er
    return "Unknown"

def isExpressDoneInGUI(run):
    dataset = "%s%s-%s/DQMIO" % (express[Dtype], eraForRun(run), getErForRun(run))
    try:
        info = dqm_get_json(serverurl, run, dataset, "Info/ProvInfo")
        done = info['runIsComplete']['value']
        return done == '1'
    except:
        return False 
    return False


####### FIXME set dataset 
#dataset =  "/StreamExpressCosmics/Commissioning2018-Express-v1/DQMIO"


def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])

##Cosmic settings...
#if options.cosmics: groupName = "Cosmics18" #TO_BE_CHANGED_EACH_YEAR
if options.cosmics: groupName = "Cosmics18CRUZET" #TO_BE_CHANGED_EACH_YEAR
##NOTE: Currently using prompt stream (not express) for central data certification
if options.cosmics: Dtype = 1

notes = {}
if options.notes:
    try:
        nfile = open(options.notes, "r");
        for l in nfile:
            m = re.match(r"\s*(\d+)\s*:?\s+(.*)", l)
            if m:
                notes[int(m.group(1))] = m.group(2)
    except IOError:
        print "Couldn't read notes file", options.notes, "Will use Tracker Prompt RECO comments instead"

lumiCache = {}; 
lumiCacheName = "lumi-by-run.txt" if not options.cosmics else "tracks-by-run.txt"

try:
    lumiFile = open(lumiCacheName, "r")
    for l in lumiFile:
        m = re.match(r"(\d+)\s+(\d+)\s+([0-9.]+).*", l)
        if m:
            if not options.cosmics:
                lumiCache[int(m.group(1))] = int(m.group(2)), float(m.group(3))
            else:
                cols = l.split()
                lumiCache[int(cols[0])] = [ int(cols[1]), int(cols[2]), cols[3], cols[4], cols[5].replace("_"," ") ] 
    print "LUMICACHE " , lumiCache
except IOError:
   pass 

runlist = {}

URL = 'http://runregistry.web.cern.ch/runregistry/' #If you want to test the logic of your query, you can access the user RR thanks to this link: https://cmswbmoffshift.web.cern.ch/cmswbmoffshift/runregistry_user/index.jsf
api = RRApi(URL, debug = False)

def getRR(whichRR, dataName):
    global groupName, runreg, runlist, options
    sys.stderr.write("Querying %s RunRegistry for %s runs...\n" % (whichRR,dataName));
    mycolumns = ['pix','strip','track','ranges','runNumber','datasetState','lastShifter']
    text = ''
    fname = "RR_%s.%s.%s.xml" % (whichRR,groupName,dataName)
    readFile = os.path.exists(fname) and options.pretend
    if os.path.exists(fname) and (time.time() - os.stat(fname).st_mtime) < 10*60 and not options.force:
        readFile = True
    if readFile:
        if options.verbose: print "  will read from %s (%.0f minutes old)" % (fname, (time.time() - os.stat(fname).st_mtime)/60)
        log = open(fname); 
        text = "\n".join([x for x in log])
    else:
        ##Query RR
        if api.app == "user":
            text = api.data(workspace = whichRR, table = 'datasets', template = 'xml', columns = mycolumns, filter = {'runNumber':'>= %d and <= %d'%(options.min,options.max),'runClassName':"like '%%%s%%'"%groupName,'datasetName':"like '%%%s%%'"%dataName})
        log = open(fname,"w"); 
        log.write(text); log.close()
    ##Get and Loop over xml data
    #print "HUGO IS UNDER TESTING"
    #print api.workspaces()
    #print api.tables('GLOBAL')
    #print api.columns('GLOBAL', 'runsummary')
    #print api.templates('GLOBAL', 'runsummary')
    print api.templates('TRACKER', 'datasets')
    #print api.tags()
    #print api.data(workspace = 'GLOBAL', table = 'runsummary', template = 'xml', columns = ['number'], filter = {'number':'>= %d and <= %d'%(options.min,options.max)})
    #print api.data(workspace = 'GLOBAL', table = 'runsummary', template = 'xml', columns = ['number'], filter = {'number':'>= %d and <= %d'%(options.min,options.max),'runClassName':"like '%%%s%%'"%groupName,'datasets.datasetName':"like '%%%s%%'"%dataName})
    #print api.data(workspace = 'GLOBAL', table = 'runsummary', template = 'xml', columns = ['number','bfield'], filter = {'runNumber':'>= %d and <= %d'%(options.min,options.max),'runClassName':'like %Collisions17%','datasetName':'like %Online%'})
    #print api.data(workspace = 'GLOBAL', table = 'runsummary', template = 'xml', columns = ['number','bfield'], filter = {"runClassName": "like '%%%s%%'"%groupName, "number": ">= %d AND <= %d" %(options.min,options.max), "datasets": {"rowClass": "org.cern.cms.dqm.runregistry.user.model.RunDatasetRowGlobal", "datasetName": "like %Online%"}}, tag= 'LATEST')
    
    dom = ''; domP = None
    domB = '';
    try:
        dom  = parseString(text)
        print "DOM " , dom
    except:
        ##In case of a non-Standard RR output (dom not set)
        print "Could not parse RR output"
    if whichRR == "GLOBAL" and dataName == "Online": 
        text_bfield = api.data(workspace = whichRR, table = 'runsummary', template = 'xml', columns = ['number','bfield'], filter = {"runClassName": "like '%%%s%%'"%groupName, "number": ">= %d AND <= %d" %(options.min,options.max), "datasets": {"rowClass": "org.cern.cms.dqm.runregistry.user.model.RunDatasetRowGlobal", "datasetName": "like %Online%"}}, tag= 'LATEST')
        #the above query is broken unfortunately... it has to be fixed
	log = open("RR_bfield.xml","w");
        log.write(text_bfield); log.close()
        try:
            domB  = parseString(text_bfield)
            print "DOMB " , domB
        except:
        ##In case of a non-Standard RR output (dom not set)
            print "Could not parse RR output"

    if os.path.exists("patches/"+fname):
        try:
            domP = parseString( "\n".join([x for x in open("patches/"+fname)]) )    
            print "Found manual patch of RR ",fname
        except:
            pass
    splitRows = 'RunDatasetRowTracker'
    if whichRR == 'GLOBAL': splitRows = 'RunDatasetRowGlobal'
    ##Protection against null return
    if dom: data = dom.getElementsByTagName(splitRows)
    else: data =[]
    if domP: dataP = domP.getElementsByTagName(splitRows)
    else: dataP =[]
    if domB: dataB = domB.getElementsByTagName('RunSummaryRowGlobal')
    else: dataB =[]
    for i in range(len(data)):
        ##Get run#
        run = int(data[i].getElementsByTagName('runNumber')[0].firstChild.data)
        if run < options.min: continue
        if run > options.max: continue
        mydata = data[i]
        for X in dataP:
            if int(X.getElementsByTagName('runNumber')[0].firstChild.data) == run:
                print "Run ",run, ": found manual patch for ",whichRR,groupName,dataName,
                mydata = X; break
        state = mydata.getElementsByTagName('datasetState')[0].firstChild.data
        shifter = mydata.getElementsByTagName('lastShifter')[0].firstChild.data
        isopen = (state  == "OPEN")
        lumis= 0
        bfield = -1
        for X in dataB:
            if int(X.getElementsByTagName('number')[0].firstChild.data) == run:
                bfield = X.getElementsByTagName('bfield')[0].firstChild.data
                break
        if run not in runlist: runlist[run] = {'ls':lumis}
        ### PIXEL
        goodp = mydata.getElementsByTagName(mycolumns[0])[0].getElementsByTagName('status')[0].firstChild.data == 'GOOD'
        commp = (mydata.getElementsByTagName(mycolumns[0])[0].getElementsByTagName('comment')[0].toxml()).replace('<comment>','').replace('</comment>','').replace('<comment/>','')
        ### STRIP
        goods = mydata.getElementsByTagName(mycolumns[1])[0].getElementsByTagName('status')[0].firstChild.data == 'GOOD'
        comms = (mydata.getElementsByTagName(mycolumns[1])[0].getElementsByTagName('comment')[0].toxml()).replace('<comment>','').replace('</comment>','').replace('<comment/>','')
        ##No tracking flag for 'Global'/'Online', cosmic data good if strips good...
        if options.cosmics:
            goodt = (goods); commt = ""
        else:
            goodt = (goods and goodp); commt = ""
        if whichRR != 'GLOBAL' and dataName != 'Online':
            ### TRACKING
            goodt = mydata.getElementsByTagName(mycolumns[2])[0].getElementsByTagName('status')[0].firstChild.data == 'GOOD'
            commt = (mydata.getElementsByTagName(mycolumns[2])[0].getElementsByTagName('comment')[0].toxml()).replace('<comment>','').replace('</comment>','').replace('<comment/>','')
        if goodt:
            verdict = "GOOD"
            if not goodp: verdict += ", px bad"
            if not goods: verdict += ", st bad"
        else:
            verdict = 'BAD'
            if goodp: verdict += ", px good" 
            if goods: verdict += ", st good" 
        if options.verbose: print "  -",run,lumis,verdict
        ##Compile comments
        comment = ""
        if commt: comment += commt
        if comms: comment += ", strip: "+comms
        if commp: comment += ", pixel: "+commp
        if isopen and shifter in NonShifters: (isopen, verdict,comment) = (True, "TODO","")
        runlist[run]['RR_'+whichRR+"_"+dataName] = [ isopen, verdict, comment ]
        if whichRR == 'GLOBAL' and dataName == 'Online':
            runlist[run]['RR_bfield'] = float(bfield)
            
        #print "runlist " , runlist[run]

#getRR("GLOBAL", "Online")
#getRR("GLOBAL", "Prompt")
getRR("TRACKER", "Express")
#getRR("TRACKER", "Prompt")
##Start running RR queries
#for work in Wkspace:
#    for reco in Recotype:
#        getRR(work,reco)
        #if options.cosmics and work == "GLOBAL" and reco == "Prompt": getRR(work,"Express")
        #else: getRR(work,reco)

print "Querying runs from DQM GUI"
ed = express[Dtype]
pd = prompt[Dtype]
pd1 = prompt1[Dtype]
pd0t = prompt0t[Dtype]
pdhi = prompthi[Dtype]
ed0t = express0t[Dtype]
edhi = expresshi[Dtype]

for n,d in (('Express',ed), ('Prompt',pd)):
    samples = dqm_get_samples(serverurl, d+yearPattern)
    for (r, d2) in samples:
        if r not in runlist: continue
        runlist[r]['GUI_'+n] = True

if Dtype == 0: #collisions-only
    for n,d in (('Express',ed), ('Prompt',pd1)):
        samples = dqm_get_samples(serverurl, d+yearPattern)
        for (r, d2) in samples:
            if r not in runlist: continue
            runlist[r]['GUI_'+n] = True

    for n,d in (('Express',ed0t), ('Prompt',pd0t)):
        samples = dqm_get_samples(serverurl, d+yearPattern)
        for (r, d2) in samples:
            if r not in runlist: continue
            runlist[r]['GUI_'+n] = True

    for n,d in (('Express',edhi), ('Prompt',pdhi)):
        samples = dqm_get_samples(serverurl, d+yearPattern)
        for (r, d2) in samples:
            if r not in runlist: continue
            runlist[r]['GUI_'+n] = True

if not options.cosmics:
    print "Getting luminosities"
    newcache = open("lumi-by-run.txt", "w");
    newcache.write("run\tls\tlumi_pb\n");
    for run in runlist.keys():
        if run not in lumiCache:
            print " - ",run
            lslumi = (-1,0)
            try:
                os.system("./lumiCalc2_wrapper.sh %d" % run)
                out = [ l for l in open("lumi.tmp","r")]
                if (len(out) <= 1): raise ValueError
                ##quick fix for multiple LS intervals...
                out[1] = out[1].replace("], [", "]; [")
                cols = out[1].strip().split(",");
                print cols
                (myrun,myls,delivered,sells,mylumi) = out[1].strip().split(",")
                myrun = myrun.split(":")[0]
                if int(myrun) == run:
                    lslumi = ( int(myls), float(mylumi)/1.0e6 )
                    if options.verbose: print "\t- %6d, %4d, %6.3f" % (run, lslumi[0], lslumi[1])
            except IOError:
                pass
            except ValueError:
                lslumi = (-1,0)

            try:
                dataset = "%s%s-%s/DQMIO" % (express[0], eraForRun(run), getErForRun(run))
                print dataset
                ei = dqm_get_json(serverurl, run, dataset, "Info/EventInfo")
                myls = ei['ProcessedLS']['nentries']
                lslumi = ( int(myls), 0 )
            except:
                pass

            lumiCache[run] = lslumi
        if lumiCache[run][0] != -1:
            newcache.write("%d\t%d\t%.3f\n" % (run, lumiCache[run][0], lumiCache[run][1]))
    newcache.close()
else:
    #print "Getting APV modes"
    #apvModeList = []; minrun = min(runlist.keys())
    #pyScript = os.environ['CMSSW_RELEASE_BASE']+"/src/CondFormats/SiStripObjects/test/SiStripLatencyInspector.py"
    #pyScript = "SiStripLatencyInspector.py"
    #modeDumpPipe = subprocess.Popen(['python', pyScript], bufsize=-1, stdout=subprocess.PIPE).stdout;
    #for line in modeDumpPipe:
    #    m = re.match(r"since = (\d+) , till = (\d+) --> (peak|deco) mode", line)
    #    if m:
    #        first, last, mode = int(m.group(1)), int(m.group(2)), m.group(3).upper() 
    #        if last >= minrun: apvModeList.append( (first, last, mode) )
    #apvModeList.sort()
    print "Getting tracks"
    newcache = open("tracks-by-run.txt", "w");
    newcache.write("run\tls\talcatracks\tmode\tmode_flag\tmode_text\n");
    for run in runlist.keys():
        if run not in lumiCache:
            print " - ",run
            dbmode = '???'
            #for (start,end,mode) in apvModeList:
            #    if run >= start and run <= end: 
            #        dbmode = mode
            #        break

            link = "http://cern.ch/erik.butz/cgi-bin/getReadOutmode.pl?RUN=" + str(run)
            f = urllib.urlopen(link)
            json_data = f.read()            
            dbmodelist = json.loads(json_data)
            dbmode = dbmodelist[0][1]
            lslumi = (-1,0,dbmode,"WAIT","from DB mode (run not in prompt GUI yet)")
            try:
                dataset = "%s%s-%s/DQMIO" % (prompt[1], eraForRun(run), getPrForRun(run))
                #FIXME#########
                dataset =  "/StreamExpressCosmics/Commissioning2018-Express-v1/DQMIO"
                print "DATASET " , dataset
                at = dqm_get_json(serverurl, run, dataset, "AlCaReco/TkAlCosmics0T/GeneralProperties")
                ei = dqm_get_json(serverurl, run, dataset, "Info/EventInfo")
                tib =dqm_get_json(serverurl, run, dataset, "SiStrip/MechanicalView/TIB")
                nlumis  = ei['ProcessedLS']['nentries']
                nalcatracks = at['Chi2Prob_ALCARECOTkAlCosmicsCTF0T']['nentries']
                ston_num = tib['Summary_ClusterStoNCorr_OnTrack__TIB']['nentries']
                ston_avg = tib['Summary_ClusterStoNCorr_OnTrack__TIB']['stats']['x']['mean']
                mode = "???"; mode_flag = 'bah'; mode_text = 'not found'
                if ston_num > 100:
                    if 28 < ston_avg and ston_avg < 35: mode, mode_flag, mode_text = "PEAK", "TODO", "from S/N plot";
                    if 18 < ston_avg and ston_avg < 24: mode, mode_flag, mode_text = "DECO", "TODO", "from S/N plot";
                if mode == dbmode:  mode, mode_flag, mode_text = dbmode, "GOOD","from both DB and S/N"
                elif mode == "???": mode, mode_flag, mode_text = dbmode, "WAIT","from DB only (S/N info is inconclusive)"
                else: mode, mode_flag, mode_text = dbmode+"?", "BAD","DB says %s, but mean S/N = %.1f suggests %s" % (dbmode,ston_avg,mode)
                lslumi = (nlumis, nalcatracks, mode, mode_flag, mode_text)
            except:
                pass
            if lslumi[1] == 0:
                try:
                    print ["line 388"]
                    dataset = "%s%s-%s/DQMIO" % (express[1], eraForRun(run), getErForRun(run))
                   # print ["line 390"]
                    print dataset
                    at = dqm_get_json(serverurl, run, dataset, "AlCaReco/TkAlCosmics0T/GeneralProperties")
                    #at = dqm_get_json(serverurl, run,"/StreamExpressCosmics/Commissioning2018-Express-v1/DQMIO" , "AlCaReco/TkAlCosmics0T/GeneralProperties")
                    #print ["line 392", at]
                    ei = dqm_get_json(serverurl, run, dataset, "Info/EventInfo")
                    #ei = dqm_get_json(serverurl, run, "/StreamExpressCosmics/Commissioning2018-Express-v1/DQMIO", "Info/EventInfo")
                    #print ["line 394", ei]
                    nlumis  = ei['ProcessedLS']['nentries']
                   # print ["line 396"]
                    nalcatracks = at['Chi2Prob_ALCARECOTkAlCosmicsCTF0T']['nentries']
                    print ["line 398"]
                    if nlumis > 0:
                        lslumi = (-nlumis,nalcatracks,dbmode,"WAIT","from DB mode (run not in prompt GUI yet)")
                        print ["line 401"]
                except:
                    pass
            print "LSLUMI cosmics ", lslumi
            lumiCache[run] = lslumi
        if lumiCache[run][0] >= 0:
            newcache.write("%d\t%d\t%d\t%s\t%s\t%s\n" % (run, 
                lumiCache[run][0], lumiCache[run][1], 
                lumiCache[run][2], lumiCache[run][3], lumiCache[run][4].replace(" ","_")))
    newcache.close()

print "Done"

html = """
<html>
<head>
  <title>Certification Status, %s (%s)</title>
  <style type='text/css'>
    body { font-family: "Candara", sans-serif; }
    td.BAD { background-color: rgb(255,100,100); }
    td.bah { background-color: rgb(255,180,80); }
    td.GOOD { background-color: rgb(100,255,100); }
    td.TODO { background-color: yellow; }
    td.WAIT { background-color: rgb(200,200,255); }
    td.Wait { background-color: rgb(200,230,255); }
    td.SKIP { background-color: rgb(200,200,200); }
    td, th { padding: 1px 5px; 
             background-color: rgb(200,200,200); 
             margin: 0 0;  }
    td.num { text-align: right; padding: 1px 10px; }
    table, tr { background-color: black; }
  </style>
</head>
<body>
<h1>Certification Status, %s (%s)</h1>
<table>
""" % (groupName, time.ctime(), groupName, time.ctime())
if not options.cosmics:
    html += "<tr><th>Run</th><th>B-field</th><th>LS</th><th>LUMI</th><th>ONLINE</th><th>EXPRESS</th><th>PROMPT</th><th>CENTRAL</th><th>NOTES</th></tr>"
else:
    html += "<tr><th>Run</th><th>B-field</th><th>LS</th><th>TRACKS<br/>ALCA</th><th>TRACK RATE<br/>ALCA</th><th>APV<br/>MODE</th><th>ONLINE</th><th>EXPRESS</th><th>PROMPT</th><th>CENTRAL</th><th>NOTES</th></tr>"

def v2c(isopen,verdict):
    if isopen: return 'TODO'
    for X,Y in [('BAD','BAD'), ('bad','bad'), ('GOOD','GOOD'), ('TODO','TODO'), ('WAIT','WAIT'), ('Wait','Wait'),('SKIP','SKIP'),('N/A','SKIP'),('STANDBY','STANDBY'),('EXCLUDED','EXCL')]:
        if X in verdict: return Y
def p2t(pair):
    (isopen, verdict, comment) = pair
    if comment:
        return "%s <span title=\"%s\">[...]</span>" % (verdict, comment)
    else:
        return verdict

allLumi_currentCRUZET=0
allLumi_currentCRAFT=0.000001
allAlcaTracks_currentCRUZET=0
allAlcaTracks_currentCRAFT=0
allLumiWait=0
allTracksWait=0
maxcosmicrunforstat = 0
allAlcaTracksPEAK=0

info_run_track_CRUZET = []; 
info_run_time_CRUZET = []; 
info_run_track_CRAFT = []; 
info_run_time_CRAFT = []; 

#allLumiB=0
#allAlcaTracksB=0

runs = runlist.keys(); runs.sort(); runs.reverse()
print "ALL RUNS: " , runs , "\n"
if options.cosmics: print "ONLINE: Global Online, EXPRESS: Trk Online, PROMPT: Trk Prompt, CENTRAL: Global ExpressStream"
else: print "ONLINE: Global Online, EXPRESS: Trk Online, PROMPT: Trk Prompt, CENTRAL: Global Prompt"
print ""
print "%-6s |  %-15s | %-15s | %-15s | %-15s | %s " % ("RUN","ONLINE","EXPRESS","PROMPT","CENTRAL","NOTES")
print "%-6s |  %-15s | %-15s | %-15s | %-15s | %s " % ("-"*6, "-"*15, "-"*15, "-"*15, "-"*15, "-"*30)
#print runlist
for r in runs:
    #if options.cosmics and lumiCache[r][3] == 'WAIT': continue #ignore cosmic runs in the waiting list (express stream)
    if options.cosmics and lumiCache[r][0] == -1: continue     #ignore irrelevant runs (?)
    R = runlist[r]
    print [' R ', R]
    All_comments=''
    online = R['RR_GLOBAL_Online'] if 'RR_GLOBAL_Online' in R else [False,'TODO','']
    (expr_t, prompt_t, central) = ([False,'WAIT',''], [False,'WAIT',''], [False,'WAIT',''])
    #if 'GUI_Express' in R:
    if not 'RR_TRACKER_Express' in R:
        if isExpressDoneInGUI(r):
            expr_t = [ False, 'TODO','' ]
    if 'RR_TRACKER_Express' in R:
        expr_t = R['RR_TRACKER_Express'] if 'RR_TRACKER_Express' in R else [False,'TODO',''];
        if options.cosmics:
            #expr_t = [ False, 'N/A','' ]
            print "COSMICS" , expr_t
        elif expr_t[1] == 'TODO' and not isExpressDoneInGUI(r):
             expr_t = [ False, 'Wait','Express not complete in GUI yet' ]
    print 'EXPRT' , expr_t         
    if not options.cosmics and (expr_t[1] == 'Wait' or expr_t[1] == 'WAIT'): continue #ignore collision runs in the waiting list
    #if 'GUI_Prompt' in R:
    if 'RR_TRACKER_Prompt' in R:
        prompt_t = R['RR_TRACKER_Prompt'] if 'RR_TRACKER_Prompt' in R else [False,'TODO',''];
        All_comments+= prompt_t[1]
        central = R['RR_GLOBAL_Prompt']  if 'RR_GLOBAL_Prompt'  in R else [False,'TODO',''];
    note = notes[r] if r in notes else All_comments
    print prompt_t
    print "%6d |  %-15s | %-15s | %-15s | %-15s | %s " % (r, online[1], expr_t[1], prompt_t[1], central[1], note)
    if not options.cosmics:
        html += "<tr><th>%d</th><td class='num'>%.1f T</td><td class='num'>%d</td><td class='num'>%.1f pb<sup>-1</sup></td>" % (r, runlist[r]['RR_bfield'] , lumiCache[r][0], lumiCache[r][1])
    else:
        if lumiCache[r][0] >= 0:
            #html += "<tr><th>%d</th><td class='num'>%.1f T</td><td class='num'>%d</td><td class='num'>%d</td><td class='num'>%.1f Hz</td>" % (r, runlist[r]['RR_bfield'], lumiCache[r][0], lumiCache[r][1], lumiCache[r][1]/lumiCache[r][0]/23.31 )
            #FIXME 
            html += "<tr><th>%d</th><td class='num'>%.1f T</td><td class='num'>%d</td><td class='num'>%d</td><td class='num'>%.1f Hz</td>" % (r, 0.018823952620, lumiCache[r][0], lumiCache[r][1], lumiCache[r][1]/lumiCache[r][0]/23.31 )
        else:
            html += "<tr><th>%d</th><td class='num'>%.1f T</td><td class='num TODO'>%d</td><td class='num TODO'>%d</td><td class='num TODO'>%.1f Hz</td>" % (r,0.018823952620, -lumiCache[r][0], lumiCache[r][1], -lumiCache[r][1]/lumiCache[r][0]/23.31)
            #html += "<tr><th>%d</th><td class='num'>%.1f T</td><td class='num TODO'>%d</td><td class='num TODO'>%d</td><td class='num TODO'>%.1f Hz</td>" % (r, runlist[r]['RR_bfield'], -lumiCache[r][0], lumiCache[r][1], -lumiCache[r][1]/lumiCache[r][0]/23.31)
        html += "<td class='%s'><span title='%s'>%s</span></td>" % (lumiCache[r][3], lumiCache[r][4], lumiCache[r][2])

    if not options.cosmics:
        for X in (online, expr_t, prompt_t, central):
            html += "<td class='%s'>%s</td>" % (v2c(X[0],X[1]), p2t(X))
    else:
        position=0
        for X in (online, expr_t , prompt_t, central):
            html += "<td class='%s'>%s</td>" % (v2c(X[0],X[1]), p2t(X))
            position=position+1
            if position == 3 and options.cosmics and abs(lumiCache[r][0]) > 1: #('BAD' not in X[1])
                #if r >= 290103 and r < 293416 and runlist[r]['RR_bfield'] < 3.6: #0T cosmics (CRUZET)  
                if   0.018823952620 < 3.6: #0T cosmics (CRUZET)               
                        allLumi_currentCRUZET=allLumi_currentCRUZET+abs(lumiCache[r][0])
                        allAlcaTracks_currentCRUZET=allAlcaTracks_currentCRUZET+abs(lumiCache[r][1])
                        info_run_track_CRUZET.append([r, lumiCache[r][1]])
                        info_run_time_CRUZET.append([r, lumiCache[r][0]])
	        #elif r >= 293416 and runlist[r]['RR_bfield'] >= 3.6: #3.8T cosmics (CRAFT)
                #FIXME
	        elif  0.018823952620 >= 3.6: #3.8T cosmics (CRAFT)
                        allLumi_currentCRAFT=allLumi_currentCRAFT+abs(lumiCache[r][0])
                        allAlcaTracks_currentCRAFT=allAlcaTracks_currentCRAFT+abs(lumiCache[r][1])
                        maxcosmicrunforstat = max(maxcosmicrunforstat, r) #should only be in the current era
                        info_run_track_CRAFT.append([r, lumiCache[r][1]])
                        info_run_time_CRAFT.append([r, lumiCache[r][0]])

                #if r >= 272118 and r <= 275418 and runlist[r]['RR_bfield'] > 3.6: #3.8T cosmics                
                #        allLumiB=allLumiB+abs(lumiCache[r][0])
                #        allAlcaTracksB=allAlcaTracksB+abs(lumiCache[r][1])

                        #if lumiCache[r][2] == 'PEAK':
                        #    allAlcaTracksPEAK = allAlcaTracksPEAK + abs(lumiCache[r][1])
                    #else:
                    #    allLumiWait=allLumiWait+abs(lumiCache[r][0])
                    #    allTracksWait=allTracksWait+abs(lumiCache[r][1])

    html += "<td>%s</td></tr>\n" % note;

html += "</table></body></html>"

print 

out = open("statusTestBibhu.%s.html" % groupName, "w")
out.write(html.encode('utf-8')) #prevent crashes when special chars somehow enter description in RR
out.close()

if options.cosmics: 
    #print "total lumi: " , allLumi , " ALCA tracks: " , allAlcaTracks , " hours: " , allLumi * 23.31 / 3600.
    #print "lumi tracks WAIT: " , allLumiWait , " " , allTracksWait

    htmlCOSMICTRACKS = """
        <!DOCTYPE html>
        <html lang="en">
          <head>
            <meta charset="utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1">
            <title>Cosmic Tracks Summary</title>

            <!-- Bootstrap -->
            <link rel="stylesheet" type="text/css" href="css/bootstrap.min.css">

            <!-- Main Style -->
            <link rel="stylesheet" type="text/css" href="css/main.css">
          </head>

        <body>
        <section id="text-about">
            2018 ALCARECO cosmic tracks (%s):
        </section>
        <section id="my-table">
            <div class="container">

                <div class="row">
                    <div class="main">

                        <div class="col-md-4 col-sm-12 col-xs-12">
                            <div class="my-table">
                                <div class="table-header">
                                    <p class="table-title">CRUZET (<a href="/event_display/RunList/TrendPlots/CRUZETBibhu.png">trend</a>, <a
					    href="/event_display/RunList/TrendPlots/CRUZETBibhu_rate.png">rate</a> and <a href="/event_display/RunList/TrendPlots/CRUZETBibhu_cumul.png">cumul</a>)</p>
                                    <p class="table-tracks"><sup>ALCA tracks</sup> 4700K <span>@ 0T</span></p>
                                </div>

                                <div class="table-details">
                                    <ul>
                                        <li>run range: 290103 - 293416</li>
                                        <li>533 hours</li>
                                    </ul>
                                </div>
                            </div>
                        </div>



                        <div class="col-md-4 col-sm-12 col-xs-12">
                            <div class="my-table">
                                <div class="table-header">
                                    <p class="table-title">CRAFT (<a href="/event_display/RunList/TrendPlots/CRAFT.png">trend</a>, <a href="/event_display/RunList/TrendPlots/CRAFT_rate.png">rate</a> and <a href="/event_display/RunList/TrendPlots/CRAFT_cumul.png">cumul</a>)</p>
                                    <p class="table-tracks"><sup>ALCA tracks</sup> %.0fK <span>@ 3.8T</span></p>
                                </div>

                                <div class="table-details">
                                    <ul>
                                        <li>run range: 293416 - %i</li>
                                        <li>%i hours and counting...</li>
                                    </ul>
                                </div>
                            </div>
                        </div>






                    </div>
                </div>

            </div>
        </section>

        </body>
        </html>

""" % (time.ctime() , allAlcaTracks_currentCRAFT / 1000. , maxcosmicrunforstat , abs(allLumi_currentCRAFT * 23.31 / 3600.) )
    outCOSMICTRACKS = open("Cosmics18testBibhu.hours.html", "w")
    outCOSMICTRACKS.write(htmlCOSMICTRACKS)
    outCOSMICTRACKS.close()

#print some info
if options.cosmics:
    print "CRUZET 2018 info: " , allAlcaTracks_currentCRUZET / 1000. , abs(allLumi_currentCRUZET * 23.31 / 3600.)
    print "CRAFT 2018 info: " , allAlcaTracks_currentCRAFT / 1000. , abs(allLumi_currentCRAFT * 23.31 / 3600.)

if options.cosmics:
    #print "B " , allAlcaTracksB / 1000. , abs(allLumiB * 23.31 / 3600.)
    
    #Add a trend plot for the number of tracks in cosmics
    #FIXME This is quite dirty, like this full script for the moment... One should add just a list of name and global variable arrays at the top of the script
    #One should work on this when having time...
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #TRENDS FOR CRUZET
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #
    ## Number of tracks per run
    #
    c_CRUZET = TCanvas("c_CRUZET","c_CRUZET",1,1,1800,800)
    c_CRUZET.SetGridx(True)
    c_CRUZET.SetGridy(True)
    
    Rleg=TLegend(0.45,0.70,0.9,0.8)
    #Rleg.SetHeader("#splitline{CRUZET tracks: "+str(allAlcaTracks_currentCRUZET / 1000.)+"K}{#splitline{Duration: "+str(abs(allLumi_currentCRUZET * 23.31 / 3600.))+" hours}{Rate: "+str(allAlcaTracks_currentCRUZET/(allLumi_currentCRUZET*23.31))+"Hz}}")
    truncated_numOfTracks=allAlcaTracks_currentCRUZET / 1000.
    truncated_numOfTracks= "%.2f" % truncated_numOfTracks
    truncated_time=allLumi_currentCRUZET * 23.31 / 3600.
    truncated_time= "%.2f" % truncated_time
    truncated_frequency = allAlcaTracks_currentCRUZET /( allLumi_currentCRUZET*23.31)
    truncated_frequency= "%.2f" % truncated_frequency
    #Rleg_cumul.SetHeader("#splitline{CRUZET tracks: "+str(allAlcaTracks_currentCRUZET / 1000.)+"K}{Duration: "+str(abs(allLumi_currentCRUZET * 23.31 / 3600.))+" hours}")
    Rleg.SetHeader("#splitline{CRUZET tracks: "+str(truncated_numOfTracks)+"K}{#splitline{Duration: "+str(truncated_time)+" hours}{Rate: "+str(truncated_frequency)+" Hz}}")
   

    Rleg.SetFillStyle(0)
    Rleg.SetBorderSize(0)
    Rleg.SetTextSize(0.08)
    
    trendCRUZET=TH1F("CRUZET","CRUZET",len(info_run_track_CRUZET),0,len(info_run_track_CRUZET))
    trendCRUZET.GetYaxis().SetTitle("# of alca tracks")
    gStyle.SetHistFillColor(4)
    
    tickSpacing = round(len(info_run_track_CRUZET)/50., 0)
    tickCounter=-1
    for i in reversed(xrange(len(info_run_track_CRUZET))):
        ibin = len(info_run_track_CRUZET)-i #careful, ibin is to plot the run numbers in the right order, while i and r are here to find the correct number of tracks associated to a run number
        trendCRUZET.SetBinContent(ibin,info_run_track_CRUZET[i][1]);
        if len(info_run_track_CRUZET) > 50:
            tickCounter += 1
	    if tickCounter % tickSpacing == 0:
                trendCRUZET.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRUZET[i][0]))
	    else:
	        trendCRUZET.GetXaxis().SetBinLabel(ibin,"")
        else:
            trendCRUZET.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRUZET[i][0]))
    
    trendCRUZET.GetXaxis().SetRange(1, len(info_run_track_CRUZET))
    trendCRUZET.GetXaxis().LabelsOption("v")
    trendCRUZET.GetYaxis().SetRangeUser(0,trendCRUZET.GetMaximum()*1.2)
    trendCRUZET.GetYaxis().SetTitleOffset(0.7)
    trendCRUZET.SetStats(0)
    trendCRUZET.Draw()
    Rleg.Draw()
    c_CRUZET.SaveAs("CRUZET.png")

    #
    ##cumulative trend
    #
    c_CRUZET_cumul = TCanvas("c_CRUZET_cumul","c_CRUZET_cumul",1,1,1800,800)
    c_CRUZET_cumul.SetGridx(True)
    c_CRUZET_cumul.SetGridy(True)
    
    Rleg_cumul=TLegend(0.18,0.75,0.48,0.85)
    truncated_numOfTracks=allAlcaTracks_currentCRUZET / 1000.
    truncated_numOfTracks= "%.2f" % truncated_numOfTracks
    truncated_time=allLumi_currentCRUZET * 23.31 / 3600.
    truncated_time= "%.2f" % truncated_time
    truncated_frequency = allAlcaTracks_currentCRUZET /( allLumi_currentCRUZET*23.31)
    truncated_frequency= "%.2f" % truncated_frequency
    #Rleg_cumul.SetHeader("#splitline{CRUZET tracks: "+str(allAlcaTracks_currentCRUZET / 1000.)+"K}{Duration: "+str(abs(allLumi_currentCRUZET * 23.31 / 3600.))+" hours}")
    Rleg_cumul.SetHeader("#splitline{CRUZET tracks: "+str(truncated_numOfTracks)+"K}{#splitline{Duration: "+str(truncated_time)+" hours}{Rate: "+str(truncated_frequency)+" Hz}}")
    Rleg_cumul.SetFillStyle(0)
    Rleg_cumul.SetBorderSize(0)
    Rleg_cumul.SetTextSize(0.08)
    
    trendCRUZET_cumul=TH1F("CRUZET_cumul","CRUZET_cumul",len(info_run_track_CRUZET),0,len(info_run_track_CRUZET))
    trendCRUZET_cumul.GetYaxis().SetTitle("Cumulative sum of alca tracks")
    gStyle.SetHistFillColor(4)
    
    cumulative_tracks=0.
    tickSpacing = round(len(info_run_track_CRUZET)/50., 0)
    tickCounter=-1
    for i in reversed(xrange(len(info_run_track_CRUZET))):
        ibin = len(info_run_track_CRUZET)-i #careful, ibin is to plot the run numbers in the right order, while i and r are here to find the correct number of tracks associated to a run number
        cumulative_tracks+=info_run_track_CRUZET[i][1]
        trendCRUZET_cumul.SetBinContent(ibin,cumulative_tracks);
        trendCRUZET_cumul.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRUZET[i][0]))
        if len(info_run_track_CRUZET) > 50:
            tickCounter += 1
	    if tickCounter % tickSpacing == 0:
                trendCRUZET_cumul.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRUZET[i][0]))
	    else:
	        trendCRUZET_cumul.GetXaxis().SetBinLabel(ibin,"")
        else:
            trendCRUZET_cumul.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRUZET[i][0]))
   
    trendCRUZET_cumul.GetXaxis().SetRange(1, len(info_run_track_CRUZET))
    trendCRUZET_cumul.GetXaxis().LabelsOption("v")
    trendCRUZET_cumul.GetYaxis().SetRangeUser(0,trendCRUZET_cumul.GetMaximum()*1.2)
    trendCRUZET_cumul.GetYaxis().SetTitleOffset(0.7)
    trendCRUZET_cumul.SetStats(0)
    trendCRUZET_cumul.SetLineColor(4)
    trendCRUZET_cumul.Draw()
    Rleg_cumul.Draw()
    c_CRUZET_cumul.SaveAs("CRUZET_cumul.png")

    #
    ## Rate trend
    #
    c_CRUZET_rate = TCanvas("c_CRUZET_rate","c_CRUZET_rate",1,1,1800,800)
    c_CRUZET_rate.SetGridx(True)
    c_CRUZET_rate.SetGridy(True)
    
    Rleg_rate=TLegend(0.16,0.79,0.46,0.94)
    #Rleg.SetHeader("#splitline{CRUZET tracks: "+str(allAlcaTracks_currentCRUZET / 1000.)+"K}{#splitline{Duration: "+str(abs(allLumi_currentCRUZET * 23.31 / 3600.))+" hours}{Rate: "+str(allAlcaTracks_currentCRUZET/(allLumi_currentCRUZET*23.31))+"Hz}}")
    truncated_numOfTracks=allAlcaTracks_currentCRUZET / 1000.
    truncated_numOfTracks= "%.2f" % truncated_numOfTracks
    truncated_time=allLumi_currentCRUZET * 23.31 / 3600.
    truncated_time= "%.2f" % truncated_time
    truncated_frequency = allAlcaTracks_currentCRUZET /( allLumi_currentCRUZET*23.31)
    truncated_frequency= "%.2f" % truncated_frequency
    #Rleg_cumul.SetHeader("#splitline{CRUZET tracks: "+str(allAlcaTracks_currentCRUZET / 1000.)+"K}{Duration: "+str(abs(allLumi_currentCRUZET * 23.31 / 3600.))+" hours}")
    Rleg_rate.SetHeader("#splitline{CRUZET tracks: "+str(truncated_numOfTracks)+"K}{#splitline{Duration: "+str(truncated_time)+" hours}{Rate: "+str(truncated_frequency)+" Hz}}")
   

    Rleg_rate.SetFillStyle(0)
    Rleg_rate.SetBorderSize(0)
    Rleg_rate.SetTextSize(0.04)
    
    trendCRUZET_rate=TH1F("CRUZET_rate","CRUZET_rate",len(info_run_track_CRUZET),0,len(info_run_track_CRUZET))
    trendCRUZET_rate.GetYaxis().SetTitle("Rate of alca tracks (Hz)")
    gStyle.SetHistFillColor(4)
    trendCRUZET_rate_processing=TH1F("CRUZET_rate_processing","CRUZET_rate_processing",len(info_run_track_CRUZET),0,len(info_run_track_CRUZET))
    
    tickSpacing = round(len(info_run_track_CRUZET)/50., 0)
    tickCounter=-1
    for i in reversed(xrange(len(info_run_track_CRUZET))):
        ibin = len(info_run_track_CRUZET)-i #careful, ibin is to plot the run numbers in the right order, while i and r are here to find the correct number of tracks associated to a run number
        if info_run_time_CRUZET[i][1] < 0:
	    trendCRUZET_rate_processing.SetBinContent(ibin,info_run_track_CRUZET[i][1]/(abs(info_run_time_CRUZET[i][1])*23.31)); #The time is negative if the run is not fully processed yet
	else:
	    trendCRUZET_rate.SetBinContent(ibin,info_run_track_CRUZET[i][1]/(info_run_time_CRUZET[i][1]*23.31)); 
	if len(info_run_track_CRUZET) > 50:
            tickCounter += 1
	    if tickCounter % tickSpacing == 0:
                trendCRUZET_rate.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRUZET[i][0]))
	    else:
	        trendCRUZET_rate.GetXaxis().SetBinLabel(ibin,"")
        else:
            trendCRUZET_rate.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRUZET[i][0]))
    
    trendCRUZET_rate.GetXaxis().SetRange(1, len(info_run_track_CRUZET))
    trendCRUZET_rate.GetXaxis().LabelsOption("v")
    trendCRUZET_rate.GetYaxis().SetRangeUser(0,trendCRUZET_rate.GetMaximum()*1.2)
    trendCRUZET_rate.GetYaxis().SetTitleOffset(0.7)
    trendCRUZET_rate.SetStats(0)
    trendCRUZET_rate.Draw()
    trendCRUZET_rate_processing.SetLineColor(8)
    trendCRUZET_rate_processing.SetFillColor(8)
    trendCRUZET_rate_processing.Draw("same")
    
    legend_info = TLegend(.45,.80,.90,.95);
    legend_info.SetHeader("Green points are being processed. Rest is prompt reco.");
    legend_info.Draw()
    legend_info.SetFillStyle(0)
    legend_info.SetBorderSize(0)
    legend_info.SetTextSize(0.04)

    Rleg_rate.Draw()
    c_CRUZET_rate.SaveAs("CRUZET_rate.png")


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #TRENDS FOR CRAFT
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #
    ## Number of tracks per run
    #

    c_CRAFT = TCanvas("c_CRAFT","c_CRAFT",1,1,1800,800)
    c_CRAFT.SetGridx(True)
    c_CRAFT.SetGridy(True)
    
    Rleg=TLegend(0.45,0.70,0.9,0.8)
    #Rleg.SetHeader("#splitline{CRAFT tracks: "+str(allAlcaTracks_currentCRAFT / 1000.)+"K}{#splitline{Duration: "+str(abs(allLumi_currentCRAFT * 23.31 / 3600.))+" hours}{Rate: "+str(allAlcaTracks_currentCRAFT/(allLumi_currentCRAFT*23.31))+"Hz}}")
    truncated_numOfTracks=allAlcaTracks_currentCRAFT / 1000.
    truncated_numOfTracks= "%.2f" % truncated_numOfTracks
    truncated_time=allLumi_currentCRAFT * 23.31 / 3600.
    truncated_time= "%.2f" % truncated_time
    truncated_frequency = allAlcaTracks_currentCRAFT /( allLumi_currentCRAFT*23.31)
    truncated_frequency= "%.2f" % truncated_frequency
    #Rleg_cumul.SetHeader("#splitline{CRAFT tracks: "+str(allAlcaTracks_currentCRAFT / 1000.)+"K}{Duration: "+str(abs(allLumi_currentCRAFT * 23.31 / 3600.))+" hours}")
    Rleg.SetHeader("#splitline{CRAFT tracks: "+str(truncated_numOfTracks)+"K}{#splitline{Duration: "+str(truncated_time)+" hours}{Rate: "+str(truncated_frequency)+" Hz}}")
   

    Rleg.SetFillStyle(0)
    Rleg.SetBorderSize(0)
    Rleg.SetTextSize(0.08)
    
    trendCRAFT=TH1F("CRAFT","CRAFT",len(info_run_track_CRAFT),0,len(info_run_track_CRAFT))
    trendCRAFT.GetYaxis().SetTitle("# of alca tracks")
    gStyle.SetHistFillColor(4)
    
    tickSpacing = round(len(info_run_track_CRAFT)/50., 0)
    tickCounter=-1
    for i in reversed(xrange(len(info_run_track_CRAFT))):
        ibin = len(info_run_track_CRAFT)-i #careful, ibin is to plot the run numbers in the right order, while i and r are here to find the correct number of tracks associated to a run number
        trendCRAFT.SetBinContent(ibin,info_run_track_CRAFT[i][1]);
        if len(info_run_track_CRAFT) > 50:
            tickCounter += 1
	    if tickCounter % tickSpacing == 0:
                trendCRAFT.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRAFT[i][0]))
	    else:
	        trendCRAFT.GetXaxis().SetBinLabel(ibin,"")
        else:
            trendCRAFT.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRAFT[i][0]))
    
    trendCRAFT.GetXaxis().SetRange(1, len(info_run_track_CRAFT))
    trendCRAFT.GetXaxis().LabelsOption("v")
    trendCRAFT.GetYaxis().SetRangeUser(0,trendCRAFT.GetMaximum()*1.2)
    trendCRAFT.GetYaxis().SetTitleOffset(0.7)
    trendCRAFT.SetStats(0)
    trendCRAFT.Draw()
    Rleg.Draw()
    c_CRAFT.SaveAs("CRAFT.png")

    #
    ##cumulative trend
    #
    c_CRAFT_cumul = TCanvas("c_CRAFT_cumul","c_CRAFT_cumul",1,1,1800,800)
    c_CRAFT_cumul.SetGridx(True)
    c_CRAFT_cumul.SetGridy(True)
    
    Rleg_cumul=TLegend(0.18,0.75,0.48,0.85)
    truncated_numOfTracks=allAlcaTracks_currentCRAFT / 1000.
    truncated_numOfTracks= "%.2f" % truncated_numOfTracks
    truncated_time=allLumi_currentCRAFT * 23.31 / 3600.
    truncated_time= "%.2f" % truncated_time
    truncated_frequency = allAlcaTracks_currentCRAFT /( allLumi_currentCRAFT*23.31)
    truncated_frequency= "%.2f" % truncated_frequency
    #Rleg_cumul.SetHeader("#splitline{CRAFT tracks: "+str(allAlcaTracks_currentCRAFT / 1000.)+"K}{Duration: "+str(abs(allLumi_currentCRAFT * 23.31 / 3600.))+" hours}")
    Rleg_cumul.SetHeader("#splitline{CRAFT tracks: "+str(truncated_numOfTracks)+"K}{#splitline{Duration: "+str(truncated_time)+" hours}{Rate: "+str(truncated_frequency)+" Hz}}")
    Rleg_cumul.SetFillStyle(0)
    Rleg_cumul.SetBorderSize(0)
    Rleg_cumul.SetTextSize(0.08)
    
    trendCRAFT_cumul=TH1F("CRAFT_cumul","CRAFT_cumul",len(info_run_track_CRAFT),0,len(info_run_track_CRAFT))
    trendCRAFT_cumul.GetYaxis().SetTitle("Cumulative sum of alca tracks")
    gStyle.SetHistFillColor(4)
    
    cumulative_tracks=0.
    tickSpacing = round(len(info_run_track_CRAFT)/50., 0)
    tickCounter=-1
    for i in reversed(xrange(len(info_run_track_CRAFT))):
        ibin = len(info_run_track_CRAFT)-i #careful, ibin is to plot the run numbers in the right order, while i and r are here to find the correct number of tracks associated to a run number
        cumulative_tracks+=info_run_track_CRAFT[i][1]
        trendCRAFT_cumul.SetBinContent(ibin,cumulative_tracks);
        trendCRAFT_cumul.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRAFT[i][0]))
        if len(info_run_track_CRAFT) > 50:
            tickCounter += 1
	    if tickCounter % tickSpacing == 0:
                trendCRAFT_cumul.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRAFT[i][0]))
	    else:
	        trendCRAFT_cumul.GetXaxis().SetBinLabel(ibin,"")
        else:
            trendCRAFT_cumul.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRAFT[i][0]))
   
    trendCRAFT_cumul.GetXaxis().SetRange(1, len(info_run_track_CRAFT))
    trendCRAFT_cumul.GetXaxis().LabelsOption("v")
    trendCRAFT_cumul.GetYaxis().SetRangeUser(0,trendCRAFT_cumul.GetMaximum()*1.2)
    trendCRAFT_cumul.GetYaxis().SetTitleOffset(0.7)
    trendCRAFT_cumul.SetStats(0)
    trendCRAFT_cumul.SetLineColor(4)
    trendCRAFT_cumul.Draw()
    Rleg_cumul.Draw()
    c_CRAFT_cumul.SaveAs("CRAFT_cumul.png")

    #
    ## Rate trend
    #
    c_CRAFT_rate = TCanvas("c_CRAFT_rate","c_CRAFT_rate",1,1,1800,800)
    c_CRAFT_rate.SetGridx(True)
    c_CRAFT_rate.SetGridy(True)
    
    Rleg_rate=TLegend(0.16,0.79,0.46,0.94)
    #Rleg.SetHeader("#splitline{CRAFT tracks: "+str(allAlcaTracks_currentCRAFT / 1000.)+"K}{#splitline{Duration: "+str(abs(allLumi_currentCRAFT * 23.31 / 3600.))+" hours}{Rate: "+str(allAlcaTracks_currentCRAFT/(allLumi_currentCRAFT*23.31))+"Hz}}")
    truncated_numOfTracks=allAlcaTracks_currentCRAFT / 1000.
    truncated_numOfTracks= "%.2f" % truncated_numOfTracks
    truncated_time=allLumi_currentCRAFT * 23.31 / 3600.
    truncated_time= "%.2f" % truncated_time
    truncated_frequency = allAlcaTracks_currentCRAFT /( allLumi_currentCRAFT*23.31)
    truncated_frequency= "%.2f" % truncated_frequency
    #Rleg_cumul.SetHeader("#splitline{CRAFT tracks: "+str(allAlcaTracks_currentCRAFT / 1000.)+"K}{Duration: "+str(abs(allLumi_currentCRAFT * 23.31 / 3600.))+" hours}")
    Rleg_rate.SetHeader("#splitline{CRAFT tracks: "+str(truncated_numOfTracks)+"K}{#splitline{Duration: "+str(truncated_time)+" hours}{Rate: "+str(truncated_frequency)+" Hz}}")
   

    Rleg_rate.SetFillStyle(0)
    Rleg_rate.SetBorderSize(0)
    Rleg_rate.SetTextSize(0.04)
    
    trendCRAFT_rate=TH1F("CRAFT_rate","CRAFT_rate",len(info_run_track_CRAFT),0,len(info_run_track_CRAFT))
    trendCRAFT_rate.GetYaxis().SetTitle("Rate of alca tracks (Hz)")
    gStyle.SetHistFillColor(4)
    trendCRAFT_rate_processing=TH1F("CRAFT_rate_processing","CRAFT_rate_processing",len(info_run_track_CRAFT),0,len(info_run_track_CRAFT))
    
    tickSpacing = round(len(info_run_track_CRAFT)/50., 0)
    tickCounter=-1
    for i in reversed(xrange(len(info_run_track_CRAFT))):
        ibin = len(info_run_track_CRAFT)-i #careful, ibin is to plot the run numbers in the right order, while i and r are here to find the correct number of tracks associated to a run number
        if info_run_time_CRAFT[i][1] < 0:
	    trendCRAFT_rate_processing.SetBinContent(ibin,info_run_track_CRAFT[i][1]/(abs(info_run_time_CRAFT[i][1])*23.31)); #The time is negative if the run is not fully processed yet
	else:
	    trendCRAFT_rate.SetBinContent(ibin,info_run_track_CRAFT[i][1]/(info_run_time_CRAFT[i][1]*23.31)); 
	if len(info_run_track_CRAFT) > 50:
            tickCounter += 1
	    if tickCounter % tickSpacing == 0:
                trendCRAFT_rate.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRAFT[i][0]))
	    else:
	        trendCRAFT_rate.GetXaxis().SetBinLabel(ibin,"")
        else:
            trendCRAFT_rate.GetXaxis().SetBinLabel(ibin,str(info_run_track_CRAFT[i][0]))
    
    trendCRAFT_rate.GetXaxis().SetRange(1, len(info_run_track_CRAFT))
    trendCRAFT_rate.GetXaxis().LabelsOption("v")
    trendCRAFT_rate.GetYaxis().SetRangeUser(0,trendCRAFT_rate.GetMaximum()*1.2)
    trendCRAFT_rate.GetYaxis().SetTitleOffset(0.7)
    trendCRAFT_rate.SetStats(0)
    trendCRAFT_rate.Draw()
    trendCRAFT_rate_processing.SetLineColor(8)
    trendCRAFT_rate_processing.SetFillColor(8)
    trendCRAFT_rate_processing.Draw("same")
    
    legend_info = TLegend(.45,.80,.90,.95);
    legend_info.SetHeader("Green points are being processed. Rest is prompt reco.");
    legend_info.Draw()
    legend_info.SetFillStyle(0)
    legend_info.SetBorderSize(0)
    legend_info.SetTextSize(0.04)

    Rleg_rate.Draw()
    c_CRAFT_rate.SaveAs("CRAFT_rate.png")
