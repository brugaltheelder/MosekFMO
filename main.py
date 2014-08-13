__author__ = 'troy'


import sys
import mosek
import numpy as np
import scipy.io as io
import scipy.sparse as sps

inf = 0.0
#inputFilename = 'pyfmoNPC5.mat'
#inputFilename = 'pyfmoMGH.mat'
inputFilename = 'pyfmoMGH6.mat'
#inputFilename = 'pyfmoMGH12.mat'
print 'Input Filename = ',inputFilename

class dao_data(object):
    def __init__(self,inputFilename):
        #Reads in map file
        matFile = io.loadmat(inputFilename)
        self.nVox,self.beams,self.dataFolder = int(matFile['nVox']), np.array(matFile['beams'].flatten()),str(matFile['dataFolder'][0])
        self.nBeams,self.nBPB,self.nDIJSPB = len(self.beams),np.array(matFile['nBPB'].flatten()),np.array(matFile['nDIJSPB'].flatten())
        self.aOver, self.aUnder, self.thresh = np.array(matFile['aOver'].flatten()),np.array(matFile['aUnder'].flatten()),np.array(matFile['thresh'].flatten())
        self.currentDose, self.objDescent,self.aperData= np.array([]),[],[]
        #self.nBeams = 1 #uncomment this to limit beams used for more efficient debugging        
        self.numVars = sum(self.nBPB[x] for x in xrange(self.nBeams)) + 3*self.nVox
        self.numCons = 3*self.nVox
        

def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()

print 'starting mosek environment'
env = mosek.Env()
task = env.Task()

print 'initializing up output stream'
task.set_Stream(mosek.streamtype.log,streamprinter)

# read in matlab file to get numvoxels and numbixels
print 'reading in data'
data = dao_data(inputFilename)

#build bounds
print 'building bounds'
#constraint bounds
bkc = [None]*data.numCons
buc = np.zeros(data.numCons)
blc = np.zeros(data.numCons)

#variable bounds
bkx = [mosek.boundkey.lo]*data.numVars
bux = np.zeros(data.numVars)
blx = np.zeros(data.numVars)


for i in range(data.nVox):
    bkc[i] = mosek.boundkey.fx
    bkx[i] = mosek.boundkey.fr
for i in range(data.nVox,2*data.nVox):
    bkc[i] = mosek.boundkey.up
    buc[i] = data.thresh[i-data.nVox]
    blc[i] = -inf
for i in range(2*data.nVox,3*data.nVox):
    bkc[i] = mosek.boundkey.lo
    buc[i] = inf
    blc[i] = data.thresh[i-2*data.nVox]


asub = [[] for i in xrange(data.numVars)]
aval = [[] for i in xrange(data.numVars)]

#populate A for z, over, under
print 'building A matrix'
for i in xrange(data.nVox):
    #populate for z
    #dose constraint
    asub[i].append(i)
    aval[i].append(-1)

    #over constraint
    asub[i].append(i+data.nVox)
    aval[i].append(1)

    #under constraint
    asub[i].append(i+2*data.nVox)
    aval[i].append(1)

    #populate for over
    asub[i+data.nVox].append(i+data.nVox)
    aval[i+data.nVox].append(-1)
    
    #populate for under
    asub[i+2*data.nVox].append(i+2*data.nVox)
    aval[i+2*data.nVox].append(1)

#populate A for beams
for b in xrange(data.nBeams):
    print '\tAdding beam ',b
    start = 3*data.nVox + sum(data.nBPB[i] for i in xrange(b))
    #read in Dij
    Dij = np.fromfile(data.dataFolder +'beam' + str(data.beams[b]) + 'bixDs.bin',dtype = np.float32).reshape((-1,3))
    DijSparse = sps.csr_matrix((Dij[:][:,2],(Dij[:][:,0]-1,Dij[:][:,1]-1)))
    for i in range(data.nBPB[b]):
        asub[start+i] = DijSparse[i,:].indices[:]
        aval[start+i] = DijSparse[i,:].data[:]
    '''
    for (i,j,d) in Dij:
        asub[start+int(i)-1].append(int(j)-1)
        aval[start+int(i)-1].append(d)
    '''

#set up mosek quad object vars (diag through over and under vars)
print 'building objective'
qsubi = np.arange(data.nVox,3*data.nVox)
qsubj = np.arange(data.nVox,3*data.nVox)
qval = np.concatenate((data.aOver,data.aUnder))

# create variables and constraints
print 'adding variables and constraints to task'
task.appendcons(data.numCons)
task.appendvars(data.numVars)


#set variable bounds
print 'setting bounds'
task.putvarboundlist(np.arange(data.numVars),bkx,blx,bux)
task.putconboundlist(np.arange(data.numCons),bkc,blc,buc)

#build A matrix
print 'setting A matrix'
for j in range(data.numVars):
    if j%1000==0:
        sys.stdout.write('\t Be patient... %d%% done   \r' % (100.0*(j+1)/data.numVars))
    task.putacol(j,asub[j],aval[j])
sys.stdout.write('\t Be patient... 100%% done       \n')

#build objective
print 'setting objectve'
task.putqobj(qsubi,qsubj,qval)

#set sense
print 'setting sense'
task.putobjsense(mosek.objsense.minimize)

#print 'writing out model'
#task.writedata('out12.lp')

print 'optimizing model'
task.optimize()

#get z zolution and print to text file for drawing
zValues = np.zeros(data.nVox,float)
task.getxxslice(mosek.soltype.itr,0,data.nVox,zValues)
intensities = np.zeros(sum(data.nBPB[x] for x in xrange(data.nBeams)) ,float)
task.getxxslice(mosek.soltype.itr,data.nVox*3,data.numVars,intensities)
allVars = np.zeros(data.numVars,float)
task.getxxslice(mosek.soltype.itr,0,data.numVars,allVars)

io.savemat('out.mat',{'dose':zValues,'intensities':intensities,'allVars':allVars})

print 'FMO complete'
