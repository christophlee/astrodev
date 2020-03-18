import numpy as np
import pandas as pd
from scipy.ndimage.filters import gaussian_filter1d
import subprocess
import sys
from math import *
import glob

def calcDerivativeAccretion (datafile):

    for rho_bin in [0,1,2]:

        output = open((datafile+'2') % rho_bin, 'w')

        data = np.loadtxt(datafile % rho_bin)

        data[:,2] *= 1.e12  # renormalize masses to center of mass bin

        diff = np.zeros(len(data))

        [z,t,m] = [0,1,2]

        lag = 2 

        for i in range (0,len(data)):

            if i >= len(data)-lag:

                diff[i] = diff[i-1]

            else:

                diff[i] = (data[i][m]-data[i+lag][m])/(data[i][t] - data[i+lag][t])


        for i in range (0,len(data)):
            
            output.write('%f %f %f %f\n' % (data[i][z],data[i][t],data[i][m],diff[i]))

        output.close()

def medianAccretionHist (datafile):

    for rho_bin in [0,1,2]:

        output = open(datafile % rho_bin,'w')

        data_base = '../../tree_code_outputs/bp_z0_centrals_m12_rho_%d_treewalk_%d.dat'

        num_fields = 3

        extracted_data = [{} for a in range(0,num_fields)]

        t = {}

        num_trees = 1100

        dot_freq = np.ceil(num_trees/10.)

        for i in range (0,num_trees):
            if (i % dot_freq == 0):
                sys.stdout.write('.')
                sys.stdout.flush()
            data_infile = np.loadtxt(data_base % (rho_bin,i)).astype(np.float64)
            for j in range (0,len(data_infile[:,1])-1):
                for k in range(0,num_fields):
                    z1 = data_infile[j][1]
                    z2 = data_infile[j+1][1]
                    # normalize accretion histories to have same final mass
                    if k is 0:
                        m = data_infile[j][11]/data_infile[0][11]
                    elif k is 1:
                        for z in [z1,z2]:
                            if z not in t.keys():
                                t[z] = np.array(subprocess.Popen("python /home/christoph/Research/age_matching_comparison/cosmocalc.py %f 67.8 0.307 0.693" % z, shell=True, stdout=subprocess.PIPE).stdout.read()[:-1].split())[0].astype(np.float64)
                        m = (data_infile[j][11]-data_infile[j+1][11])/(t[z1]-t[z2])
                    elif k is 2:
                        m = (data_infile[j][11]-data_infile[j+1][11])/(t[z1]-t[z2])/data_infile[0][11]

                    #m = data_infile[j][46]
                    #m = data_infile[j][36]
                    #m = data_infile[j][12]/data_infile[j][13]
                    if z1 in extracted_data[k].keys():
                        extracted_data[k][z1].append(m)
                    else:
                        extracted_data[k][z1] = [m]

        print ''

        medians = np.zeros((len(extracted_data[0].keys()),num_fields+1))
        medians[:,0] = extracted_data[0].keys()
        medians[:,0].sort()

        for i in range(0,len(medians[:,0])):
            values = [extracted_data[k][medians[i][0]] for k in range(0,num_fields)]
            for field in values:
                field.sort()
            median = [(values[k])[int(len(values[k])/2.)] for k in range(0,num_fields)]
            medians[i,1:] = median
            #print str(medians[i][0]) + " , " + str(medians[i][1])

            output.write("%f %f %f %f %f\n" % (medians[i][0],t[medians[i][0]],medians[i][1],medians[i][2],medians[i][3]))

        output.close()

    return


def genAccretionHistPlot (datafile):

    rho_bin = 2

    output = open(datafile % rho_bin,'w')

    data_base = 'bp_z0_centrals_m12_rho_%d_treewalk_%d.dat'

    output.write("plot '"+data_base % (rho_bin,0) + "' u 2:(log10($12)) w l lw 1 lt 1 lc rgb 'black'")

    for i in range(1,100):
        output.write(", '" + data_base % (rho_bin,i))
        output.write("' u 2:(log10($12)) w l lw 1 lt 1 lc rgb 'black'")

    output.write('\n')

    return


def smooth_data (datafile):

    datafiles = []

    # construct array of datafiles

    datafile = "paper1/plotgen/data/bp_z0_centrals_hmbins_rw%d_%s_%d_med.bin"
    for hmbin in [0,1,2,3]:
        for field in ['rs','rs_jiang']: #['amhalf','almm','tf']:#'cnfw','lambdap','mar','cklypin','lambda','vmax','bsr','p','p500']: #['amhalf']:
            for rho in [0,1,2,4,8,16]:
                datafiles.append(datafile % (rho,field,hmbin))

    #datafile = 'data/bp_z0_centrals_bsr_%s_hmbins_%d_med.bin'
    #for hmbin in [0,1,2,3]:
    #    for field in ['tf','rho','p','almm']:#'cnfw','lambdap']:
    #        datafiles.append(datafile % (field,hmbin))

    #datafile = 'data/bp_z0_centrals_hmbins_rw%s_%s_p_%d_med.bin'
    #for rho in ['0','1','2','4','8','16']:
    #    for field in ['xoff']:
    #        for hmbin in [0,1,2,3]:
    #            datafiles.append(datafile % (rho,field,hmbin))

    #datafile = 'paper2/plotgen/data/bp_z0_centrals_hmbins_%s_mdist_%d_hist1D.bin'
    #for hmbin in ['0','0_rw1','1_rw2','2_rw4','3_rw8']:
    #    for rho in [0,1,2,3,4]:
    #        datafiles.append(datafile % (hmbin,rho))

    #datafile = 'paper2/plotgen/data/bp_z0_centrals_hmbins_%d_bsr_tf_m_%d_almm_%d_%s_hist1D_count.dat'
    #for hmbin in [0,1,2,3]:
    #    for field in ['tf','rho','p','almm','cnfw','lambdap','bsr','xoff','tu']:
    #        for tf in [0,1]:
    #            for almm in [0,1]:
    #                    datafiles.append(datafile % (hmbin,tf,almm,field))
    #        datafiles.append('paper2/plotgen/data/bp_z0_centrals_hmbins_%d_ns_%s_hist1D_count.dat' % (hmbin,field))

    #datafile = "paper2/plotgen/data/bp_z0_centrals_hmbins_0_bsr_tf_m_%d_almm_%dm_mtree2_track_%d.dat"
    #datafile = "paper2/plotgen/data/bp_z0_centrals_hmbins_1_pd_%s_mtree_track_%d.dat"
    #datafile = "paper2/plotgen/data/bp_z0_centrals_hmbins_1_pd_*_mtree_track_*.dat"

    #datafile = "paper2/plotgen/halo_tracks/bp_z0_centrals_hmbins_*_pd_*_mtree_track_*.dat"

    #datafiles = glob.glob(datafile)

    #datafile = "paper2/plotgen/halo_tracks/bp_z0_centrals_all_random_*_mtree2_track_*.dat"

    #datafiles = datafiles + glob.glob(datafile)

    #datafile = "paper2/plotgen/halo_tracks/bp_z0_centrals_1e12_*_mtree2_track_*.dat"

    #datafiles = datafiles + glob.glob(datafile)

    #print datafiles

    #for track in range (0,21):
    #    datafiles.append(datafile % (1,0,track))
    #    datafiles.append(datafile % (0,1,track))

    #for datafile in datafiles:
    #    print datafile
    #    try:
    #        data = np.loadtxt(datafile)
    #    except IOError as e:
    #        print e
    #    data[:,10] = gaussian_filter1d(data[:,10],1.0)
    #    data[:,11] = gaussian_filter1d(data[:,11],1.0)
    #    data[:,12] = gaussian_filter1d(data[:,12],1.0)
    #    data[:,16] = gaussian_filter1d(data[:,16],1.0)
    #    data[:,35] = gaussian_filter1d(data[:,35],1.0)
    #    data[:,43] = gaussian_filter1d(data[:,43],1.0)
    #    data[:,45] = gaussian_filter1d(data[:,45],1.0)
    #    data[:,56] = gaussian_filter1d(data[:,56],1.0)

    #    print "writing prolateness to column ", data.shape[1]

    #    data2 = np.zeros((data.shape[0],data.shape[1]+1))
    #    data2[:,:-1] = data
    #    data2[:,-1] = 1-np.sqrt(data[:,46]*data[:,46]+data[:,47]*data[:,47])/sqrt(2.)
    #    data2[:,-1] = gaussian_filter1d(data2[:,-1],1.0)

    #    np.savetxt (datafile+'2',data2,fmt="%.6f")

    #return

    
    #for hmbin in ['0_rw1','1_rw2','2_rw4','3_rw8']:
    #    for field in ['cnfw','lambdap']:
    #        for rho in [0,1]:
    #            datafiles.append(datafile % (hmbin,rho,field))

    #for hmbin in ['0','1','2','3']:
    #    for field in ['cnfw','lambdap']:
    #        for i in [0,1]:
    #            datafiles.append(datafile % (hmbin,field,i))

    #datafile = 'data/bp_z0_centrals_hmbins_%s_bsr_%s_%d_hist1D.bin'
    #for hmbin in ['0','1','2','3']:
    #    for field in ['tf','rho','p','almm','cnfw','lambdap']:
    #        for i in [0,1,2,3]:
    #            datafiles.append(datafile % (hmbin,field,i))

    for datafile in datafiles:

        data = np.loadtxt(datafile)

        if len(data) == 0:
            continue

        # this would throw except if only one line of data is present
        try:
            (a,b) = data.shape
        except ValueError as e:
            np.savetxt(datafile+'2',np.array([data]))
            continue
        
        # select columns to smooth

        # smooth 1D histograms (pdfs)
        if ("_hist1D" in datafile):
            data[:,1] = gaussian_filter1d(data[:,1],1.0)
            try:
                data[:,4] = gaussian_filter1d(data[:,4],1.0)
            except IndexError as e:
                print "File does not contain 5 columns: "+datafile

        # smooth median output files
        elif ("_med" in datafile):
            # for now, we want to smooth medians, +- CI, 80-20 (cols 1,2,3,5,6)
            data[:,1] = gaussian_filter1d(data[:,1],1.0)
            data[:,2] = gaussian_filter1d(data[:,2],1.0)
            data[:,3] = gaussian_filter1d(data[:,3],1.0)
            data[:,5] = gaussian_filter1d(data[:,5],1.0)
            data[:,6] = gaussian_filter1d(data[:,6],1.0)

        np.savetxt(datafile+'2',data)

    return


def smooth_prog_hist (datafile):

    datafiles = []

    #datafile = 'data/bp_z0_centrals_hmbins_%s_rho_%d_mtree_med.dat'
    #for hmbin in ['0_rw1','1_rw2','2_rw4','3_rw8']:
    #    for rho in [0,1,2,3]:
    #        datafiles.append(datafile % (hmbin,rho))


    hmbins = ['0','1','2','3']
    hmbin = ['hmbins_'+hm for hm in hmbins]
    ns_txt = ['sub_mtree','sub_mtree','mtree','mtree']

    datafile = 'paper2/plotgen/data/bp_z0_centrals_hmbins_%d_bsr_tf_m_%d_almm_%d_%s_med.dat'
    for hm in [0,1,2,3]:
        for field in ['mtree']:
            for tf in [0,1]:
                for almm in [0,1]:
                        datafiles.append(datafile % (hm,tf,almm,field))
            datafiles.append('paper2/plotgen/data/bp_z0_centrals_hmbins_%d_ns_%s_med.dat' % (hm,ns_txt[hm]))

    #datafile = 'data/bp_z0_centrals_hmbins_%s_%d_mtree_med.dat'

    #for hmbin in ['0_bsr','1_bsr','2_bsr','3_bsr']:
    #    for bsr in [0,1,2]:
    #        datafiles.append(datafile % (hmbin,bsr))

    for datafile in datafiles:

        try:
            data = np.loadtxt(datafile)
        except IOError as e:
            print e
            continue

        if len(data) == 0:
            continue

        # now, we want to smooth all columns except for the first column (z).
        # additional, any trailing zeros should be excepted from the smoothing.
        for i in range (1,len(data[0])):

            # let's find where the trailing zeros end. start at the end and work forward
            c = -1
            while (data[c,i] == 0):
                c -= 1

                # whole column is zero
                if c < -1*len(data):
                    break
            c += 1

            # c now tells us how much of the array to include in the smoothing
            #print "Smoothing up to element %d" % c

            if hmbin[0] in datafile:
                sigma = 2
            elif hmbin[1] in datafile:
                sigma = 3
            elif hmbin[2] in datafile:
                sigma = 5.0
            elif hmbin[3] in datafile:
                sigma = 6.0

            if (c == 0):
                data[:,i] = gaussian_filter1d(data[:,i],sigma)
            else:
                data[:c,i] = gaussian_filter1d(data[:c,i],sigma)

        np.savetxt(datafile+'2',data)

    return


def normalize_prog_hist (datafile):

    datafile = 'paper2/plotgen/data/bp_z0_centrals_hmbins_%d_bsr_tf_m_%d_almm_%d_mtree_med.dat2'
    datafile_ns = 'paper2/plotgen/data/bp_z0_centrals_hmbins_%d_ns_%s_med.dat2'

    ns_txt = ['sub_mtree','sub_mtree','mtree','mtree']

    for hmbin in [0,1,2,3]:

        data = []

        datafiles = []

        #####################
        # build custom datafile names here. change as needed.
        # datafile[0] will be used to normalize datafiles[1:].
        datafiles.append(datafile_ns % (hmbin,ns_txt[hmbin]))

        for tf in [0,1]:
            for almm in [0,1]:
                datafiles.append(datafile % (hmbin,tf,almm))
        #####################

        for df in datafiles:

            try:
                data.append(np.loadtxt(df))
            except IOError as e:
                print e

        for (c1,c2) in [(1,6),(16,21),(21,26),(41,46),(56,61),(61,66)]:

            for d in range(1,len(data)):
                data[d][:,c1:c2] /= data[0][0,c1]

            data[0][:,c1:c2] /= data[0][0,c1]

        for d in range(0,len(datafiles)):
            try:
                np.savetxt(datafiles[d]+'n', data[d])
            except IndexError as e:
                print e
    
#    old version used for paper1
#
#    datafile = 'data/bp_z0_centrals_hmbins_%s_rho_%d_mtree_med.dat2'
#
#    for hmbin in ['0_rw1','1_rw2','2_rw4','3_rw8']:
#    #for hmbin in ['0','1','2','3']:
#
#        data = []
#
#        for rho in [0,1,2,3]:
#
#            try:
#                data.append(np.loadtxt(datafile % (hmbin,rho)))
#            except IOError as e:
#                print e
#
#        for (c1,c2) in [(1,6),(16,21),(21,26),(41,46),(56,61),(61,66)]:
#
#            data[0][:,c1:c2] /= data[1][0,c1]
#            data[2][:,c1:c2] /= data[1][0,c1]
#            try:
#                data[3][:,c1:c2] /= data[1][0,c1]
#            except IndexError as e:
#                print e
#
#            data[1][:,c1:c2] /= data[1][0,c1]
#
#        for rho in [0,1,2,3]:
#            try:
#                np.savetxt((datafile+'n') % (hmbin,rho) , data[rho])
#            except IndexError as e:
#                print e

    return


def reconstructHist1D (file_name):

    datafiles = []

    #for hmbin in [0,1,2,3]:
    #    for field in ['amhalf']:
    #        for density in [0,1,2,3]:
    #            datafiles.append(file_name % (hmbin,field,density)

    for hmbin in [0,1,2,3]:
            for nbin in [0,1,2,3]:
                datafiles.append(file_name % (hmbin,nbin))

    for datafile in datafiles:

        data = np.loadtxt(datafile)

        output = open(datafile+'2','w')

        for row in data:
            if row[1] > 0.0:
                output.write('%f %f\n' % (row[0],row[1]))

        output.close()

        data = np.loadtxt(datafile+'2')

        total = 0
        # let's renormalize hist
        try:
            for i in range(1,len(data)):
                total += (data[i][0]-data[i-1][0])*(np.minimum(data[i][1],data[i-1][1])+0.5*np.fabs(data[i][1]-data[i-1][1]))
        except ValueError:
            print "Value Error!!!"

        data[:,1] /= total

        data[:,1] = gaussian_filter1d(data[:,1],2.0)

        np.savetxt(datafile+'2',data)

    return


def rescale_xaxis (datafile, z, percentiles=False):

    for field in ['cnfw','lambdap','mar']:
    
        for mvir in [0,1,2,3]:
        
            file_name = data_file % (z,1,field,mvir) 
            base_scale = np.loadtxt(file_name)
    
            if len(base_scale):
            
                if not percentiles:
                
                    xmin = np.min(base_scale[:,0])
                    xmax = np.max(base_scale[:,0])

                else:

                    xmin = 0.
                    xmax = 100.
    
            np.savetxt(file_name+'2',base_scale)
            
            for scale in [1,2,4,6]:
            
                if not percentiles:
                    if scale is 1:
                        continue

                file_name = data_file % (z,scale,field,mvir)
                data = np.loadtxt(file_name)
    
                if len(data):
                    xmin2 = np.min(data[:,0])
                    xmax2 = np.max(data[:,0])
                    data[:,0] = (data[:,0]-xmin2)*(xmax-xmin)/(xmax2-xmin2) + xmin
    
                np.savetxt(file_name+'2',data)


def rescale_xaxis_stacked (datafile, percentiles=False):

    #for field in ['cnfw','lambdap','mar','mar2']:
    for field in ['cnfw','mar','mar2']:
    
        for mvir in [0,1,2,3,4]:
        
            file_name = data_file % (mvir,field,1) 
            base_scale = np.loadtxt(file_name)
    
            if len(base_scale):
            
                if not percentiles:
                
                    xmin = np.min(base_scale[:,0])
                    xmax = np.max(base_scale[:,0])

                else:

                    xmin = 0.
                    xmax = 100.
    
            np.savetxt(file_name+'2',base_scale)
            
            for scale in [1,2,4,6]:
            
                if not percentiles:
                    if scale is 1:
                        continue

                file_name = data_file % (mvir,field,scale)
                data = np.loadtxt(file_name)
    
                if len(data):
                    xmin2 = np.min(data[:,0])
                    xmax2 = np.max(data[:,0])
                    data[:,0] = (data[:,0]-xmin2)*(xmax-xmin)/(xmax2-xmin2) + xmin
    
                np.savetxt(file_name+'2',data)

def separate_mtrees (datafile):

    data = np.loadtxt(datafile)

    datafile = datafile[:-4]+"_%d"+datafile[-4:]

    of_num = 0
    start = 0
    max_files = 40

    for i in range(0,len(data)):
        if data[i][0] == 1.00231:
            if i > 0:

                for j in [10,11,12,35,46,47]:
                    data[start:i,j] = gaussian_filter1d(data[start:i,j],1.0)
                
                print "Writing file %d ..." % of_num
                np.savetxt(datafile % of_num,data[start:i])
                start = i
                of_num += 1
                
                if (of_num == max_files):
                    break

    # get last tree in file
    if (of_num < max_files):
        for j in [10,11,12,35,46,47]:
            data[start:,j] = gaussian_filter1d(data[start:,j],1.0)
        np.savetxt(datafile % of_num, data[start:])

def msmh_relation ():
    def f(x):
        r1 = -log10(pow(10.,a*x)+1.)
        r2 = d*pow(log10(1.+exp(x)),g)
        r2 = r2/(1.+exp(pow(10.,-x)))
        return r1+r2
    def Ms(x):
        return le+lm1+f(log10(x)-lm1)-f(0)
    a = -1.84754     #-1.412
    d = 3.85767      #3.508
    g = 0.44557      #0.316
    le = -1.959715   #-1.777
    lm1 = 11.483698  #11.514
    return Ms

    #h=0.678
    #x = np.arange(10,15,0.1) #msun/h
    #y = [msmh_relation()(10**a)/h for a in x] #msun/h
    #y2 = [(10.**(msmh_relation()(10**a))/(h*h))/((10.**a)/h) for a in x] #msun/h / msun/h
    #
    #sharc2=np.loadtxt('data/msmh_SHARC_2.dat')
    #sharc=np.loadtxt('data/msmh_SHARC.dat')
    #sharc3=np.loadtxt('data/msmh_SHARC_3.dat')
    #
    #pl.clf()
    #pl.plot(np.log10((10**x)/h),y2) #msun / ms/mh
    ##pl.plot(sharc[:,0],(10**sharc[:,1])/(10**sharc[:,0]))
    #pl.plot(sharc[:,0],(10**sharc[:,1])/(10**sharc[:,0])) #msun / mh/mh
    #pl.yscale('log')

data_dir = '/home/christoph/Research/age_matching_comparison/paper1/plotgen/data/'

#data_file_base = 'bp_z%s_centrals_hmbins_%d_%s_%d_med.bin'
#data_file_base = 'bp_centrals_stacked2_mmstar_%d_%s_med_%d.bin'
#data_file_base = 'bp_z0_centrals_cut_%s_pbins_%d_hist1D.bin'
data_file_base = 'bp_z0_centrals_hmbins_%d_%s_rho_%d_hist1D.dat'

data_file = data_dir + data_file_base

#rescale_xaxis(data_file,'2',percentiles=True)
#rescale_xaxis_stacked(data_file, percentiles=True)
#reconstructHist1D('data/bp_z0_centrals_hmbins_%d_bsr_almm_%d_hist1D.bin')
#genAccretionHistPlot('bp_z0_centrals_m12_rho_%d.gnu')
#medianAccretionHist('../../tree_code_outputs/bp_z0_centrals_m12_rho_%d_med.gnu')
#calcDerivativeAccretion('../../tree_code_outputs/bp_z0_centrals_m12_rho_%d_med.gnu')
#smooth_data('data/bp_z0_centrals_hmbins_rw%d_%s_p_s_%d_med.bin')
#smooth_data('data/bp_z0_centrals_hmbins_rw%d_%s_%d_med.bin')
#smooth_data('data/bp_z0_centrals_bsr_%s_hmbins_%d_med.bin')
#smooth_data('data/bp_z0_centrals_hmbins_%s_bsr_rho_%d_%s_med.dat')
#smooth_data('data/bp_z0_centrals_hmbins_%s_bsr_p_%s_%d_med.bin')
#normalize_prog_hist('data/bp_z0_centrals_hmbins_%s_rho_%d_mtree_med.dat2')
#normalize_prog_hist('data/bp_z0_centrals_hmbins_%s_bsr_%d_mtree_med.dat2')

#separate_mtrees('data/bp_z0_centrals_hmbins_0_bsr_0_med_mtree.dat')
#separate_mtrees('data/bp_z0_centrals_hmbins_0_bsr_1_med_mtree.dat')
#separate_mtrees('data/bp_z0_centrals_hmbins_0_bsr_2_med_mtree.dat')
#############################

#smooth_prog_hist("")
#normalize_prog_hist("")
smooth_data("")

