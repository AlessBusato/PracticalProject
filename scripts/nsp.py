import numpy as np
from numpy import linalg as LA
from scipy.stats import pearsonr
from scipy.spatial import distance
from matplotlib.lines import Line2D
NB_ROI = 360

def prediction_model(sc_mat, coupling):
    '''Prediction of stable FC based on SC using Gaussian linear diffusion model'''
    Q = np.linalg.inv(np.identity(NB_ROI) + coupling * sc_mat)
    cov = Q @ Q.T
    predicted_fc = np.zeros((NB_ROI,NB_ROI))
    for i in range(0,NB_ROI):
        for j in range(0,NB_ROI):
            predicted_fc[i,j] = cov[i,j]/(cov[i,i]*cov[j,j])**0.5
    return predicted_fc

def hierarchichal_clustering(fc_mat):
    '''Compute spectral clusters and return the number of clusters and their size at every level'''
    fc_mat[fc_mat<0]=0
    _, modes = LA.eigh(fc_mat)
    modes = np.fliplr(modes).T

    clus_num = [] #Number of clusters for each order
    clus_size = [] #Elements are arrays of cluster size for each order
    clusters = [] #Elements are arrays of positions for each clusters in the previous order
    clus_num.append(1) #First mode only has one module
    clus_size.append([360]) #First mode has one module of size 360

    clusters.append((modes[1]<0).nonzero()) #Initialize clusters of 2nd mode
    clusters.append((modes[1]>=0).nonzero())
    clus_num.append(2) #The 2nd mode always has 2 clusters (1 pos and 1neg)

    for order in range(2,NB_ROI):
        x = (modes[order]>=0).nonzero() #Get indices for positive values
        y = (modes[order]<0).nonzero() #Get indices for negative values
        length = []
        new_clus = []
        for cluster in clusters: #Find clusters of mode N
            length.append(np.size(cluster)) #Find length of clusters at order-1
            if np.size(cluster)<=1:
                new_clus.append(cluster)
            else:
                if np.size(np.intersect1d(x,cluster)) >= 1 : new_clus.append(np.intersect1d(x,cluster))
                if np.size(np.intersect1d(y,cluster)) >= 1 : new_clus.append(np.intersect1d(y,cluster))
        clus_size.append(length)
        clusters = new_clus #Update new cluster
        clus_num.append(len(clusters))
    #Find length of clusters at order=NB_ROI
    length=[]
    for cluster in clusters:
        length.append(np.size(cluster))
    clus_size.append(length)
    return clus_size, clus_num

def segint_component(fc, clus_size, clus_num):
    '''Given the spectral cluster, compute the integration and segregation components'''
    clus_num = np.divide(clus_num,360)
    fc = (fc + fc.T)/2
    fc[fc<0] = 0
    evals, modes = LA.eigh(fc)
    evals[evals<0] = 0
    evals = np.flip(np.square(evals))
    p=np.zeros(NB_ROI)
    hf=np.zeros(NB_ROI)    #Compute correction p
    for i in range(0,NB_ROI):
        c = 1/clus_num[i]
        p[i] = np.sum(abs(np.asarray(clus_size[i])-c)/NB_ROI)
    for i in range(0,NB_ROI):
        hf[i] = evals[i] * clus_num[i] * (1-p[i])
    hin = hf[0]/NB_ROI
    hse = np.sum(hf[1:])/NB_ROI
    return hin, hse

def gendata(threshold=0.1):
    '''Compute stable FC for each individual and then compute spectral partitions and int/seg components'''
    emp_stable = np.genfromtxt('../Data/WR/stable_FC.csv', delimiter=',')
    l_hin = []
    l_hse = []
    dist_to_tfc = []
    true_fc = []
    mean_corr = []
    '''Load True FC'''
    for filename in sorted(os.listdir('../Data/WR/FC/'), key=lambda f:int(''.join(filter(str.isdigit, f)))):
        fc = np.genfromtxt('../Data/WR/FC/'+filename, delimiter=',')
        fc[fc<0] = 0
        true_fc.append(fc)
    '''Compute the integration and segregation components for every SC'''
    for i,filename in enumerate(sorted(os.listdir('../Data/WR/SC_match/'), key=lambda f:int(''.join(filter(str.isdigit, f))))):
        print(filename)
        print(i)
        tfc = true_fc[i] #True FC associated to the SC (this is ensured by the sorted() in the for loop)
        indiv_hin = []
        indiv_hse = []
        ldist = []
        mean = []
        net = Connectome(atlas='../Data/MMP_atlas.txt', edges='../Data/WR/SC_match/'+filename)
        net.density_threshold(threshold)
        mat = net.get_normalizedlaplacian()
        for coupling in [13]: #range(1,31)
            fc = prediction_model(mat, coupling*5) #Compute stable fc
            print(coupling*5)
            clus_size, clus_num = hierarchichal_clustering(fc) #Compute spectral partitions
            hin, hse = segint_component(fc, clus_size, clus_num) #Compute seg/int components
            fc[fc<0] = 0
            dist = distance.euclidean(fc.flatten(), emp_stable.flatten()) #Compute distance between empirical FC(tfc) and simulated fc(fc) for c=coupling*5
            ldist.append(dist)
            mean.append(np.mean(fc))
            indiv_hin.append(hin)
            indiv_hse.append(hse)
        l_hin.append(indiv_hin)
        l_hse.append(indiv_hse)
        print(ldist)
        dist_to_tfc.append(ldist)
        mean_corr.append(mean)
    np.savetxt('../Data/Hierarchical_SegInt/thresholded/meancorr_'+str(threshold*100)+'.csv', np.asarray(mean_corr), delimiter=',')
    np.savetxt('../Data/Hierarchical_SegInt/thresholded/Hin_'+str(threshold*100)+'.csv', np.asarray(l_hin), delimiter=',')
    np.savetxt('../Data/Hierarchical_SegInt/thresholded/Hse_'+str(threshold*100)+'.csv', np.asarray(l_hse), delimiter=',')
    np.savetxt('../Data/Hierarchical_SegInt/thresholded/dist_'+str(threshold*100)+'.csv', np.asarray(dist_to_tfc), delimiter=',')


def plot_balance_thresholds():
    '''Plot the Hse/Hin balance across thresholds for c=65'''
    plt.rcParams['font.size'] = 34
    plt.rcParams['figure.figsize'] = 15, 10
    hin = np.genfromtxt('../Data/Hierarchical_SegInt/thresholded/Hin_all.csv', delimiter=',')
    hse = np.genfromtxt('../Data/Hierarchical_SegInt/thresholded/Hse_all.csv', delimiter=',')
    bal = hin-hse
    hse_sd = np.std(hse, axis=0)
    bal_av = np.mean(bal, axis=0)
    bal_sd = np.std(bal, axis=0)
    hin_av = np.mean(hin, axis=0)
    hse_av = np.mean(hse, axis=0)
    hin_sd = np.std(hin, axis=0)
    hse_sd = np.std(hse, axis=0)
    print(bal_sd)
    plt.rcParams['font.size'] = 34
    plt.rcParams['figure.figsize'] = 15, 10
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.arange(5,100,5), bal_av, linewidth=3, color='blue')
    ax.fill_between(np.arange(5,100,5), (bal_av-bal_sd), (bal_av+bal_sd), color='blue', alpha=.2)
    ax.set_xlabel('Network density (%)')
    ax.set_ylabel('Hb')
    plt.axvline(50, color='black', linestyle='--', linewidth=3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.show()


def main():
    # t = np.arange(0.05,1,0.05)
    # print(t)
    # for thresh in t:
    #     gendata(thresh)
    # plot_balance_thresholds()
    '''*********************************Plotting**************************************************'''
    mean = np.genfromtxt('../Data/Hierarchical_SegInt/meancorr_full4sess.csv', delimiter=',')
    mean_top = np.genfromtxt('../Data/Hierarchical_SegInt/meancorr_full4sess_th.csv', delimiter=',')
    true_fc = np.genfromtxt('../Data/WR/stable_FC.csv', delimiter=',')
    distance = np.genfromtxt('../Data/Hierarchical_SegInt/dist_fullpos.csv', delimiter=',')
    dist = np.mean(distance, axis=0)
    dist = dist/(np.amax(dist))
    dist_sd = np.std(dist)
    opt_c = (np.argmin(dist)+1)*5 #optimal coupling based on dist
    true_mean = np.mean(true_fc)
    sim_mean = np.mean(mean, axis=0)
    sim_sd = np.std(mean, axis=0)
    thresh_mean = np.mean(mean_top, axis=0)
    thresh_sd = np.std(mean_top, axis=0)
    plt.rcParams['font.size'] = 30
    plt.rcParams['figure.figsize'] = 15, 10
    hinf = np.genfromtxt('../Data/Hierarchical_SegInt/Hin_full4sess.csv', delimiter=',')
    hsef = np.genfromtxt('../Data/Hierarchical_SegInt/Hse_full4sess.csv', delimiter=',')
    hin = np.genfromtxt('../Data/Hierarchical_SegInt/Hin_full4sess_th.csv', delimiter=',')
    hse = np.genfromtxt('../Data/Hierarchical_SegInt/Hse_full4sess_th.csv', delimiter=',')
    hinf_av = np.mean(hinf, axis=0)
    hsef_av = np.mean(hsef, axis=0)
    hin_av = np.mean(hin, axis=0)
    hse_av = np.mean(hse, axis=0)
    hinf_sd = np.std(hinf, axis=0)
    hsef_sd = np.std(hsef, axis=0)
    hin_sd = np.std(hin, axis=0)
    hse_sd = np.std(hse, axis=0)
    print(hse_sd)
    hin_d = hinf_av-hin_av #full network - thresholded network
    hse_d = hsef_av-hse_av
    print(hinf_av)
    print(hsef_av)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    '''********************************Distance******************************************'''
    # ax.plot(np.arange(5,155,5), dist, linewidth=3, color='yellowgreen', label='distance')
    # ax.fill_between(np.arange(5,155,5), (dist-dist_sd), (dist+dist_sd), color='yellowgreen', alpha=.2)
    '''********************************Mean******************************************'''
    # ax.plot(np.arange(5,155,5), sim_mean, linewidth=3, color='blue', label='Full network average weight')
    # ax.fill_between(np.arange(5,155,5), (sim_mean-sim_sd), (sim_mean+sim_sd), color='blue', alpha=.2)
    # ax.plot(np.arange(5,155,5), thresh_mean, linewidth=3, color='red', label='10% network average weight')
    # ax.fill_between(np.arange(5,155,5), (thresh_mean-thresh_sd), (thresh_mean+thresh_sd), color='red', alpha=.2)
    # plt.axhline(true_mean, color='yellowgreen', linestyle='--', linewidth=3)
    '''********************************Integration/Segregation******************************************'''
    ax.plot(np.arange(5,155,5), hin_d, linewidth=3, color='blue', label='Integration')
    ax.plot(np.arange(5,155,5), hse_d, linewidth=3, color='red', label='Segregation')
    # ax.plot(np.arange(5,155,5), hinf_av, linewidth=3, color='blue', label='Full network integration')
    # ax.plot(np.arange(5,155,5), hin_av, linewidth=3, color='blue', linestyle='dashed', label='10% network integration')
    # ax.plot(np.arange(5,155,5), hsef_av, linewidth=3, color='red', label='Full network segregation')
    # ax.plot(np.arange(5,155,5), hse_av, linewidth=3, color='red', linestyle='dashed', label='10% network segregation')
    # ax.fill_between(np.arange(5,155,5), (hin_av-hin_sd), (hin_av+hin_sd), color='blue', alpha=.2)
    # ax.fill_between(np.arange(5,155,5), (hse_av-hse_sd), (hse_av+hse_sd), color='red', alpha=.2)
    # ax.fill_between(np.arange(5,155,5), (hinf_av-hinf_sd), (hinf_av+hinf_sd), color='blue', alpha=.2)
    # ax.fill_between(np.arange(5,155,5), (hsef_av-hsef_sd), (hsef_av+hsef_sd), color='red', alpha=.2)
    plt.axhline(0, color='black', linestyle='--', linewidth=3)

    plt.axvline(opt_c, color='yellowgreen', linestyle='--', linewidth=3)
    # plt.legend()
    ax.set_xlabel('Coupling')
    ax.set_ylabel(r'$\Delta$H')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.show()


if __name__== "__main__":
    main()
