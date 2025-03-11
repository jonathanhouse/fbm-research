import matplotlib.pyplot as plt
import numpy as np
from plot import plot, msd_fit, plot_times, plot_binned, log_fit, gen_fit, signed_square, signed_sqrt, log_fit2, unscaled_dis_plot, scaled_cdis_plot, scaled_idis_plot, unscaled_idis_plot, plot_msd_with_fits, cdis_multi_plot, idis_multi_plot, fbm_cdis_multi_plot,plot_msds_by_prefactor,fbm_idis_multi_plot,file_output_rebinned,fbm_integrated_density
from numpy.polynomial.polynomial import Polynomial as Poly
from decimal import Decimal
#import tikzplotlib as tz # see github repo to use tikz via venv
from data_files import DataFile

def suptitle_gen1(x,fig):
    title = "2*NBIN/L=" + "{:.2E}".format(Decimal(2*float(x.params["NBINS"]))) + "/" + "{:.2E}".format(Decimal(x.params["L"])) + "=" + str(2*float(x.params["NBINS"])/float(x.params["L"]))
    title += ", gamma=" + str(float(x.params["GAMMMA"])) 
    title += ", nconf=" + str(x.params["NCONF"]) 
    title += ", NT=$2^{" + str(np.log(float(x.params["NT"]))/np.log(2)) + "}$"

    fig.suptitle(title)


grab_cdis = False
grab_idis = False
grab_cdis_sup = False
grab_idis_fbm = False


BM_23 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.0/walks_per_set=64/nt=2^23",grab_dis=False)
BM_25 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.0/walks_per_set=64/nt=2^25",grab_dis=grab_cdis)
BM_27 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.0/walks_per_set=64/nt=2^27/cdis",grab_dis=grab_cdis)

BM_F26 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10 +  finite-space/gamma=1.0/L=1K/nt=2^20; sig=0.1",grab_dis=False)

idis_27 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.0/walks_per_set=64/nt=2^27/idis(dt=9K)",grab_dis=grab_idis)
idis_25 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.0/walks_per_set=64/nt=2^25/idis(dt=9K,nsets=480)",grab_dis=grab_idis)
idis_23 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.0/walks_per_set=64/nt=2^23/idis(dt=9K,nsets=480)",grab_dis=grab_idis)
idis_14 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.0/walks_per_set=64/nt=2^14/idis(dt=10,nsets=1080)",grab_dis=grab_idis)
idis_18 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.0/walks_per_set=64/nt=2^18/idis(dt=100,nsets=1080)",grab_dis=grab_idis)


cdis_14 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.0/walks_per_set=64/nt=2^14/cdis",grab_dis=grab_cdis)
# cdis_18 =DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.0/walks_per_set=64/nt=2^18/cdis",grab_dis=False)
cdis_23_480 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.0/walks_per_set=64/nt=2^23/cdis(nsets=480)",grab_dis=grab_cdis)
cdis_20 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.0/walks_per_set=64/nt=2^20",grab_dis=grab_cdis)

cdis_sup_14 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=0.5/walks_per_set=64/nt=2^14", grab_dis=grab_cdis_sup)
cdis_sup_18 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=0.5/walks_per_set=64/nt=2^18", grab_dis=grab_cdis_sup)
cdis_sup_23 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=0.5/walks_per_set=64/nt=2^23", grab_dis=grab_cdis_sup)
cdis_sup_25 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=0.5/walks_per_set=64/nt=2^25/L=3M", grab_dis=grab_cdis_sup)
cdis_sup_27 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=0.5/walks_per_set=64/nt=2^27/L=8M", grab_dis=grab_cdis_sup)
cdis_sup_20 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=0.5/walks_per_set=64/nt=2^20",grab_dis=grab_cdis_sup)

msd_g05 = DataFile("final data + figures/msd/gamma=0.5",grab_dis=False)
msd_g03 = DataFile("final data + figures/msd/gamma=0.3",grab_dis=False)
msd_g10 = DataFile("final data + figures/msd/gamma=1.0",grab_dis=False)
msd_g13 = DataFile("final data + figures/msd/gamma=1.3",grab_dis=False)
msd_g23 = DataFile("final data + figures/msd/gamma=2:3",grab_dis=False)

msd_pre_25 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.3",grab_dis=False)
msd_pre_16 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.3/prefactors/weight=1:16",grab_dis=False)
msd_pre_100 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.3/prefactors/weight=1.0",grab_dis=False)
msd_pre_400 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.3/prefactors/weight=4.0",grab_dis=False)
msd_pre_32 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=1.3/prefactors/weight=1:32",grab_dis=False)

fbm_idis_14 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=0.5/walks_per_set=64/nt=2^14/idis(dt=10)",grab_dis=grab_idis_fbm)
fbm_idis_18 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=0.5/walks_per_set=64/nt=2^18/idis(dt=100)",grab_dis=grab_idis_fbm)
fbm_idis_23 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=0.5/walks_per_set=64/nt=2^23/nset=240(dt=9K)",grab_dis=grab_idis_fbm)
fbm_idis_25 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=0.5/walks_per_set=64/nt=2^25/L=3M/idis(dt=9K)",grab_dis=grab_idis_fbm)
fbm_idis_27 = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=0.5/walks_per_set=64/nt=2^27/L=8M/dt=9K",grab_dis=grab_idis_fbm)


test_dis = DataFile("code/test_data",grab_dis=True)
theoretical_dis = DataFile("data/parallel walkers/procs as sets/gaussian-binned/fbm-on/bin:len=10/gamma=0.5/theoretical-fbm-cdis",grab_dis=True)

bm_dis = DataFile("code/test_data/normal-bm",grab_dis=True)

pure_g2 = DataFile("../bfbm/data/fall-23-runs/pure fbm/gamma=0.2",grab_dis=False)
pure_g8 = DataFile("../bfbm/data/22-23 data/cutoff_data/Gamma08/L10000",grab_dis=False)
pure_g10 = DataFile("../bfbm/data/22-23 data/cutoff_data/Gamma10/L1000",grab_dis=False)
pure_g12 = DataFile("../bfbm/data/22-23 data/cutoff_data/Gamma12/L1000/NT=2^22",grab_dis=False)
pure_g16 = DataFile("../bfbm/data/fall-23-runs/pure fbm/gamma=1.6",grab_dis=False)
pure_g15 = DataFile("../bfbm/data/22-23 data/cutoff_data/Gamma15/L300",grab_dis=False)
pure_g4 = DataFile("../bfbm/data/22-23 data/cutoff_data/Gamma04/L100000/NT=2^22",grab_dis=False)

idis_a5_14 = DataFile("code/figure_data/gamma=1.5/NT=2^14/pdis",grab_dis=True,dis_prefix='pdis')
idis_a5_18 = DataFile("code/figure_data/gamma=1.5/NT=2^18/pdis",grab_dis=True,dis_prefix='pdis')
idis_a5_23 = DataFile("code/figure_data/gamma=1.5/NT=2^23/pdis",grab_dis=True,dis_prefix='pdis')
idis_a5_25 = DataFile("code/figure_data/gamma=1.5/NT=2^25/pdis",grab_dis=True,dis_prefix='pdis')
idis_a5_27 = DataFile("code/figure_data/gamma=1.5/NT=2^27/pdis",grab_dis=True,dis_prefix='pdis')

NSETS = 16
NWALKS_PER_SET = 128
NCONF = NWALKS_PER_SET * NSETS
LEN_PER_BIN = 0.1

PARTIAL_SIM_NORM = NCONF*LEN_PER_BIN
NORM_CORRECTION = 16.0 / 128.0

### FINAL PREFACTOR-VARIED MSD
# colors = ['blue', 'red', 'green', 'orange', 'purple']
# markers = ['s', 'o', '^', 'D', 'x']
# msds = [msd_pre_400,msd_pre_100,msd_pre_25,msd_pre_16,msd_pre_32]
# plot_msds_by_prefactor(msds,colors,markers)


### FINAL ALPHA-VARIED MSD
# intv1 = [1,60412582]
# intv2 = [10679537, 60412582]
# intv3 = [1,5544]
# # intvs = [intv1, intv1, intv1, intv2, intv2]
# intvs=[intv3,intv3,intv3,intv3,intv3]
# colors = ['blue', 'red', 'green', 'orange', 'purple']
# markers = ['s', 'o', '^', 'D', 'x']
# # msds = [msd_g03, msd_g05, msd_g23, msd_g10, msd_g13]
# msds = [pure_g4,pure_g8, pure_g10,pure_g12,pure_g15]
# plot_msd_with_fits(msds, intvs, colors, markers, False)
# tz.save("pure_fbm_msd.tex")


### FINAL NOISE-DOMINATED INTEGRATED DISTRIBUTIONS 
# fbm_cdis_list = [cdis_sup_27, cdis_sup_25,cdis_sup_23]
# unscaled_fbm_cdis = [cdis_sup_27, cdis_sup_25,cdis_sup_23, cdis_sup_20, cdis_sup_14]
# fits = [ {"x" : theoretical_dis.dis["x_2"]/1e5, "P" : theoretical_dis.dis["P_2"]}, {"x" : theoretical_dis.dis["x_3"]/1e5, "P" : theoretical_dis.dis["P_3"]} ] 
# x_test1 = np.linspace(-40e5,40e5,250)
# x_test2 = np.linspace(-20e5,20e5,250)
# fits = [ {"x" : x_test1/1e5, "P" : fbm_integrated_density(x_test1,2**27,1.5)}, {"x" : x_test2/1e5, "P" : fbm_integrated_density(x_test2,2**25,1.5)} ] 
# fbm_cdis_multi_plot(fbm_cdis_list, 2**27, unscaled_dis_list=unscaled_fbm_cdis,fits=fits)
# plt.show()

### FINAL NOISE-DOMINATED INSTANTANEOUS DISTRIBUTIONS 
# fbm_idis_list = [fbm_idis_27,fbm_idis_25, fbm_idis_23]
# unscaled_fbm_idis = [fbm_idis_27, fbm_idis_25,fbm_idis_23, fbm_idis_18, fbm_idis_14]
# fbm_idis_multi_plot(fbm_idis_list, 2**27, unscaled_dis_list=unscaled_fbm_idis, to_fit_list=[fbm_idis_27,fbm_idis_25])

### FINAL FORCE-DOMINATED INTEGRATED DISTRIBUTIONS 
# dis_list = [BM_27, BM_25, cdis_23_480]
# unscaled_dis_list = [BM_27, BM_25, cdis_23_480, cdis_20, cdis_14]
# cdis_multi_plot(dis_list, 2**27, unscaled_dis_list=unscaled_dis_list, extra_plots=[bm_dis], to_fit_list=[BM_27, BM_25])
# plt.show()


### FINAL FORCE-DOMINATED INSTANTANEOUS DISTRIBUTIONS 
#idis_list = [idis_27, idis_25, idis_23, idis_18, idis_14]
#idis_multi_plot(idis_list, 2**27, unscaled_dis_list=idis_list)
#plt.show()

# tz.save("idis_twopanel.tex")









