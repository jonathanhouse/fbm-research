from os import walk
import numpy as np
from plot import ordered_binning

loc = '/Users/jonhouse/Desktop/research/soft_fbm'

class DataFile: 

    def __init__(self,path,grab_avx=True,grab_dis=True,grab_log=False,grab_out=False,grab_full=True, dis_prefix='dis'):
        self.path = path
        self.avx = {}
        self.dis = {}
        self.log = {}
        self.full = {}
        self.gen_out = None # general purpose output 
        self.tags = {}
        self.params = {}

        if (grab_out):
            self.gen_out = self.make_data_dict('out')
            return

        if (grab_avx):
            self.avx = self.make_data_dict('avx')
        if (grab_dis):
            self.dis = self.make_data_dict(dis_prefix)
            if "P(x)" in list(self.dis.keys()):
                self.dis["P(|x|)"] = 0.5 * (self.dis["P(x)"] + self.dis["P(x)"][::-1]) # hard code symmetrization of P(x)
        if (grab_log):
            self.log = self.make_data_dict('log')
        if (grab_full):
            if (bool(self.avx)):
                self.full.update(self.avx)
            if (bool(self.dis)):
                self.full.update(self.dis)
            if (bool(self.log)):
                self.full.update(self.log)
        print("Finished loading " + path)

        

    def file_parser(self, filename):
        run_read = open(loc + '/' + self.path + filename,'r').readlines() # reads each line as a string 
        num_lines = len(run_read)

        if (filename == ""):
            return run_read

        type = filename[0:3] # some files might have longer identifiers so this is just arbitrary

        offset = 0
        for n in range(num_lines):
            line = run_read[n].split()
            for z in range(len(line)):
                if (line[z][-1] == '=' and (z != len(line) - 1)): # if the param ends in an '=" and isn't the long string of them
                    self.params[line[z][0:-1]] = line[z+1]
                if (line[z] == '=================================='):
                    offset = n 
                    break
        
        header = run_read[offset + 1].split()
        header_len = len(header)

        data_spaces = num_lines - (offset + 2)
        file_data = np.empty((header_len,data_spaces),float)

        self.tags[type] = header
        data_dict = {}

        for k in range(data_spaces):
            splits = run_read[k + offset + 2].split()
            for i in range(len(splits)):
                try:
                    file_data[i,k] = splits[i]
                except ValueError: 
                    file_data[i,k] = 0
            #file_data[:,k] = run_read[k + offset + 2].split()

        for t in range(len(header)):
            data_dict[header[t]] = file_data[t,:]

        return data_dict



    def make_data_dict(self,dtype):
        if (dtype == 'out'):
            return self.file_parser(str()) # don't need to find file because we just include it in path

        full_path = loc + '/' + self.path
        save_filename = ''
        for dirpaths,dirnames,files in walk(full_path):
            for f in files:
                type = f[0:len(dtype)] 
                if type == dtype:
                    save_filename = f

        if(save_filename == ''):
            return None
        else:
            return self.file_parser('/' + save_filename)
        
    def trim_dis_file(self,binsize,dis_type):
        prog_title = "Mean-density FBM Gaussian Binned (procs as sets)"
        desc = ''
        if(dis_type == 'pdis'):
            desc = 'Instantaneous probability density distribution'
        if(dis_type == 'idis'):
            desc = 'Integrated density distribution'



        low_bound = 0
        high_bound = 0
        col_len = 20
        ep = 1e-10
        high_bound = len(self.dis["ibin"]) - 1
        for l in range(len(self.dis["P(|x|)"])):
            if(self.dis["P(|x|)"][l] > ep):
                low_bound = l
                break
        for h in range(len(self.dis["P(|x|)"])):
            if(self.dis["P(|x|)"][len(self.dis["P(|x|)"]) - 1 - h] > ep):
                high_bound = len(self.dis["P(|x|)"]) - 1 - h
                break

        NT = int(self.params["NT"])
        filename = f'{dis_type}{NT}.dat'
        with open(filename, 'w') as file:
            # Write the header for clarity
            file.write("Program " + prog_title + '\n')
            file.write(desc + '\n')

            for key in list(self.params.keys()):
                file.write(f"{key}= {self.params[key]}\n")
            file.write(f"REBIN_BINWIDTH_FACTOR= {binsize}\n")

            last_avx_idx = len(self.avx["time"]) - 1
            last_avx_t = int(self.avx["time"][last_avx_idx])
            last_avx_msd = self.avx["<r^2>"][last_avx_idx]
            x_rms = np.sqrt(last_avx_msd)

            file.write(f"x_rms(t_last={last_avx_t})= {x_rms}\n")
            file.write("==================================\n")

            key_list = []
            if(dis_type == 'pdis'):
                key_list = ['ibin','x','x/x_rms','P(|x|)','x_rms*P']
            if(dis_type == 'idis'):
                key_list =  ['ibin','x','x/x_rms','P(|x|)','(x_rms/t_last)P']


            header=''
            for key in key_list:
                header += f" {key:<15}"
            file.write(header + '\n')

            y_prime, bin_prime, x_prime = ordered_binning(df=self,binsize=int(binsize), low_bound=low_bound, high_bound=high_bound)

            indices_prime = bin_prime + int(self.params["NBINS"]) # -NBINS + NBINS = 0, NBINS + NBINS=2*NBINS
            y_prime = y_prime / binsize # keep normalization

            data_out = []
            if(dis_type == 'pdis'):
                 data_out = [bin_prime,x_prime,x_prime/x_rms,y_prime,(x_rms)*y_prime]
            if(dis_type == 'idis'):
                 data_out = [bin_prime,x_prime,x_prime/x_rms,y_prime,(x_rms/last_avx_t)*y_prime]

            if(len(y_prime) != len(x_prime)):
                print("Rebinned to unmatching lengths")
                return

            for i,bin in enumerate(bin_prime):
                row = ''
                for d in data_out:
                    row += f" {d[i]:13.6E}  "
                file.write(row + '\n')


            # file.write(f"{'x':<col_len}{'P(|x|)'}\n")
            
        print(f"Output saved to {filename}")    

        




