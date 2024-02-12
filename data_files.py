from os import walk
import numpy as np

loc = '/Users/jonhouse/Desktop/research/soft_fbm'
class DataFile: 

    def __init__(self, path, type='base'):
        self.path = path
        self.weight, self.gamma, self.length, self.nconf, self.nt, self.nbin = 0,0,0,0,0,0
        self.series = ""
        self.avx, self.dis, self.log, self.cor = None, None, None, None
        if(type == 'base'):
            self.avx, self.dis, self.log,self.cor = self.get_data(self.path)
        if(type == 'grad'):
            self.avx, self.dis, self.log = self.get_grad_data(self.path)
        if(type == 'fast hosking'):
            self.avx, self.dis, self.log = self.get_fast_hosking_data(self.path)
        if(type == 'avx'):
            self.avx = self.get_avx(self.path)
        if(type == 'cor'):
            self.cor = self.get_cor(self.path)

        if(isinstance(type,list)):
            self.markers = self.parse_gen_data(type[0])



        
    def parse_gen_data(self,signal):
        run_read = open(loc + '/' + self.path,'r').readlines()
        N = len(run_read)
        ret_array = np.array([])
        
        for i in range(N):
            line = run_read[i].split()
            if (len(line) != 0):
                if line[0] == signal:
                    ret_array = np.append(ret_array,float(line[1]))
        
        return ret_array


    def parse_data(self, file_name,log_type):

        log_info = {'avx':{'signal':'time'},'dis':{'signal':'ibin'},'log':{'signal':')=P(ln(x))/x'},'cor':{'signal':'time'}}
        run_read = open(file_name,'r').readlines()
        N = len(run_read)
        offset = 0

        for i in range(N):
            line = run_read[i].split()
            if log_type == 'avx' and line[0] == 'weight':
                self.weight = float(line[3])

            if log_type == 'avx' and line[0] == 'GAMMMA=':
                self.gamma = float(line[1])

            if log_type == 'avx' and line[0] == 'L=':
                self.length = float(line[1])

            if log_type == 'avx' and line[0] == 'NCONF=':
                self.nconf = int(line[1])

            if log_type == 'avx' and line[0] == "NT=":
                self.nt = int(line[1])

            if line[0] == log_info[log_type]['signal']:
                offset = i + 1
                break

        N = N - offset

        if log_type == 'avx':
            run_data = np.empty((3,N),float)
        if log_type == 'dis':
            run_data = np.empty((7,N),float)
        if log_type == 'log':
            run_data = np.empty((7,N),float)
        if log_type == 'cor':
            run_data = np.empty((8,N),float)

        for i in range(N):
            run_data[:,i] = run_read[i+offset].split()

        data_dict = {}
        if log_type == 'avx':
            data_dict['t'] =        run_data[0,:]
            data_dict['<r>'] =      run_data[1,:]
            data_dict['<r^2>'] =    run_data[2,:]

        if log_type == 'dis':
            data_dict['ibin'] =     run_data[0,:]
            data_dict['x'] =        run_data[1,:]
            data_dict['x/L'] =      run_data[2,:]
            data_dict['P(x)'] =     run_data[3,:]
            data_dict['P(x)*L'] =   run_data[4,:]
            data_dict['P(|x|)'] =   run_data[5,:]
            data_dict['L/2-x'] =    run_data[6,:]

            self.nbin = data_dict['ibin'][0]*-1 # half bins are on left side, and indexing starts negative

        if log_type == 'log':
            data_dict['ilogbin'] =  run_data[0,:]
            data_dict['ln(x)'] =    run_data[1,:]
            data_dict['x'] =        run_data[2,:]
            data_dict['P(ln(x))'] = run_data[3,:]
            data_dict['P(x)'] =     run_data[4,:]
            data_dict['x/L'] =      run_data[5,:]
            data_dict['P(x)*L'] =   run_data[6,:]

        if log_type == 'cor': 
            data_dict['t'] = run_data[0,:]
            data_dict['<r^2>'] = run_data[1,:]
            data_dict['<xix_pos^2>'] = run_data[2,:]
            data_dict['<f_grad_pos^2>'] = run_data[3,:]
            data_dict['<xix_pos*f_grad_pos>'] = run_data[4,:]
            data_dict['<|xix_pos|>'] = run_data[5,:]
            data_dict['<|grad_pos|>'] = run_data[6,:]
            data_dict['norm(<xix_pos*f_grad_pos>)'] = run_data[7,:]

        return data_dict

    def get_data(self, folder):
        path = loc + '/' + folder
        avx, dis, log, cor = None, None, None, None
        for dirpaths,dirnames,files in walk(path):
            for f in files:
                type = f[0:3] # first three characters of the data log 
                if type == 'avx':
                    avx = self.parse_data(path + '/' + f,'avx')
                elif type == 'dis':
                    dis = self.parse_data(path + '/' + f, 'dis')
                elif type == 'log':
                    log = self.parse_data(path + '/' + f, 'log')
                elif type == 'cor':
                    cor = self.parse_data(path + '/' + f, 'cor')

        return avx, dis, log, cor


    def get_avx(self,folder):
        path = loc + '/' + folder
        avx = None
        for dirpaths,dirnames,files in walk(path):
            for f in files:
                type = f[0:3] # first three characters of the data log 
                if type == 'avx':
                    avx = self.parse_data(path + '/' + f,'avx')

        return avx
    
    def get_cor(self,folder):
        path = loc + '/' + folder
        cor = None
        for dirpaths,dirnames,files in walk(path):
            for f in files:
                type = f[0:3] # first three characters of the data log 
                if type == 'cor':
                    avx = self.parse_data(path + '/' + f,'cor')

        return cor
    


    def get_label(self, label):
            if label == 'gamma':
                return "$\gamma$=" + str(self.gamma)
            if label == 'nconf':
                return "nconf=" + str(self.nconf)
            if label == 'length':
                return "L=" + str(self.length)
            if label == 'weight': 
                return "weight=" + str(self.weight)
            if label == 'nt':
                return "nt=" + str(self.nt)
            if label == 'nbin':
                return "nbin=" + str(self.nbin)
            if label == "series":
                return "\n" + str(self.series)
            if label == 't':
                # we can use the path name as this specifies the time interval we measure our distribution on 
                intv_start = self.path.find("[") 
                intv_end = self.path[intv_start:].find("]")
                return "t=" + str(self.path[intv_start:intv_start+intv_end+1])

    def parse_grad_data(self, file_name,log_type):

        log_info = {'avx':{'signal':'time'},'dis':{'signal':'ibin'},'log':{'signal':'ilogbin'}}
        run_read = open(file_name,'r').readlines()
        N = len(run_read)
        offset = 0

        for i in range(N):
            line = run_read[i].split()
            if log_type == 'avx' and line[0] == 'weight':
                self.weight = float(line[3])

            if log_type == 'avx' and line[0] == 'GAMMMA=':
                self.gamma = float(line[1])

            if log_type == 'avx' and line[0] == 'L=':
                self.length = float(line[1])

            if log_type == 'avx' and line[0] == 'NCONF=':
                self.nconf = int(line[1])

            if log_type == 'avx' and line[0] == "NT=":
                self.nt = int(line[1])

            if line[0] == log_info[log_type]['signal']:
                offset = i + 1
                break

        N = N - offset

        if log_type == 'avx':
            run_data = np.empty((4,N),float)
        if log_type == 'dis':
            run_data = np.empty((7,N),float)
        if log_type == 'log':
            run_data = np.empty((7,N),float)

        for i in range(N):
            run_data[:,i] = run_read[i+offset].split()

        data_dict = {}
        if log_type == 'avx':
            data_dict['t'] =        run_data[0,:]
            data_dict['<nwalkers>'] =      run_data[1,:]
            data_dict['<r>'] =    run_data[2,:]
            data_dict['<r^2>'] =     run_data[3,:]

        if log_type == 'dis':
            data_dict['ibin'] =     run_data[0,:]
            data_dict['x'] =        run_data[1,:]
            data_dict['x/L'] =      run_data[2,:]
            data_dict['P(x)'] =     run_data[3,:]
            data_dict['P(x)*L'] =   run_data[4,:]
            data_dict['P(|x|)'] =   run_data[5,:]
            data_dict['L/2-x'] =    run_data[6,:]

            self.nbin = data_dict['ibin'][0]*-2 # half bins are on left side, and indexing starts negative

        if log_type == 'log':
            data_dict['ilogbin'] =  run_data[0,:]
            data_dict['ln(x)'] =    run_data[1,:]
            data_dict['x'] =        run_data[2,:]
            data_dict['P(ln(x))'] = run_data[3,:]
            data_dict['P(x)'] =     run_data[4,:]
            data_dict['x/L'] =      run_data[5,:]
            data_dict['P(x)*L'] =   run_data[6,:]

        return data_dict

    def get_grad_data(self, folder):
        path = loc + '/' + folder
        avx, dis, log = None, None, None
        for dirpaths,dirnames,files in walk(path):
            for f in files:
                type = f[0:3] # first three characters of the data log 
                if type == 'avx':
                    avx = self.parse_grad_data(path + '/' + f,'avx')
                elif type == 'dis':
                    dis = self.parse_grad_data(path + '/' + f, 'dis')
                elif type == 'log':
                    log = self.parse_grad_data(path + '/' + f, 'log')


        return avx, dis, log
    

    def parse_fast_hosking_data(self, file_name,log_type):

        log_info = {'avx':{'signal':'time'},'dis':{'signal':'ibin'},'log':{'signal':'ilogbin'}}
        run_read = open(file_name,'r').readlines()
        N = len(run_read)
        offset = 0

        for i in range(N):
            line = run_read[i].split()

            if log_type == 'avx' and line[0] == 'GAMMMA=':
                self.gamma = float(line[1])

            if log_type == 'avx' and line[0] == 'L=':
                self.length = float(line[1])

            if log_type == 'avx' and line[0] == 'NCONF=':
                self.nconf = int(line[1])

            if log_type == 'avx' and line[0] == "NT=":
                self.nt = int(line[1])

            if line[0] == log_info[log_type]['signal']:
                offset = i + 1
                break

        N = N - offset

        if log_type == 'avx':
            run_data = np.empty((8,N),float)
        if log_type == 'dis':
            run_data = np.empty((7,N),float)
        if log_type == 'log':
            run_data = np.empty((7,N),float)

        for i in range(N):
            run_data[:,i] = run_read[i+offset].split()

        data_dict = {}
        if log_type == 'avx':
            data_dict['t'] =                run_data[0,:]
            data_dict['configs'] =          run_data[1,:]
            data_dict['<nwalkers>'] =       run_data[2,:]
            data_dict['<std.dev.nw>'] =     run_data[3,:]
            data_dict['<r>'] =              run_data[4,:]
            data_dict['<r^2>'] =            run_data[5,:]
            data_dict['std.dev.r'] =        run_data[6,:]
            data_dict['std.dev.r^2'] =      run_data[7,:]

        if log_type == 'dis':
            data_dict['ibin'] =     run_data[0,:]
            data_dict['x'] =        run_data[1,:]
            data_dict['x/L'] =      run_data[2,:]
            data_dict['P(x)'] =     run_data[3,:]
            data_dict['P(x)*L'] =   run_data[4,:]
            data_dict['P(|x|)'] =   run_data[5,:]
            data_dict['L/2-x'] =    run_data[6,:]

            self.nbin = data_dict['ibin'][0]*-2 # half bins are on left side, and indexing starts negative

        if log_type == 'log':
            data_dict['ilogbin'] =  run_data[0,:]
            data_dict['ln(x)'] =    run_data[1,:]
            data_dict['x'] =        run_data[2,:]
            data_dict['P(ln(x))'] = run_data[3,:]
            data_dict['P(x)'] =     run_data[4,:]
            data_dict['x/L'] =      run_data[5,:]
            data_dict['P(x)*L'] =   run_data[6,:]

        return data_dict

    def get_fast_hosking_data(self, folder):
        path = loc + '/' + folder
        avx, dis, log = None, None, None
        for dirpaths,dirnames,files in walk(path):
            for f in files:
                type = f[0:3] # first three characters of the data log 

                if type == 'avx':
                    avx = self.parse_fast_hosking_data(path + '/' + f,'avx')
                elif type == 'dis':
                    dis = self.parse_fast_hosking_data(path + '/' + f, 'dis')
                elif type == 'log':
                    log = self.parse_fast_hosking_data(path + '/' + f, 'log')


        return avx, dis, log