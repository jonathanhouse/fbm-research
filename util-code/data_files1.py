from os import walk
import numpy as np

loc = '/Users/jonhouse/Desktop/research/soft_fbm'

class DataFile1: 

    def __init__(self,path,grab_avx=True,grab_dis=True,grab_log=False,grab_out=False,grab_full=True):
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
            self.dis = self.make_data_dict('dis')
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




