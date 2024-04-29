import numpy as np
from os import walk

loc = '/Users/jonhouse/Desktop/research/soft_fbm'

def parse_data(file_name,log_type):

    log_info = {'avx':{'signal':'time'},'dis':{'signal':'ibin'},'log':{'signal':')=P(ln(x))/x'}}
    run_read = open(file_name,'r').readlines()
    N = len(run_read)
    offset = 0

    for i in range(N):
        if run_read[i].split()[0] == log_info[log_type]['signal']:
            offset = i + 1
            break

    N = N - offset

    if log_type == 'avx':
        run_data = np.empty((3,N),float)
    if log_type == 'dis':
        run_data = np.empty((7,N),float)
    if log_type == 'log':
        run_data = np.empty((7,N),float)

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

    if log_type == 'log':
        data_dict['ilogbin'] =  run_data[0,:]
        data_dict['ln(x)'] =    run_data[1,:]
        data_dict['x'] =        run_data[2,:]
        data_dict['P(ln(x))'] = run_data[3,:]
        data_dict['P(x)'] =     run_data[4,:]
        data_dict['x/L'] =      run_data[5,:]
        data_dict['P(x)*L'] =   run_data[6,:]

    return data_dict

def get_data(folder):
    path = loc + '/' + folder
    avx, dis = None, None
    for dirpaths,dirnames,files in walk(path):
        for f in files:
            type = f[0:3] # first three characters of the data log 
            if type == 'avx':
                avx = parse_data(path + '/' + f,'avx')
            elif type == 'dis':
                dis = parse_data(path + '/' + f, 'dis')

    return avx, dis

def parse_tempered_data(file_name,log_type):

    log_info = {'avx':{'signal':'time'},'dis':{'signal':'ibin'},'log':{'signal':')=P(ln(x))/x'}}
    run_read = open(file_name,'r').readlines()
    N = len(run_read)
    offset = 0

    for i in range(N):
        if run_read[i].split()[0] == log_info[log_type]['signal']:
            offset = i + 1
            break

    N = N - offset

    if log_type == 'avx':
        run_data = np.empty((3,N),float)
    if log_type == 'dis':
        run_data = np.empty((7,N),float)
    if log_type == 'log':
        run_data = np.empty((7,N),float)

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

    if log_type == 'log':
        data_dict['ilogbin'] =  run_data[0,:]
        data_dict['ln(x)'] =    run_data[1,:]
        data_dict['x'] =        run_data[2,:]
        data_dict['P(ln(x))'] = run_data[3,:]
        data_dict['P(x)'] =     run_data[4,:]
        data_dict['x/L'] =      run_data[5,:]
        data_dict['P(x)*L'] =   run_data[6,:]

    return data_dict

def get_tempered_data(folder):
    path = loc + '/' + folder
    for dirpaths,dirnames,files in walk(path):
        for f in files:
            type = f[0:3] # first three characters of the data log 
            if type == 'avx':
                avx = parse_tempered_data(path + '/' + f,'avx')
            elif type == 'dis':
                dis = parse_tempered_data(path + '/' + f, 'dis')
            elif type == 'log':
                log = parse_tempered_data(path + '/' + f, 'log')

    return avx, dis, log

def get_err_data(folder):
    path = loc + '/' + folder
    for dirpaths,dirnames,files in walk(path):
        for f in files:
            type = f[0:3] # first three characters of the data log 
            if type == 'avx':
                avx = parse_err_data(path + '/' + f,'avx')
            elif type == 'dis':
                dis = parse_err_data(path + '/' + f, 'dis')
            elif type == 'log':
                log = parse_err_data(path + '/' + f, 'log')

    return avx, dis, log   


def parse_err_data(file_name,log_type):

    log_info = {'avx':{'signal':'time'},'dis':{'signal':'ibin'},'log':{'signal':'n(x))/x'}}
    run_read = open(file_name,'r').readlines()
    N = len(run_read)
    offset = 0

    for i in range(N):
        if run_read[i].split()[0] == log_info[log_type]['signal']:
            offset = i + 1
            break

    N = N - offset

    if log_type == 'avx':
        run_data = np.empty((7,N),float)
    if log_type == 'dis':
        run_data = np.empty((9,N),float)
    if log_type == 'log':
        run_data = np.empty((9,N),float)

    for i in range(N):
        run_data[:,i] = run_read[i+offset].split()

    data_dict = {}
    if log_type == 'avx':
        data_dict['t'] =            run_data[0,:]
        data_dict['<nwalker>'] =    run_data[1,:]
        data_dict['<r>'] =          run_data[2,:]
        data_dict['<r^2>'] =        run_data[3,:]
        data_dict['<r^4>'] =        run_data[4,:]
        data_dict['var'] =          run_data[5,:]
        data_dict['err'] =          run_data[6,:]

    if log_type == 'dis':
        data_dict['ibin'] =     run_data[0,:]
        data_dict['x'] =        run_data[1,:]
        data_dict['x/L'] =      run_data[2,:]
        data_dict['P(x)'] =     run_data[3,:]
        data_dict['P(x)*L'] =   run_data[4,:]
        data_dict['P(|x|)'] =   run_data[5,:]
        data_dict['L/2-x'] =    run_data[6,:]
        data_dict['var'] =      run_data[7,:]
        data_dict['err'] =      run_data[8,:]

    if log_type == 'log':
        data_dict['ilogbin'] =  run_data[0,:]
        data_dict['ln(x)'] =    run_data[1,:]
        data_dict['x'] =        run_data[2,:]
        data_dict['P(ln(x))'] = run_data[3,:]
        data_dict['P(x)'] =     run_data[4,:]
        data_dict['x/L'] =      run_data[5,:]
        data_dict['P(x)*L'] =   run_data[6,:]
        data_dict['var'] =      run_data[7,:]
        data_dict['err'] =      run_data[8,:]

    return data_dict


def get_fast_data(folder):
    path = loc + '/' + folder
    for dirpaths,dirnames,files in walk(path):
        for f in files:
            type = f[0:3] # first three characters of the data log 
            if type == 'avx':
                avx = parse_fast_data(path + '/' + f,'avx')
            elif type == 'dis':
                dis = parse_fast_data(path + '/' + f, 'dis')
            elif type == 'log':
                log = parse_fast_data(path + '/' + f, 'log')

    return avx, dis, log   


def parse_fast_data(file_name,log_type):

    log_info = {'avx':{'signal':'time'},'dis':{'signal':'ibin'},'log':{'signal':')=P(ln(x))/x'}}
    run_read = open(file_name,'r').readlines()
    N = len(run_read)
    offset = 0

    for i in range(N):
        if run_read[i].split()[0] == log_info[log_type]['signal']:
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
        data_dict['t'] =            run_data[0,:]
        data_dict['configs'] =      run_data[1,:]
        data_dict['<nwalker>'] =    run_data[2,:]
        data_dict['std.dev.nw'] =   run_data[3,:]
        data_dict['<r>'] =          run_data[4,:]
        data_dict['<r^2>'] =        run_data[5,:]
        data_dict['std.dev.r'] =    run_data[6,:]
        data_dict['std.dev.r^2'] =  run_data[7,:]

    if log_type == 'dis':
        data_dict['ibin'] =     run_data[0,:]
        data_dict['x'] =        run_data[1,:]
        data_dict['x/L'] =      run_data[2,:]
        data_dict['P(x)'] =     run_data[3,:]
        data_dict['P(x)*L'] =   run_data[4,:]
        data_dict['P(|x|)'] =   run_data[5,:]
        data_dict['L/2-x'] =    run_data[6,:]

    if log_type == 'log':
        data_dict['ilogbin'] =  run_data[0,:]
        data_dict['ln(x)'] =    run_data[1,:]
        data_dict['x'] =        run_data[2,:]
        data_dict['P(ln(x))'] = run_data[3,:]
        data_dict['P(x)'] =     run_data[4,:]
        data_dict['x/L'] =      run_data[5,:]
        data_dict['P(x)*L'] =   run_data[6,:]


    return data_dict