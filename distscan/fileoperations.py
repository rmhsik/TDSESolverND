import numpy as np
def getlinesnumber(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    return len(lines)

def getdata(filepath,number):
    data = []
    linecounter = 0
    with open(filepath, "r") as f:
        for k in range(number):
            line = f.readline()
            if not line:
                break
            data.append([linecounter, np.array(line.strip('\n').split(' '), dtype=np.double)] )
            linecounter +=1

    if len(data)%number != 0 :
        offset = number - len(data)
        for i in range(offset):
            data.append([-1 ,np.array([np.nan,np.nan])] )
    return data

def removedataline(filepath, linenumber):
    with open(filepath,'r') as f:
        lines = f.readlines()

    with open(filepath,'w') as f:
        for number, line in enumerate(lines):
            if number != linenumber:
                f.write(line)

def LoadComplexData(file,**genfromtext_args):
    """
    Load complex data in the C++ format in numpy.
    """
    array_as_strings = np.genfromtxt(file,dtype=str,**genfromtext_args)
    complex_parser = np.vectorize(lambda x: complex(*eval(x)))
    return complex_parser(array_as_strings)
