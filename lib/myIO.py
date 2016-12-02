import numpy as np

def wrtContour(path,X,Y,Z,hdr=['X','Y','Z']):
    fh=open(path,'w')
    shdr=hdr[0]+'\t'+hdr[1]+'t'+hdr[2]+'\n'
    for j in range(X.shape[1]):
        for i in range(X.shape[0]):
            line=('%.8e' % X[i,j])+'\t'+('%.8e' % Y[i,j])+'\t'+('%.8e' % Z[i,j])+'\n'
            fh.write(line)
    fh.close()
    