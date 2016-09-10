### Contains general tools for python###


def listAlloc(n,m=0):
    l=[]
    if (m<1):
        for i in range(n):
            l.append()
    else:
        for i in range(n+1):
            l.append([])
            for j in range(m+1):
                l[i].append(None)
            
    return l
