from lib.grid import *

nbkx=3;nbky=2;blocks=[];xpos=-0.5

lxibk=([120,320,360])
letbk=([120,360])

corners=[]
lines=[]
for j in range(nbky):
    corners.append([])
    lines.append([])
    for i in range(nbkx):
        corners[j].append([])
        lines[j].append([])
        
bkx=0;bky=0

corners[bky][bkx].append(corner(-7.0,-7.0))
corners[bky][bkx].append(corner(-3.0,-7.0))
corners[bky][bkx].append(corner(-7.0,0.0))
corners[bky][bkx].append(corner(-0.5,0.0))

lines[bky][bkx].append(line(lxibk[bkx],corners[bky][bkx][0],corners[bky][bkx][1],0))
lines[bky][bkx].append(line(lxibk[bkx],corners[bky][bkx][2],corners[bky][bkx][3],0))
lines[bky][bkx].append(line(letbk[bky],corners[bky][bkx][0],corners[bky][bkx][2],1))
lines[bky][bkx].append(line(letbk[bky],corners[bky][bkx][1],corners[bky][bkx][3],1))

blocks.append(block(lines[bky][bkx]))

blocks[0].be.c1.dxdy=1.0
blocks[0].be.c0.dxdy=0.0
blocks[0].update()

bkx=1;bky=0

corners[bky][bkx].append(corners[bky][bkx-1][1].clone())
corners[bky][bkx].append(corner(0.5,-7.0))
corners[bky][bkx].append(corners[bky][bkx-1][3].clone())
corners[bky][bkx].append(corner(0.5,0.0))

lines[bky][bkx].append(line(lxibk[bkx],corners[bky][bkx][0],corners[bky][bkx][1],0))
lines[bky][bkx].append(line(lxibk[bkx],corners[bky][bkx][2],corners[bky][bkx][3],0,naca=-21.0))
lines[bky][bkx].append(lines[bky][bkx-1][3])
lines[bky][bkx].append(line(letbk[bky],corners[bky][bkx][1],corners[bky][bkx][3],1))

blocks.append(block(lines[bky][bkx]))
blocks[1].bw.c1.dydet=blocks[1].bn.c0.dxdxi*(np.sqrt(2.0)/2.0)
blocks[1].be.c1.dydet=blocks[1].bw.c1.dydet
blocks[1].update()
dxn=blocks[1].bn.c1.dxdxi
dxs=blocks[1].bs.c1.dxdxi
blocks[0].bn.c1.dydx=0.0
blocks[0].update()

bkx=2;bky=0

corners[bky][bkx].append(corners[bky][bkx-1][1])
corners[bky][bkx].append(corner(12.0,-7.0))
corners[bky][bkx].append(corners[bky][bkx-1][3])
corners[bky][bkx].append(corner(12,0.0))

lines[bky][bkx].append(line(lxibk[bkx],corners[bky][bkx][0],corners[bky][bkx][1],0))
lines[bky][bkx].append(line(lxibk[bkx],corners[bky][bkx][2],corners[bky][bkx][3],0))
lines[bky][bkx].append(lines[bky][bkx-1][3])
lines[bky][bkx].append(line(letbk[bky],corners[bky][bkx][1],corners[bky][bkx][3],1))

blocks.append(block(lines[bky][bkx]))
blocks[2].bn.c0.dxdxi=dxn
blocks[2].bs.c0.dxdxi=dxs
blocks[2].update()
