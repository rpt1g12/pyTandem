#%%
from lib.myGrid import *
from lib.myPlots import *

fig2,ax2=getFig('Fig2')
ax2.set_aspect('equal')

fig,ax=getFig('Fig1')
ax.set_aspect('equal')

mbk=3
nbkx=3;nbky=2;blocks=[]

lxibk=([120,320,360])
letbk=([120,360])

corners=[]
lines=[]
for j in range(nbky):
    for i in range(nbkx):
        corners.append([])
        lines.append([])

def replot(axes=None,doFit=False):
    if (axes==None):
        axes=plt.gca()
    for l in axes.lines:
        l.getPoints()
    if (doFit):
        fit(axes)

#%%        
bkx=0;bky=0;bkn=(bky*bkx)+bkx

corners[bkn].append(corner(-7.0,-7.0))
corners[bkn].append(corner(-3.5,-7.0))
corners[bkn].append(corner(-7.0,0.0))
corners[bkn].append(corner(-0.5,0.0))

lines[bkn].append(line(lxibk[bkx],corners[bkn][0],corners[bkn][1],0))
lines[bkn].append(line(lxibk[bkx],corners[bkn][2],corners[bkn][3],0))
lines[bkn].append(line(letbk[bky],corners[bkn][0],corners[bkn][2],1))
lines[bkn].append(line(letbk[bky],corners[bkn][1],corners[bkn][3],1))

#%%  
bkx=1;bky=0;bkn=(bky*bkx)+bkx

corners[bkn].append(corners[bkn-1][1].clone())
corners[bkn].append(corner(0.5,-7.0))
corners[bkn].append(corners[bkn-1][3].clone())
corners[bkn].append(corner(0.5,0.0))

lines[bkn].append(line(lxibk[bkx],corners[bkn][0],corners[bkn][1],0))
lines[bkn].append(line(lxibk[bkx],corners[bkn][2],corners[bkn][3],0,naca=-21.0))
lines[bkn].append(lines[bkn-1][3])
lines[bkn].append(line(letbk[bky],corners[bkn][1],corners[bkn][3],1))

#%%
bkx=2;bky=0;bkn=(bky*bkx)+bkx

corners[bkn].append(corners[bkn-1][1].clone())
corners[bkn].append(corner(12.0,-7.0))
corners[bkn].append(corners[bkn-1][3].clone())
corners[bkn].append(corner(12.0,0.0))

lines[bkn].append(line(lxibk[bkx],corners[bkn][0],corners[bkn][1],0))
lines[bkn].append(line(lxibk[bkx],corners[bkn][2],corners[bkn][3],0))
lines[bkn].append(lines[bkn-1][3])
lines[bkn].append(line(letbk[bky],corners[bkn][1],corners[bkn][3],1))
#%%
for i in range(mbk):
    blocks.append(block(lines[i]))
    for l in lines[i]:
        l.set_marker('o')
        l.set_markerfacecolor('Red')
        ax.add_line(l)
           
#%%
blocks[1].bn.c0.dxdxi=0.0075*0.05       
blocks[1].bw.c1.dydet=blocks[1].bn.c0.dxdxi*np.sqrt(2)*0.5
blocks[1].be.c1.dydet=blocks[1].bw.c1.dydet
blocks[1].bs.c0.dxdxi=blocks[1].bn.c0.dxdxi
blocks[1].bs.c1.dxdxi=blocks[1].bn.c1.dxdxi

blocks[1].updateLines()

#%%
blocks[0].bw.c1.dydet=blocks[0].be.c1.dydet
blocks[0].bn.c1.dxdxi=blocks[1].bn.c0.dxdxi
blocks[0].bs.c1.dxdxi=blocks[1].bs.c0.dxdxi

blocks[0].be.c0.dxdy=0.0
blocks[0].be.c1.dxdy=1.0

blocks[0].updateLines()

blocks[0].bw.c0.dydet=max(blocks[0].bw.dydet)
blocks[0].be.c0.dydet=max(blocks[0].be.dydet)
blocks[0].bs.c0.dxdxi=max(blocks[0].bs.dxdxi)
blocks[0].bn.c0.dxdxi=max(blocks[0].bn.dxdxi)

i=0;err=1.0
while (err>0.0001):
    blocks[0].bs.c0.dxdxi=max(blocks[0].bs.dxdxi)
    blocks[0].updateLines()
    err=abs(blocks[0].bs.c0.dxdxi/max(blocks[0].bs.dxdxi)-1.0)
    i+=1;print(i)

i=0;err=1.0
while (err>0.0001):
    blocks[0].bw.c0.dydet=max(blocks[0].bw.dydet)
    blocks[0].updateLines()
    err=abs(blocks[0].bw.c0.dydet/max(blocks[0].bw.dydet)-1.0)
    i+=1;print(i)

#%%
blocks[2].bn.c0.dxdxi=blocks[1].bn.c1.dxdxi
blocks[2].bs.c0.dxdxi=blocks[1].bs.c1.dxdxi
blocks[2].be.c1.dydet=blocks[2].bw.c1.dydet

blocks[2].updateLines()

blocks[2].bw.c0.dydet=max(blocks[2].bw.dydet)
blocks[2].be.c0.dydet=max(blocks[2].be.dydet)
blocks[2].bs.c1.dxdxi=max(blocks[2].bs.dxdxi)
blocks[2].bn.c1.dxdxi=max(blocks[2].bn.dxdxi)
#%%
cll(ax2)
for b in blocks:
    b.draw(ax2,update=True)
replot(ax,True)

#%%    
#fh=open('Mblk.xyz','wb')
#hdr=[len(blocks)]
#for b in blocks:
#    hdr.append(b.lxi)
#    hdr.append(b.let)
#    hdr.append(b.lze)
#hdr=np.int32(np.array(hdr))
#fh.write(hdr)
#for b in blocks:
#    fh.write(np.float32(np.transpose(b.x).copy(order='C')))
#    fh.write(np.float32(np.transpose(b.y).copy(order='C')))
#    fh.write(np.float32(np.transpose(b.z).copy(order='C')))
#fh.close()

#%%
#i=-1
#for b in blocks:
#    i+=1    
#    b.update()
#    b.draw(ax2)
#    b.write(i)
#fit(ax2)
#replot(ax,True)

