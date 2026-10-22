# Figure 1: model structure schematic (Python 3, matplotlib). Run: python3 figure1_schematic.py
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10})
INK="#1f2328"
C={"animal":("#f6dccb","#a4521f"),"env":("#dcebd3","#3c7a2e"),"human":("#d3e5f5","#2b6cb0"),"ext":("#e9e0f5","#6b46a6")}
TR="#3d4247"; SH="#3c7a2e"; EXT="#6b46a6"
fig,ax=plt.subplots(figsize=(10,7.4)); ax.set_xlim(0,10); ax.set_ylim(-1.9,8.2); ax.axis("off")
N={}
def box(k,x,y,t,kind,w=1.9,h=0.72):
    fc,ec=C[kind]
    ax.add_patch(FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.02,rounding_size=0.08",fc=fc,ec=ec,lw=1.4,zorder=3))
    ax.text(x,y,t,ha="center",va="center",fontsize=9.5,color=INK,zorder=4,linespacing=1.15)
    N[k]=(x,y,w,h)
def edge(k,side):
    x,y,w,h=N[k]
    return {"r":(x+w/2,y),"l":(x-w/2,y),"t":(x,y+h/2),"b":(x,y-h/2)}[side] if isinstance(side,str) else side
def arr(a,sa,b,sb,col=TR,ls="-",lw=1.3,rad=0.0,both=False,ms=12):
    p=FancyArrowPatch(edge(a,sa),edge(b,sb),arrowstyle="<|-|>" if both else "-|>",mutation_scale=ms,color=col,lw=lw,
                      linestyle=ls,connectionstyle=f"arc3,rad={rad}",shrinkA=2,shrinkB=2,zorder=2)
    ax.add_patch(p)
AX=3.7
box("BP",AX,7.6,"Backyard poultry\n(SEIR)","animal")
box("WB",AX,5.95,"Wild birds\n(SEIR)","animal")
box("CPL",AX,4.3,"Commercial poultry,\nlarge-scale (SEIR)","animal")
box("CPS",AX,2.65,"Commercial poultry,\nsmall-scale (SEIR)","animal")
box("X",AX,0.75,"Reassortment hosts\n(SEIR)","animal")
box("E",6.6,2.95,"Environmental\nreservoir","env",w=1.7)
box("W",9.0,5.95,"Poultry workers\n(SEIRD)","human")
box("H",9.0,1.55,"General population\n(SEIRD)","human")
box("IMP",1.0,5.95,"Seasonal wild-bird\nimportation","ext",w=1.75)
box("RS",1.0,1.6,"Restocking of\ndepopulated\npremises","ext",w=1.6,h=0.95)
# animal-animal (bidirectional, shared coefficient)
arr("WB","t","BP","b",both=True)
arr("WB","b","CPL","t",both=True)
arr("WB",(AX-0.95,5.75),"CPS",(AX-0.95,2.85),both=True,rad=0.45)
# to workers
arr("BP",(AX+0.95,7.5),"W",(8.6,6.31))
arr("WB","r","W","l")
arr("CPL",(AX+0.95,4.45),"W",(8.3,5.59),lw=2.6)
arr("CPS",(AX+0.95,2.85),"W",(8.55,5.59))
arr("X",(AX+0.95,0.9),"W",(8.85,5.59),rad=0.38)
# to general population
arr("BP",(AX+0.95,7.75),"H",(9.85,1.91),rad=-0.38)
arr("WB",(AX+0.95,5.8),"H",(8.1,1.85),rad=-0.32)
arr("X",(AX+0.95,0.65),"H","l")
arr("W","b","H","t")
# environment -> humans
arr("E",(7.45,3.15),"W",(8.75,5.59))
arr("E",(7.45,2.75),"H",(8.4,1.91))
# shedding: bracket collecting every infectious host
bx=5.05
ax.plot([bx,bx],[0.75,7.6],color=SH,lw=1.2,ls=(0,(2,2)),zorder=1)
for y in [7.6,5.95,4.3,2.65,0.75]:
    ax.plot([AX+0.95,bx],[y-0.22,y-0.22] if False else [y,y],color="none")
for y in [7.6,5.95,4.3,2.65,0.75]:
    ax.add_patch(FancyArrowPatch((AX+0.95,y-0.2),(bx,y-0.2),arrowstyle="-",color=SH,lw=1.2,linestyle=(0,(2,2)),zorder=1))
ax.add_patch(FancyArrowPatch((bx,2.95),(5.75,2.95),arrowstyle="-|>",mutation_scale=12,color=SH,lw=1.6,zorder=2))
ax.add_patch(FancyArrowPatch((9.6,5.59),(7.25,3.31),arrowstyle="-|>",mutation_scale=12,color=SH,lw=1.2,linestyle=(0,(2,2)),connectionstyle="arc3,rad=-0.3",zorder=2))
ax.add_patch(FancyArrowPatch((9.6,1.91),(7.25,2.59),arrowstyle="-|>",mutation_scale=12,color=SH,lw=1.2,linestyle=(0,(2,2)),connectionstyle="arc3,rad=0.3",zorder=2))
# one-year extension
arr("IMP","r","WB","l",col=EXT,ls=(0,(5,3)),lw=1.5)
arr("RS",(1.8,1.75),"CPS",(AX-0.95,2.45),col=EXT,ls=(0,(5,3)),lw=1.5)
arr("RS",(1.0,2.08),"CPL",(AX-0.95,4.1),col=EXT,ls=(0,(5,3)),lw=1.5,rad=-0.25)
# legend
h=[Line2D([0],[0],color=TR,lw=1.3,marker=">",markersize=6,label="Transmission / exposure (non-zero β)"),
   Line2D([0],[0],color=TR,lw=2.4,label="Dominant occupational pathway ($\\beta_{CPLW}$)"),
   Line2D([0],[0],color=SH,lw=1.2,ls=(0,(2,2)),label="Viral shedding into the environment"),
   Line2D([0],[0],color=EXT,lw=1.5,ls=(0,(5,3)),label="One-year extension (Eqs. 16–17)")]
ax.legend(handles=h,loc="upper center",bbox_to_anchor=(0.5,0.12),ncol=2,frameon=False,fontsize=9)
ax.text(5,-1.75,"Double-headed arrows: transmission in both directions with a shared coefficient. Each population also transmits within itself (not drawn).\nCommercial poultry and reassortment hosts have no direct pathway to the general population; no pathway links wild birds to reassortment hosts.",
        ha="center",va="center",fontsize=8,color="#57606a")
fig.savefig("Figure1_model_structure.png",dpi=300,bbox_inches="tight",facecolor="white")
fig.savefig("Figure1_model_structure.pdf",bbox_inches="tight",facecolor="white")
