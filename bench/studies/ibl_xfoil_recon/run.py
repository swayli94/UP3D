"""侦察图：松耦合 IBL vs XFOIL。★ 侦察产物，写 scratchpad（未提交）。"""
import sys,csv,pickle,os; sys.path.insert(0,"/home/lrz/codes/UP3D")
import numpy as np
SC=os.path.dirname(os.path.abspath(__file__))
CACHE=os.path.join(SC,"d13_data.pkl")
M,A,RE=0.5,2.0,3.0e6
ref={"upper":[],"lower":[]}
for r in csv.DictReader(open("cases/reference_data/naca0012_viscous_xfoil/"
                             "delta_star_cf_alpha2_m05_xtr005.csv")):
    ref[r["surface"]].append((float(r["x_c"]),float(r["dstar_over_c"]),float(r["cf"])))
ref={k:np.array(sorted(v)) for k,v in ref.items()}

if os.path.exists(CACHE):
    data=pickle.load(open(CACHE,"rb")); print("  用缓存")
else:
    from pyfp3d.mesh.reader import read_mesh
    from pyfp3d.mesh.wake_cut import cut_wake
    from pyfp3d.post.surface import wall_force_coefficients
    from pyfp3d.viscous import closures as C
    from pyfp3d.viscous.coupling import (CouplingConfig, build_airfoil_case,
                                         make_picard_lifting_driver, run_loose_coupling)
    data={}
    for lvl in ("coarse","medium"):
        mc,wc=cut_wake(read_mesh(f"cases/meshes/naca0012_2.5d/{lvl}.msh"))
        cfg=CouplingConfig(re_chord=RE,m_inf=M,alpha_deg=A,n_outer_max=10)
        case=build_airfoil_case(mc.nodes,mc.elements,mc.boundary_faces["wall"],cfg)
        res=run_loose_coupling(make_picard_lifting_driver(mc,wc,M,A),case,cfg)
        st=case.stations
        xcn=np.asarray(st.xc)[np.asarray(st.station_of)]
        yn=np.asarray(st.xy)[np.asarray(st.station_of)][:,1]
        sn=np.asarray(st.side_node)
        mid=(xcn>0.3)&(xcn<0.7); up=+1 if np.mean(yn[mid&(sn==1)])>0 else -1
        dz=float(np.ptp(mc.nodes[:,2]))
        f=wall_force_coefficients(mc.nodes,mc.elements,mc.boundary_faces["wall"],
            np.asarray(res.phi),alpha_deg=A,u_inf=1.0,s_ref=dz,m_inf=M)
        d={"cl":float(f["cl"]),"n_outer":int(res.n_outer),"converged":bool(res.converged)}
        for k,sv in (("upper",up),("lower",-up)):
            m=sn==sv; x=xcn[m]; o=np.argsort(x)
            d[k]=(x[o],res.outs[m,C.OUT_DS1][o],res.outs[m,C.OUT_CF1][o])
        data[lvl]=d
        print(f"  {lvl}: cl={d['cl']:.4f} 外迭代 {d['n_outer']}")
    pickle.dump(data,open(CACHE,"wb"))

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig,ax=plt.subplots(2,2,figsize=(12.5,8.4))
CL={"coarse":"C0","medium":"C1"}
# (1) δ*
for k,ls in (("upper","-"),("lower","--")):
    ax[0,0].plot(ref[k][:,0],ref[k][:,1],"k"+ls,lw=1.8,
                 label=f"XFOIL {k}" if True else None)
    for lvl in data:
        x,ds,_=data[lvl][k]
        ax[0,0].plot(x,ds,ls,color=CL[lvl],lw=1.2,
                     label=f"IBL {lvl} {k}")
ax[0,0].set_xlabel("x/c"); ax[0,0].set_ylabel(r"$\delta^*/c$")
ax[0,0].set_title(r"(1) displacement thickness $\delta^*/c$")
ax[0,0].legend(fontsize=7,ncol=2); ax[0,0].grid(alpha=.3)

# (2) 比值 —— 本轮的发现
g=np.linspace(0.05,1.0,80)
for k,ls in (("upper","-"),("lower","--")):
    for lvl in data:
        x,ds,_=data[lvl][k]
        r_=np.interp(g,x,ds)/np.interp(g,ref[k][:,0],ref[k][:,1])
        ax[0,1].plot(g,r_,ls,color=CL[lvl],lw=1.4,label=f"{lvl} {k}")
ax[0,1].axhline(1.0,color="k",lw=1.2); ax[0,1].axhline(0.5,color="r",ls=":",lw=1)
ax[0,1].set_xlabel("x/c"); ax[0,1].set_ylabel(r"$\delta^*_{\rm IBL}/\delta^*_{\rm XFOIL}$")
ax[0,1].set_title("(2) THE FINDING: monotone downstream drift, not scatter")
ax[0,1].set_ylim(0.3,1.3); ax[0,1].legend(fontsize=7,ncol=2); ax[0,1].grid(alpha=.3)

# (3) c_f
for k,ls in (("upper","-"),("lower","--")):
    ax[1,0].plot(ref[k][:,0],ref[k][:,2],"k"+ls,lw=1.8,label=f"XFOIL {k}")
    for lvl in data:
        x,_,cf=data[lvl][k]
        ax[1,0].plot(x,cf,ls,color=CL[lvl],lw=1.2,label=f"IBL {lvl} {k}")
ax[1,0].set_xlabel("x/c"); ax[1,0].set_ylabel(r"$c_f$"); ax[1,0].set_ylim(0,0.012)
ax[1,0].set_title(r"(3) skin friction $c_f$   (trip at x/c = 0.05, both surfaces)")
ax[1,0].legend(fontsize=7,ncol=2); ax[1,0].grid(alpha=.3)

# (4) cl
vals=[("XFOIL inviscid",0.2921,"0.5"),("IBL coarse",data["coarse"]["cl"],"C0"),
      ("IBL medium",data["medium"]["cl"],"C1"),("XFOIL viscous",0.2691,"k")]
ax[1,1].barh([v[0] for v in vals],[v[1] for v in vals],
             color=[v[2] for v in vals],height=.55)
for i,(n,v,_) in enumerate(vals):
    ax[1,1].text(v+0.001,i,f"{v:.4f}",va="center",fontsize=9)
ax[1,1].axvline(0.2691,color="k",ls=":",lw=1); ax[1,1].axvline(0.2921,color="0.5",ls=":",lw=1)
ax[1,1].set_xlim(0.25,0.305); ax[1,1].set_xlabel("cl")
ax[1,1].set_title("(4) lift: refinement moves AWAY from XFOIL, toward inviscid")
ax[1,1].grid(alpha=.3,axis="x")
fig.suptitle("Reconnaissance: loose-coupled IBL vs XFOIL 6.99  "
             "(NACA0012, M 0.5, Re 3e6, alpha 2, x_tr 0.05)",fontsize=11)
fig.tight_layout()
out=os.path.join(SC,"d13_recon.png"); fig.savefig(out,dpi=125)
print(f"  ✓ {out}")
