"""Compare representations using identical events, not different cosmologies."""

import csv
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from spacetime.diagrams import event_coordinates, visible_ring
from spacetime.model import ClosedDust

FIGURES = ROOT / "paper/figures"
RESULTS = ROOT / "results"
BLUE, ORANGE, PURPLE, GRAY = "#2166ac", "#b35806", "#8e3b7a", "#9ca8b2"
EO = 1.8
B = ClosedDust()
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "svg.hashsalt": "geometry-diagram-comparison", "savefig.dpi": 180})


def save(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES/f"{name}.pdf", bbox_inches="tight",
                metadata={"CreationDate": None, "ModDate": None})
    svg = FIGURES/f"{name}.svg"
    fig.savefig(svg, bbox_inches="tight", metadata={"Date": None})
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines())+"\n")
    fig.savefig(FIGURES/f"{name}.png", bbox_inches="tight")
    plt.close(fig)


def comparison():
    eta = np.linspace(.035, EO, 700)
    chi = np.linspace(-1.85, 1.85, 400)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9.6))
    fig.subplots_adjust(hspace=.36, wspace=.28, bottom=.14, top=.90)
    specs = [
        ("comoving", "A  Cosmic time and comoving position", r"$\chi$ (rad)", r"$ct/A_{\max}$"),
        ("proper", "B  Cosmic time and proper separation", r"$D/A_{\max}= (A/A_{\max})\chi$", r"$ct/A_{\max}$"),
        ("conformal", "C  Conformal spacetime", r"$\chi$ (rad)", r"$\eta$"),
        ("nested", "D  Nested layers: the same spatial patch", r"$X/A_o$", r"$Y/A_o$"),
    ]
    for ax, (key, title, xlabel, ylabel) in zip(axes.flat, specs):
        for epoch in [.4, .9, 1.4, EO]:
            xy = event_coordinates(epoch, chi)[key]
            ax.plot(*xy.T, color=GRAY, lw=.9, alpha=.65)
        for location in [-1.5, -.9, 0, .9, 1.5]:
            xy = event_coordinates(eta, location)[key]
            ax.plot(*xy.T, color=PURPLE if location==.9 else GRAY,
                    lw=1.2 if location==.9 else .8, ls="--")
        for sign in [-1, 1]:
            xy = event_coordinates(eta, sign*(EO-eta))[key]
            ax.plot(*xy.T, color=BLUE, lw=2.2)
        for epoch, location, label, filled in [(.9,.9,"E",True),(EO,.9,"N",False)]:
            xy = event_coordinates(epoch, location)[key]
            ax.plot(*xy, marker="o", mec=PURPLE, mfc=PURPLE if filled else "white", ms=7, zorder=5)
            ax.annotate(label, xy, xytext=(6,8), textcoords="offset points", color=PURPLE, weight="bold")
        xy = event_coordinates(EO, 0)[key]
        ax.plot(*xy, "*", color="black", ms=11, zorder=5)
        ax.annotate("O", xy, xytext=(-15,8), textcoords="offset points", weight="bold")
        ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
        if key in ["comoving", "conformal"]:
            ax.set_xlim(-1.9,1.9)
        if key=="proper": ax.set_xlim(-1.2,1.2)
        if key in ["proper", "comoving"]: ax.set_ylim(-.015,.46)
        if key=="conformal":
            ax.set_ylim(-.05,2.05)
            ax.set_aspect("equal", adjustable="box")
        if key=="nested":
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(-.40,1.23); ax.set_ylim(-1.1,1.1)
    handles=[Line2D([0],[0],color=BLUE,lw=2,label="Photons arriving at O"),
             Line2D([0],[0],color=PURPLE,ls="--",label="One galaxy's history"),
             Line2D([0],[0],color=GRAY,lw=1,label="Equal cosmic-time layers")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5,.055), ncol=3, frameon=False)
    fig.suptitle("One universe, one reception event, four drawings", fontsize=17)
    fig.text(.5,.018,"E = galaxy when it emitted the light received at O.   N = the same galaxy now; N is not seen at O.",
             ha="center", fontsize=11)
    save(fig,"07-diagram-comparison")


def slices_and_rings():
    fig=plt.figure(figsize=(11.5,5.7))
    ax=fig.add_subplot(121,projection="3d")
    polar=fig.add_subplot(122)
    longitude=np.linspace(0,2*np.pi,70)
    colatitude=np.linspace(0,np.pi,25)
    uu,vv=np.meshgrid(colatitude,longitude)
    eta=np.linspace(.08,EO,300)
    radius=B.A(eta)/B.A(EO)
    for epoch in [.6,1.2,EO]:
        a=B.A(epoch)/B.A(EO)
        ax.plot_wireframe(a*np.cos(uu),a*np.sin(uu)*np.cos(vv),a*np.sin(uu)*np.sin(vv),
                          color=GRAY,alpha=.17,rstride=6,cstride=4,linewidth=.55)
        theta=np.linspace(-np.pi,np.pi,500)
        polar.plot(a*np.cos(theta),a*np.sin(theta),color=GRAY,lw=.9)
        if epoch<EO:
            ring=visible_ring(epoch,longitude)
            ax.plot(*ring.T,color=ORANGE,lw=2.5)
            for sign in [-1,1]:
                xy=event_coordinates(epoch,sign*(EO-epoch))["nested"]
                polar.plot(*xy,"o",color=ORANGE,ms=6)
    for bearing in np.linspace(0,2*np.pi,5)[:-1]:
        ax.plot(radius*np.cos(EO-eta),radius*np.sin(EO-eta)*np.cos(bearing),
                radius*np.sin(EO-eta)*np.sin(bearing),color=BLUE,lw=1.6)
    for sign in [-1,1]:
        xy=event_coordinates(eta,sign*(EO-eta))["nested"]
        polar.plot(*xy.T,color=BLUE,lw=2)
    ax.scatter([1],[0],[0],color="black",marker="*",s=75)
    ax.text(1.12,0,.12,"O",weight="bold")
    ax.set(xlabel=r"$X_0/A_o$",ylabel=r"$X_1/A_o$",zlabel=r"$X_2/A_o$")
    ax.set(xticks=[-1,0,1],yticks=[-1,0,1],zticks=[-1,0,1])
    ax.set_box_aspect((1,1,1)); ax.view_init(elev=22,azim=35)
    ax.set_title("A section of the nested spatial slices",pad=15)
    polar.plot(1,0,"*",color="black",ms=11)
    polar.annotate("O",(1,0),xytext=(8,8),textcoords="offset points")
    polar.annotate("Only these intersections\nare visible from O",(.36,.25),xytext=(-.98,.72),
                   arrowprops={"arrowstyle":"->","color":ORANGE},color=ORANGE,fontsize=10)
    polar.set(xlabel=r"$X_0/A_o$",ylabel=r"$X_1/A_o$",aspect="equal",xlim=(-1.12,1.18),ylim=(-1.05,1.06),
              title="The corresponding meridian section")
    fig.subplots_adjust(bottom=.23,top=.83,wspace=.28)
    fig.suptitle("A time layer is larger than the set of events we see on it",fontsize=16)
    fig.text(.75,.105,"Gray circle: section of a whole spatial slice.\nOrange dots: section of its visible subset.",ha="center",fontsize=10)
    fig.text(.5,.025,"The full universe is S³. Meshes show S² sections; orange rings are sections of the visible S² subsets.",
             ha="center",fontsize=10)
    save(fig,"08-slices-and-visible-subsets")


def stress_tests():
    fig,axes=plt.subplots(1,3,figsize=(14,4.8))
    eta=np.geomspace(.004,EO,500)
    axes[0].loglog(eta/EO,eta/EO,color=BLUE,label="Conformal time")
    axes[0].loglog(eta/EO,B.A(eta)/B.A(EO),color=ORANGE,label="Nested radius")
    axes[0].loglog(eta/EO,B.time(eta)/B.time(EO),color=PURPLE,label="Cosmic time")
    axes[0].set(xlabel=r"Conformal epoch $\eta/\eta_o$",ylabel="Fraction of full linear axis span",
                title="Early epochs occupy very little radius",xlim=(.002,1),ylim=(1e-8,1))
    axes[0].legend(frameon=False,fontsize=9,loc="lower right")
    eta=np.linspace(.001,2*np.pi-.001,500)
    axes[1].plot(eta[:250],B.A(eta[:250]),color=BLUE,label="Expansion")
    axes[1].plot(eta[249:],B.A(eta[249:]),color=ORANGE,label="Contraction")
    pair=np.array([np.pi/2,3*np.pi/2])
    axes[1].plot(pair,[.5,.5],"o",color=PURPLE)
    axes[1].axhline(.5,color=GRAY,ls=":",lw=1)
    for point,label in zip(pair,["P","Q"]): axes[1].annotate(label,(point,.5),xytext=(6,8),textcoords="offset points",color=PURPLE)
    axes[1].set(xlabel=r"Conformal time $\eta$",ylabel=r"$A/A_{\max}$",title="Two epochs can have the same radius",
                xticks=[0,np.pi/2,np.pi,3*np.pi/2,2*np.pi],xticklabels=["0",r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"])
    axes[1].legend(frameon=False,fontsize=9)
    theta=np.linspace(0,2*np.pi,600)
    axes[2].plot(np.cos(theta),np.sin(theta),color=GRAY,lw=.8)
    axes[2].plot(.5*np.cos(theta),.5*np.sin(theta),color=BLUE,lw=3)
    axes[2].plot(.5*np.cos(theta),.5*np.sin(theta),color=ORANGE,lw=2,ls="--")
    xy=.5*np.array([np.cos(.8),np.sin(.8)])
    axes[2].plot(*xy,"o",color=PURPLE,ms=8)
    axes[2].annotate("P and Q coincide",xy,xytext=(-.95,.85),
                     arrowprops={"arrowstyle":"->","color":PURPLE},color=PURPLE)
    axes[2].set(aspect="equal",xlabel=r"$X/A_{\max}$",ylabel=r"$Y/A_{\max}$",
                title="The nested drawing loses time ordering",xlim=(-1.1,1.1),ylim=(-1.1,1.1))
    fig.subplots_adjust(top=.80,bottom=.18,wspace=.40)
    fig.suptitle("Two limitations of using expansion as the drawing's clock",fontsize=16)
    fig.text(.5,.025,"Left: coordinate compression, not measured reader performance.   Right: the full recollapsing dust history, beyond the draft's expanding interval.",ha="center",fontsize=10)
    save(fig,"09-diagram-stress-tests")


def diagnostics():
    RESULTS.mkdir(exist_ok=True)
    eta_early=2*np.arcsin(np.sqrt(B.A(EO)/1100))
    ratios={"nested_radius":1/1100,"conformal_time":eta_early/EO,
            "cosmic_time":float(B.time(eta_early)/B.time(EO))}
    pair=event_coordinates([np.pi/2,3*np.pi/2],[.8,.8])
    report={"scope":"Analytic/graphical comparison; no human participants or learning-outcome data.",
            "eta_observer":EO,"illustrative_1_plus_z":1100,
            "early_eta":float(eta_early),"early_axis_span_fractions":ratios,
            "early_pixels_for_500_pixel_linear_span":{k:500*v for k,v in ratios.items()},
            "null_tangent_angle_to_radial_degrees_at_reception":float(np.degrees(EO/2)),
            "recollapse_nested_event_separation":float(np.linalg.norm(pair["nested"][1]-pair["nested"][0])),
            "recollapse_cosmic_time_separation_Amax_over_c":float(pair["comoving"][1,1]-pair["comoving"][0,1])}
    (RESULTS/"diagram-comparison.json").write_text(json.dumps(report,indent=2)+"\n")
    with (RESULTS/"diagram-events.csv").open("w",newline="") as stream:
        writer=csv.writer(stream,lineterminator="\n")
        writer.writerow(["event","eta","chi_rad","ct_over_Amax","D_over_Amax","X_over_Ao","Y_over_Ao"])
        for name,eta,chi in [("E",.9,.9),("N",EO,.9),("O",EO,0)]:
            e=event_coordinates(eta,chi)
            writer.writerow([name,eta,chi,e["comoving"][1],e["proper"][0],*e["nested"]])
    generated=ROOT/"paper/generated/diagram-numbers.tex"
    generated.write_text("% Generated by scripts/compare_diagrams.py\n"
        +f"\\newcommand{{\\NestedPixels}}{{{500*ratios['nested_radius']:.2f}}}\n"
        +f"\\newcommand{{\\ConformalPixels}}{{{500*ratios['conformal_time']:.2f}}}\n"
        +f"\\newcommand{{\\CosmicPixels}}{{{500*ratios['cosmic_time']:.3f}}}\n"
        +f"\\newcommand{{\\NullAngle}}{{{np.degrees(EO/2):.2f}}}\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    comparison()
    slices_and_rings()
    stress_tests()
    diagnostics()
