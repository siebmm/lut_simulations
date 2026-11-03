"""Apply the explicit COMSOL-to-SAFT explanation consistently to active papers."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXT = HERE.parent/'comsol_chapter_20261101'


def edit(path, old, new):
    text=path.read_text(encoding='utf-8')
    if new in text:return
    assert old in text, (path.name,old[:100])
    path.write_text(text.replace(old,new,1),encoding='utf-8',newline='\n')


def workflow(scope):
    return r'''
\begin{table}[htbp]\centering\small
\caption{Computational workflow. Heat transfer is solved locally for source calibration; the production weld simulations solve elastic propagation, and SAFT is applied to exported recordings.}
\label{SCOPE:workflow}
\begin{tabularx}{\linewidth}{@{}p{.20\linewidth}XX@{}}\toprule
Stage & Implementation & Output used next\\\midrule
Source calibration & COMSOL Heat Transfer in Solids and Solid Mechanics; heating drives thermal expansion & A mechanical surface traction calibrated against outgoing velocity waveforms\\
Weld propagation & COMSOL Elastic Waves, Time Explicit; full/partial geometry and independent source shots & Surface-velocity histories containing direct, reflected, converted and scattered waves\\
Image reconstruction & Custom C++ SAFT; Python study and figure scripts & Path-specific amplitude images and separate measurement-support maps\\\bottomrule
\end{tabularx}
\end{table}
'''.replace('SCOPE',scope)


def main():
    c=HERE/'chapter_comsol_core.tex'
    edit(c,r'\chapter{Laser-ultrasonic imaging of full and partial welds with a common SAFT workflow}',r'\chapter{COMSOL modelling and SAFT imaging of full and partial welds}')
    edit(c,'A calibrated laser-equivalent source drives a two-dimensional elastic steel model; five independent source positions and 337 or 331 ideal surface-velocity receivers provide the acquisition.',
         'COMSOL Multiphysics 6.4 simulates a calibrated laser-equivalent source in a two-dimensional elastic steel model. Five source positions and 337 or 331 ideal surface-velocity receivers provide time histories for a separate C++ SAFT reconstruction.')
    edit(c,'The model, calibration and reconstruction are introduced before the results. The detailed simulation records and workflow investigations remain in \\corefullappendix{} and \\corepartialappendix{}. The present chapter concentrates on useful images, the physical reasons for choosing their paths and the acquisition arrangements required to sustain them.',
         'The numerical experiment has three stages: calibrate the source in COMSOL, simulate its propagation through each weld in COMSOL, and focus the exported recordings in a separate SAFT implementation. COMSOL predicts the interacting wave field; SAFT tests whether selected arrivals are consistent with scattering at an image pixel. The detailed simulation records and workflow investigations remain in \\corefullappendix{} and \\corepartialappendix{}.'+workflow('core'))
    edit(c,r'\subsection{Elastic calculation and recorded quantities}'+'\nThe model solves small-strain linear elastodynamics.',
         r'\subsection{COMSOL wave-propagation model and signal export}'+'\nThe production calculations use the \\emph{Elastic Waves, Time Explicit} interface in COMSOL Multiphysics 6.4 (build 429). Each weld is a fixed elastic geometry during the 50 $\\mu$s acoustic calculation; welding heat flow, deposition and melt evolution are outside this model. The underlying physics is small-strain linear elastodynamics.')
    edit(c,'The COMSOL elastic-wave formulation uses fourth-order velocity and strain fields and explicit fourth-order Runge--Kutta integration.',
         'COMSOL advances fourth-order velocity and strain fields using explicit fourth-order Runge--Kutta integration; the displacement equations above state the equivalent continuum physics.')
    edit(c,'A mechanical surface load is used to represent the elastic consequence of laser heating. Its amplitude was calibrated independently of the weld images against a coupled thermoelastic calculation on an 8 mm by 4 mm flat patch, reduced by symmetry.',
         'A mechanical surface load represents the elastic consequence of laser heating. Its amplitude was calibrated in COMSOL independently of the weld images, using \\emph{Heat Transfer in Solids} coupled one way through thermal expansion to \\emph{Solid Mechanics} on an 8 mm by 4 mm flat patch, reduced by symmetry. A second mechanical field applies the candidate traction on the same mesh, allowing the outgoing waveforms to be compared.')
    edit(c,r'\section{A common basic SAFT workflow}',r'\section{SAFT reconstruction outside COMSOL}')
    edit(c,'A receiver provides a time trace; SAFT assigns a score to each trial image pixel.',
         'COMSOL exports the normal velocity at every receiver after each independent source shot. The C++ reconstruction reads these traces and assigns a score to every trial image pixel. Choosing an imaging path changes its delay law and accepted channels; it does not rerun COMSOL or remove other wave modes from the simulated recording.')
    edit(c,'Pulse screening addresses a different issue. An otherwise valid focusing time can coincide with a strong direct L, S or Rayleigh arrival.',
         'The simulated recording contains strong motion even without a flaw: direct bulk waves, surface waves and boundary returns are part of the elastic solution. SAFT can map such motion to a false scattering position when a predicted delay samples the same packet across receivers. This explains source-related image background without attributing every feature to a solver error. Pulse screening addresses the direct part of this ambiguity: an otherwise valid focusing time can coincide with a strong L, S or Rayleigh arrival.')
    edit(c,'The 0.1 mm image step samples a wavefield', 'The finite-element mesh controls the accuracy of COMSOL wave propagation; the SAFT grid only samples the reconstructed image. The 0.1 mm image step samples a wavefield')

    m=EXT/'chapter_comsol.tex'
    edit(m,r'\chapter{Numerical feasibility of laser-ultrasonic weld imaging}',r'\chapter{COMSOL modelling of laser-ultrasonic weld inspection}')
    edit(m,'This study addresses the acoustic part of that question. It simulates laser-equivalent generation and ideal normal-surface-velocity reception, then reconstructs images from the recordings.',
         'COMSOL Multiphysics 6.4 is used to model the acoustic part of this inspection: laser-equivalent excitation, elastic propagation, interaction with the weld and LOF, and ideal normal-surface-velocity reception. A separate C++ SAFT implementation reconstructs images from the exported recordings.')
    edit(m,'Appendices~\\ref{comsol:app_full} and~\\ref{comsol:app_partial} contain the complete standalone simulation papers.',
         workflow(r'\comsolScope')+'\nAppendices~\\ref{comsol:app_full} and~\\ref{comsol:app_partial} contain the complete standalone simulation papers.')
    edit(m,'A flat 8 mm by 4 mm patch, reduced by symmetry, received Gaussian absorbed heating',
         'In COMSOL, Heat Transfer in Solids supplies the temperature field to Solid Mechanics through thermal expansion. A flat 8 mm by 4 mm patch, reduced by symmetry, receives Gaussian absorbed heating')

    detail=EXT/'comsol_model_detail.tex'
    edit(detail,r'\subsection{Governing equations, boundary conditions and numerical solution}',r'\subsection{COMSOL physics, boundary conditions and numerical solution}')
    edit(detail,'The mechanical calculation solves small-strain, linear elastodynamics in the cross-section.',
         'The production models use COMSOL Multiphysics 6.4 (build 429), with the Elastic Waves, Time Explicit interface. Each full or partial weld is a fixed geometry during the acoustic calculation; the simulation does not model deposition, melting or a transient welding temperature field. The mechanical calculation solves small-strain, linear elastodynamics in the cross-section.')
    edit(detail,'Exposed plate, backwall and groove surfaces','Exposed weld, side, backwall and groove surfaces')
    edit(detail,'not the elastic properties of weld metal relative to plate','not the elastic properties of weld metal relative to surrounding steel')
    edit(detail,'Quantity & Specification\\\\\\midrule','Quantity & Specification\\\\\\midrule\nSoftware / production interface & COMSOL Multiphysics 6.4 (build 429); Elastic Waves, Time Explicit\\\\')
    edit(detail,'Reconstruction studies reuse those records;','The exported receiver histories are the input to the separate SAFT calculation. Reconstruction studies reuse those records;')

    model=EXT/'paper_model.tex'
    intro=r'''\section{Computational approach: COMSOL simulation and SAFT analysis}
The study separates source generation, wave propagation and imaging. COMSOL Multiphysics 6.4 (build 429) solves the finite-element models. A local Heat Transfer in Solids calculation drives thermal expansion in Solid Mechanics and calibrates an equivalent mechanical traction. The production weld models then apply that traction in Elastic Waves, Time Explicit. Their surface-velocity exports are processed by a custom C++ SAFT implementation, with Python scripts organizing comparisons and publication figures.

This sequence separates what is simulated from what is inferred. COMSOL predicts the interacting wave field, including structural waves and flaw scattering. SAFT does not solve the wave equation again: it aligns recorded samples under a chosen source-to-pixel-to-receiver hypothesis. A bright image pixel therefore requires checks on competing arrivals and measurement support. The full and partial welds are fixed elastic geometries, not time-evolving simulations of deposition or melting.
'''+workflow(r'\comsolScope')+'\n'
    edit(model,r'\section{Specimen and simulated measurement}',intro+r'\section{Specimen and simulated measurement}')
    edit(model,r'\section{How the wave field is simulated}',r'\section{COMSOL model of elastic-wave propagation}')
    edit(model,'The saved production models use fourth-order velocity and strain fields and explicit fourth-order Runge--Kutta integration.',
         'The Elastic Waves, Time Explicit interface advances fourth-order velocity and strain fields using explicit fourth-order Runge--Kutta integration. The displacement equations above describe the equivalent elastic physics; the production solution directly supplies velocity.')
    edit(model,'Heat conduction is coupled one way to linear thermal expansion.',
         'The COMSOL Heat Transfer in Solids interface is coupled one way through thermal expansion to Solid Mechanics. A separate Solid Mechanics field applies the candidate dipole on the same mesh.')
    edit(model,'Weld metal and plate share','Weld metal and surrounding steel share')
    edit(model,'material contrast between weld and plate','material contrast between weld and surrounding steel')
    method=EXT/'paper_method.tex'
    edit(method,r'\section{Reconstructing an image from the recordings}',r'\section{SAFT reconstruction from the COMSOL recordings}')
    edit(method,'The calibrated elastic simulation supplies the complete wave field; these imaging hypotheses select parts of that field for focusing.',
         'COMSOL supplies complete receiver waveforms containing all simulated arrivals. The separate C++ reconstruction tests different timing hypotheses on those same recordings; choosing a pure-mode path does not suppress other modes in COMSOL. The image grid is a postprocessing grid, distinct from the finite-element mesh.')
    for kind,name in [('full','a full weld'),('partial','a partly filled weld')]:
        p=EXT/(kind+'_simulation_paper.tex')
        old=r'\chapter{'+('Full-weld' if kind=='full' else 'Partial-weld')+' laser-ultrasonic imaging}'
        edit(p,old,r'\chapter{COMSOL simulation and SAFT imaging of '+name+'}')
    edit(EXT/'full_simulation_paper.tex','The calibrated two-dimensional elastic model supplies complete waveforms; synthetic aperture focusing tests direct and reflected propagation paths.',
         'COMSOL Multiphysics 6.4 supplies the waveforms from a calibrated two-dimensional elastic model; a separate C++ synthetic aperture focusing calculation tests direct and reflected propagation paths.')
    edit(EXT/'partial_simulation_paper.tex','Five source positions and 331 receivers are simulated in paired flawed and pristine specimens; three sources lie on the exposed floor of a 9 mm fill.',
         'COMSOL Multiphysics 6.4 simulates five source positions and 331 receivers in paired flawed and pristine specimens; three sources lie on the exposed floor of a 9 mm fill. A separate C++ SAFT implementation reconstructs images from the exported velocity histories.')
    print('Updated COMSOL roles, named interfaces, workflow and interpretation in four active papers.')


if __name__=='__main__':main()
