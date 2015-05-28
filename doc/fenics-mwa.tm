<TeXmacs|1.0.7.14>

<style|article>

<\body>
  <doc-data|<doc-title|Minimally Invasive Tumour Ablation
  Model>|<\doc-author-data|<author-name|Sheldon Hall>>
    \;
  </doc-author-data>>

  <no-indent>MITA models user guide.

  \;

  <no-indent>Copyright (C) 2014 Sheldon Hall (sheldon.hall@eng.ox.ac.uk)

  \;

  <no-indent>Permission is granted to copy, distribute and/or modify this
  document under the terms of the GNU Free Documentation License, Version 1.3
  or any later version published by the Free Software Foundation; with no
  Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts. A copy
  of the license is included with this manual. If not, see
  (http://www.gnu.org/licenses/).

  <new-page>

  <\table-of-contents|toc>
    <vspace*|1fn><with|font-series|bold|math-font-series|bold|1<space|2spc>Getting
    Started> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-1><vspace|0.5fn>

    <with|par-left|1.5fn|1.1<space|2spc>Mesh Generation
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-2>>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|2<space|2spc>Code
    Structure> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-3><vspace|0.5fn>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|3<space|2spc>Code
    Components> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-5><vspace|0.5fn>

    <with|par-left|1.5fn|3.1<space|2spc>Apparent heat capacity form of the
    Pennes equation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-6>>

    <with|par-left|3fn|3.1.1<space|2spc>Thermal conductivity
    <with|mode|math|k> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-8>>

    <with|par-left|3fn|3.1.2<space|2spc>Steady-state temperature
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-9>>

    <with|par-left|1.5fn|3.2<space|2spc>Electric Potential Solver (RFA)
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-11>>

    <with|par-left|3fn|3.2.1<space|2spc>Electrical conductivity
    <with|mode|math|\<sigma\>> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-12>>

    <with|par-left|3fn|3.2.2<space|2spc>Control System
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-13>>

    <with|par-left|6fn|Impedance <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-14><vspace|0.15fn>>

    <with|par-left|3fn|3.2.3<space|2spc>Computation of useful quantities.
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-15>>

    <with|par-left|1.5fn|3.3<space|2spc>Vector Helmholtz Solver (MWA)
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-16>>

    <with|par-left|3fn|3.3.1<space|2spc>Scalar Form
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-17>>

    <with|par-left|3fn|3.3.2<space|2spc>Scalar Variational Form
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-18>>

    <with|par-left|3fn|3.3.3<space|2spc>Expanded integrals
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-19>>

    <with|par-left|3fn|3.3.4<space|2spc>Expanding complex-valued functions
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-20>>

    <with|par-left|3fn|3.3.5<space|2spc>Vector Variational Form With
    Inhomogeneous Material Properties <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-21>>

    <with|par-left|3fn|3.3.6<space|2spc>First-Order Absorbing Boundary
    Condition (Scattering) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-22>>

    <with|par-left|3fn|3.3.7<space|2spc>Waveguide Port Boundary Condition
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-23>>

    <with|par-left|3fn|3.3.8<space|2spc>Verification Cases (Axisymmetric)
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-24>>

    <with|par-left|6fn|Coaxial cable analytic solution
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-25><vspace|0.15fn>>

    <with|par-left|6fn|Concentric Cylinders
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-26><vspace|0.15fn>>

    <with|par-left|6fn|Annular Slot Antenna Test Case
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-27><vspace|0.15fn>>

    <with|par-left|3fn|3.3.9<space|2spc>Dependencies of Dielectric Tissue
    Properties <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-28>>

    <with|par-left|3fn|3.3.10<space|2spc>Probes
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-29>>

    <with|par-left|3fn|3.3.11<space|2spc>Code Parameters
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-30>>

    <with|par-left|1.5fn|3.4<space|2spc>Cell Death Model (RFA & MWA)
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-32>>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|4<space|2spc>Test
    Cases> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-33><vspace|0.5fn>

    <with|par-left|1.5fn|4.1<space|2spc>Solidification due to line heat sink
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-34>>

    <with|par-left|1.5fn|4.2<space|2spc>Empirical SAR
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-35>>

    <with|par-left|1.5fn|4.3<space|2spc>Simple Coaxial Antenna
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-36>>

    <with|par-left|1.5fn|4.4<space|2spc>Cell death test
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-37>>

    <with|par-left|1.5fn|4.5<space|2spc>Coaxial cable
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-38>>

    <with|par-left|1.5fn|4.6<space|2spc>Concentric cylinder test
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-39>>

    <with|par-left|1.5fn|4.7<space|2spc>Concentric cylinder analytic solution
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-40>>

    <with|par-left|1.5fn|4.8<space|2spc>Method of Manufactured solutions
    steady-state Pennes equation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-41>>

    <with|par-left|1.5fn|4.9<space|2spc>Method of Manufactured solutions
    transient nonlinear Pennes equation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-42>>

    <with|par-left|1.5fn|4.10<space|2spc>Transient Laser Heating Benchmark
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-43>>

    <with|par-left|1.5fn|4.11<space|2spc>Single Tine RFA
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-44>>

    <with|par-left|1.5fn|4.12<space|2spc>Single Tine RFA convergence tests
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-45>>

    <with|par-left|1.5fn|4.13<space|2spc>Annular Slot Antenna
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-46>>

    <with|par-left|1.5fn|4.14<space|2spc>Temperature Dependent Electrical
    Properties (MWA) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-47>>

    <with|par-left|1.5fn|4.15<space|2spc>AMICA Prototype 2003
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-48>>

    <with|par-left|1.5fn|4.16<space|2spc>AMICA Prototype 2003 Comparison
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-49>>

    <with|par-left|1.5fn|4.17<space|2spc>AMICA Prototype 2011
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-50>>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|5<space|2spc>Implementation
    Notes> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-51><vspace|0.5fn>

    <with|par-left|1.5fn|5.1<space|2spc>Parameter Dependencies
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-52>>

    <with|par-left|1.5fn|5.2<space|2spc>Nedelec Elements in MWA
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-53>>

    <with|par-left|1.5fn|5.3<space|2spc>Strang Splitting
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-54>>

    <with|par-left|1.5fn|5.4<space|2spc>git Repository Details
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-55>>

    <vspace*|1fn><with|font-series|bold|math-font-series|bold|Bibliography>
    <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
    <no-break><pageref|auto-56><vspace|0.5fn>
  </table-of-contents>

  <new-page>

  Minimally invasive tumour ablation (MITA) therapies involve the placement
  of one or more applicators in close proximity to a tumour under image
  guidance. This can be done either percutaneously or laparoscopically to
  avoid open surgery and the related risks and comorbities. A variety of
  modalities are available: microwave ablation, radiofrequency ablation,
  cryoablation, laser ablation and irreversible electroporation. Once the
  applicator(s) have been positioned they are resonsible for generating or
  removing heat or inducing a pulsed electric field to produce a coagulation
  zone encompassing the tumour.

  The mathematical models describing these MITA therapies can all be solved
  in FEniCS, and due to the similarities in the underlying equations and
  solution methods, it is advantageous to collect them in a single
  application. Models of MITA therapies can be broken into three main
  constituents: a model of the action of the applicator (eg: microwave
  antenna); a model of the heat transfer in the tissue (eq: bioheat
  equation); and a model for cell death as a result of the procedure. These
  three components form a coupled nonlinear system of equations, commonly
  seen in multiphysics applications.

  The use of the bioheat equation is shared by most MITA models, and
  therefore it makes sense to develop a single model, applicable to all, that
  can be reused without having to repeat all of the code verification and
  validation for each modality.

  <section|Getting Started>

  This code relies on the programs git, FEniCS, Python, GMSH and MatLab being
  installed and configured in advance. This is most easily achieved on a
  system running Ubuntu 12.04 LTS, which is the code development platform.
  The code is contained in a git repository (distributed version control
  http://git-scm.com/) and the latest release can be obtained by cloning the
  repository to a local machine using the command line:

  <\pseudo-code>
    cd ``directory''

    git clone https://github.com/sheldonkhall/MITA-model.git
  </pseudo-code>

  <no-indent>where ``directory'' is the desired installation directory. Once
  this has been done the install script should be run by issuing the command:

  <\pseudo-code>
    ./install.sh
  </pseudo-code>

  <no-indent>The install script currently uses GMSH to generate the meshes
  for the test cases and converts them to the FEniCS file format using
  dolfin-convert.

  In order to individually run the test cases provided, the following command
  should be issued on the command line from the installation directory:

  <\pseudo-code>
    python ``filename''.py
  </pseudo-code>

  <no-indent>where ``filename'' should be replaced with one of the test cases
  described in Sec. <reference|sec:test_case>. All tests can be run at once
  using the script run_tests.sh to check an installation.

  <subsection|Mesh Generation>

  All of the test cases include mesh files created using GMSH, an open-source
  finite element grid generator. To identify relevant boundaries and
  subdomains for the solver GMSH can be used to tag certain parts of the
  domain with integers. This is done by defining physical surfaces and lines
  before generating the mesh. Parameters can then be set in the input file to
  inform the solver where to apply boundary conditions and set material
  properties.

  <section|Code Structure>

  The most complex component of the code, utilised in modelling all
  modalities, is the effective heat capacity form of Pennes equation. This is
  a nonlinear PDE solver that is coupled to the other components, e.g. cell
  death model, in specific combinations to form the models of each modality.
  For this reason the core component of this code is the solver for Pennes
  equation, which calls the other components as required. The other
  components can all be called in isolation for static solutions, but in
  order to obtain a transient solution to the full system of equations the
  framework in the Pennes equation solver must be utilised. This is
  illustrated in Fig. <reference|fig:code-schem>.

  <big-figure|<image|code-schem.eps|527px|373px||>|<label|fig:code-schem>Schematic
  of the code highlighting the importance of the Pennes equation solver. The
  information being passed between the components is given explicitly:
  <math|T> is the temperature, <math|V> is the cell viability and SAR the
  specific absorption rate.>

  <section|Code Components>

  In this section the determining equations for each model component are
  detailed for reference. The relevant code parameters are also defined.

  <subsection|Apparent heat capacity form of the Pennes equation>

  In order to compute the evolution of the system over time, including the
  dominant physical processes, the apparent heat capacity form of the Pennes
  equation is solved. This is given by:

  <\equation*>
    \<rho\>c<rsub|app><around*|(|T|)><frac|\<partial\>T|\<partial\>t>=<with|math-font-series|bold|\<nabla\>>\<cdot\>k<around*|(|T|)><with|math-font-series|bold|\<nabla\>>T-\<omega\><rsub|b><around*|(|T|)>
    c<rsub|b><around*|[|T-T<rsub|0>|]>+Q<around*|(|T|)>,
  </equation*>

  where <math|T> is the temperature in K, <math|k > is the thermal
  conductivity, <math|\<omega\><rsub|b>> is the blood perfusion in
  kg/s/m<math|<rsup|3>>, <math|c<rsub|b>> is the specific heat capacity of
  blood in J/kg/K, <math|T<rsub|0>> is the arterial blood temperature taken
  to be 310K and <math|Q> is the heat source generated by the probe in
  W/m<math|<rsup|3>>. The effective specific heat, <math|\<rho\>c<rsub|app>>,
  is given by:

  <\equation*>
    \<rho\>c<rsub|app>=<choice|<tformat|<table|<row|<cell|\<rho\>
    c>|<cell|T\<leq\>T<rsub|l>>>|<row|<cell|<frac|\<rho\> c+\<rho\><rsub|vap>
    c<rsub|vap>|2>+\<rho\><rsub|w> L C <frac|1|T<rsub|u>-T<rsub|l>><rsub|>>|<cell|T<rsub|l>\<less\>T\<leq\>T<rsub|u>>>|<row|<cell|\<rho\><rsub|vap>
    c<rsub|vap>>|<cell|T<rsub|u>\<less\>T>>>>>,
  </equation*>

  where <math|\<rho\>> and <math|\<rho\><rsub|vap>> are the density
  (kg/m<math|<rsup|3>>) of normal and vapourised tissue respectively,
  <math|c> and <math|c<rsub|vap>> are the specific heat (J/kg/K) of normal
  and vapourised tissue respectively, <math|\<rho\><rsub|w>> is the density
  of water, <math|L> is the latent heat of vapourisation (J/kg), <math|C> is
  the tissue water fraction. Water is assumed to be released from the tissue
  over a range of temperatures defined by
  <math|T\<in\><around*|[|T<rsub|l>,T<rsub|u>|]>>, as is expected in
  mixtures.

  To set the various parameters in this model see Tab.
  <reference|tav:pennes>.

  <small-table|<tabular|<tformat|<table|<row|<cell|Symbol>|<cell|Variable
  Name>|<cell|Default Value>>|<row|<cell|<math|\<rho\>
  c>>|<cell|rho_c_t>|<cell|1060 x 3411>>|<row|<cell|<math|\<rho\><rsub|vap>
  c<rsub|vap>>>|<cell|rho_c_v>|<cell|4.4e5>>|<row|<cell|<math|\<rho\><rsub|w>
  L>>|<cell|Lh>|<cell|2260.0e3>>|<row|<cell|<math|C>>|<cell|Cliq>|<cell|0.8>>|<row|<cell|<math|T<rsub|u>>>|<cell|Tu>|<cell|373>>|<row|<cell|<math|T<rsub|l>>>|<cell|Tl>|<cell|363>>|<row|<cell|<math|c<rsub|b>>>|<cell|c>|<cell|3640>>>>>|<label|tav:pennes>>

  <no-indent>All other parameters are dependent on a state variables and are
  treated separately.

  <subsubsection|Thermal conductivity <math|k>>

  The thermal conductivity is set to a constant value thp.k by default, which
  is commonly used in the literature <cite|Trujillo2013> due to the limited
  impact on the results. For the purpose of the sensitivity analysis a linear
  temperature dependence has been introduced, which is selected by setting
  thp.k_method = ``linear''. This sets

  <\equation*>
    k=k<rsub|0>+\<Delta\>k<around*|[|T-T<rsub|0>|]>,
  </equation*>

  where <math|k<rsub|0>> is the baseline thermal conductivity,
  <math|\<Delta\>k> is the change in <math|k > per Kelvin and
  <math|T<rsub|0>> is the reference temperature at which <math|k<rsub|0>> has
  been measured (set to thp.T0 by default).

  When considering phase change the functional form of <math|k> must be valid
  at <math|T\<gtr\>373>K. In order to do this a piecewise continuous function
  is used

  <\equation*>
    k<around*|(|T|)>=<choice|<tformat|<table|<row|<cell|k<rsub|0>+\<Delta\>k<around*|[|T-T<rsub|0>|]>>|<cell|T\<leq\>373K>>|<row|<cell|k<rsub|0>+\<Delta\>k<around*|[|373-T<rsub|0>|]>>|<cell|T\<gtr\>373K>>>>>.
  </equation*>

  This is selected by setting thp.k_method = ``linear_limited''.

  <subsubsection|Steady-state temperature>

  In order to compute a steady-state solution of the Pennes equation, which
  is useful for benchmarking, a few parameter values must be set. These are
  summarised in the table below:

  <small-table|<tabular|<tformat|<table|<row|<cell|Parameter>|<cell|Value>>|<row|<cell|thp.rho_c_t>|<cell|0.>>|<row|<cell|thp.rho_c_v>|<cell|0.>>|<row|<cell|dt>|<cell|a\<gtr\>0.>>|<row|<cell|dt_min>|<cell|a>>|<row|<cell|dt_max>|<cell|a>>|<row|<cell|tmax>|<cell|a>>|<row|<cell|t_out>|<cell|a>>|<row|<cell|thp.perf_model>|<cell|``constant''>>|<row|<cell|thp.k_method>|<cell|``constant''>>|<row|<cell|thp.em_method>|<cell|``custom''>>|<row|<cell|thp.stop_on_me>|<cell|False>>>>>|>

  There is an example of this in the mms-heat-only-test.py file.

  <subsection|Electric Potential Solver (RFA)>

  In order to determine the heat deposited by the probe, a simplified form of
  Maxwell's equations can be solved. Due to the large difference in the
  timescales of the electrical and thermal problems a quasi-static
  approximation can be made and a solution found in terms of the electric
  potential <math|V> determined by:

  <\equation*>
    <with|math-font-series|bold|\<nabla\>>\<cdot\>\<sigma\><around*|(|T|)><with|math-font-series|bold|\<nabla\>>V=0,
  </equation*>

  where <math|\<sigma\>> is the electrical conductivity. This is subject to
  dirichlet boundary conditions on the probe surface and on the external
  surfaces. The probe voltage is set in emp.Vprobe, and the external surface
  voltage in emp.V0 which establishes a potential difference similar to that
  between the probe and ground pads in a patient. This drives a current
  through the tissue resulting in resistive heating. The generated heat can
  be computed using:

  <\equation*>
    SAR=\<sigma\><around*|(|T|)><around*|\||<with|math-font-series|bold|\<nabla\>>V|\|><rsup|2>.
  </equation*>

  <subsubsection|Electrical conductivity <math|\<sigma\>>>

  The electrical conductivity is set to a constant value by default and can
  be adjusted by setting emp.cond. A more realistic linear function of
  temperature is also available:

  <\equation*>
    \<sigma\><around*|(|T|)>=\<sigma\><rsub|0>+\<Delta\>\<sigma\><around*|[|T-T<rsub|0>|]>,
  </equation*>

  where <math|\<sigma\><rsub|0>> is the baseline conductivity,
  <math|\<Delta\>\<sigma\> > is the change in the conduvitvity per Kelvin and
  <math|T<rsub|0>> is the reference temperature at which
  <math|\<sigma\><rsub|0>> has been measured (set to thp.T0 by default). This
  can be selected by setting emp.cond_model = ``linear''.

  When considering phase change the functional form of <math|\<sigma\>> must
  be valid at <math|T\<gtr\>373>K. In order to do this a piecewise continuous
  function is used

  <\equation*>
    \<sigma\><around*|(|T|)>=<choice|<tformat|<table|<row|<cell|\<sigma\><rsub|0>+\<Delta\>\<sigma\><around*|[|T-T<rsub|0>|]>>|<cell|T\<leq\>T<rsub|u>>>|<row|<cell|<around*|[|\<sigma\><rsub|vap>-<around*|[|\<sigma\><rsub|0>+\<Delta\>\<sigma\><around*|[|T<rsub|u>-T<rsub|0>|]>|]>|]><around*|[|T-T<rsub|u>|]>/5+\<sigma\><rsub|0>+\<Delta\>\<sigma\><around*|[|T<rsub|u>-T<rsub|0>|]>>|<cell|T<rsub|u>\<less\>T\<leq\>T<rsub|u>+5>>|<row|<cell|\<sigma\><rsub|vap>>|<cell|T<rsub|u>+5\<less\>T>>>>>,
  </equation*>

  where This is selected by setting thp.k_method = ``nonlinear''.

  <subsubsection|Control System>

  Commercial RFA systems are controlled by adjusting the applied voltage to
  ensure either: the deposited power remains constant; the maximum
  temperature is never exceeded; or the impedance does not cross a threshold
  indicating extensive vapourisation <cite|Haemmerich2003a>.

  <paragraph|Impedance>

  When the impedance specified in emp.imp_max is exceeded the calculation of
  the SAR is stopped for a period of emp.imp_t_off <math|s>
  <cite|Trujillo2013,Trujillo2012>. The flag imp_on reflects whether or not
  the SAR will be computed for time <math|t>.

  <subsubsection|Computation of useful quantities.>

  The impedance of the tissue between the probe and the ground pads is an
  ideal indicator of the presence of vapourisation near the probe surface.
  For this reason it is used to control the power deposition during RFA. An
  impedance threshold is chosen; above which the generator is halted to allow
  the tissue to cool and any generated gas to dissipate from the near
  vicinity of the probe. In the EQS approximation being used to compute the
  potential, the resistance can be computed using:

  <\equation*>
    R<around*|(|t|)>=<frac|V<around*|(|t|)>|P<rsub|total><around*|(|t|)>>,
  </equation*>

  where <math|R<around*|(|t|)>> is the resistance, <math|V<around*|(|t|)>> is
  the applied voltage and <math|P<rsub|total<around*|(|t|)>>> the power
  deposited in the tissue. <math|P<rsub|total><around*|(|t|)>> is computed
  from the SAR by integrating over the computational domain(Kroger):

  <\equation*>
    P<rsub|total><around*|(|t|)>=<big|int>SAR<around*|(|t,\<b-r\>|)><space|1spc>dV.
  </equation*>

  In the axisymmetric cylindrical coordinate system commonly used for these
  computations the equivalent form of the integral is:

  <\equation*>
    P<rsub|total><around*|(|t|)>=2\<pi\><big|int><big|int><big|int>SAR<around*|(|t,\<b-r\>|)><space|1spc>r<space|1spc>dr<space|1spc>dz.
  </equation*>

  <subsection|Vector Helmholtz Solver (MWA)>

  The literature, mostly composed of COMSOL models, indicates that in the
  axisymmetric case, in cylindrical coordinates, Maxwell's Equations simplify
  significantly. The problem reduces to that of finding a scalar field in two
  spacial dimensions <cite|Gas2012>. Unfortunately in 3D there is no such
  simplification and the Vector Helmholtz or curl curl form must be solved.
  The advantage of dealing with a microwave system driven at a single
  frequency with ``simple'' materials however, is that system of equations
  can be solved in the frequency domain. This removes the complications
  related with computing the eigenvalues and eigenmodes.

  <subsubsection|Scalar Form>

  Suppose <math|\<b-H\>> defines the magnetic field over a domain
  <math|\<Omega\>>. We assume <math|\<b-H\><around*|(|\<b-r\>,t|)>=\<b-H\><around*|(|r,z,t|)>>
  is axisymmetric and that <math|\<b-H\><around*|(|\<b-r\>,t|)>=<around*|(|0,<wide|H|~><rsub|\<phi\>>,0|)>>
  in cylindrical coordinates. We further assume that <math|\<b-H\>> is
  time-harmonic, that is we can define an <math|H<rsub|\<phi\>>> such that
  <math|<wide|H|~><rsub|\<phi\>><around*|(|\<b-r\>,t|)>=H<rsub|\<phi\>><around*|(|\<b-r\>|)>*cos<around*|(|\<omega\>*t+\<alpha\>|)>>.
  We take <math|\<alpha\>=0>, without loss of generality. In the frequency
  domain, it is possible to model a lossy medium by considering a
  complex-valued permittivity - <math|\<varepsilon\><rsub|r>=\<varepsilon\><rprime|'><rsub|r>-j\<varepsilon\><rprime|''><rsub|r>>
  - this accounts for the non-local nature of the dispersive model. However,
  this means we must also consider <math|H<rsub|\<phi\>>> as a map to the
  complex domain.

  From Gas <cite|Gas2012> (of the COMSOL family),

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\>*\<times\><around*|[|<around*|(|\<varepsilon\><rsub|r>-j*<frac|\<sigma\>|\<omega\>*\<varepsilon\><rsub|0>>|)><rsup|-1>*\<nabla\>\<times\>\<b-H\>|]>-\<mu\><rsub|r>*k<rsub|0><rsup|2>*\<b-H\>>|<cell|=>|<cell|0<eq-number><label|strong-form>>>>>
  </eqnarray*>

  Then,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\>\<times\>\<b-H\>>|<cell|=>|<cell|<around*|[|<around*|(|\<um\><frac|\<partial\>|\<partial\>
    z>H<rsub|\<phi\>>|)>*\<b-e\><rsub|r>+<around*|(|<frac|1|r>*<frac|\<partial\>|\<partial\>
    r>*<around*|(|r*H<rsub|\<phi\>>|)>|)>*\<b-e\><rsub|z>|]>,<eq-number><label|eqn-curlH>>>>>
  </eqnarray*>

  Using <reference|eqn-curlH> we can expand,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\>\<times\><wide|\<varepsilon\><rsub|r>|^><around*|(|\<nabla\>\<times\>\<b-H\>|)>>|<cell|=>|<cell|<around*|[|<frac|1|r>*<frac|\<partial\>|\<partial\>
    \<phi\>> <around*|(|<frac|<wide|\<varepsilon\><rsub|r>|^>|r>*<frac|\<partial\>|\<partial\>
    r>*<around*|(|r*H<rsub|\<phi\>>|)>|)>*\<b-e\><rsub|r>+<around*|[|<frac|\<partial\>|\<partial\>
    z> <around*|(|\<um\><wide|\<varepsilon\><rsub|r>|^><frac|\<partial\>|\<partial\>
    z>H<rsub|\<phi\>>|)>-<frac|\<partial\>|\<partial\>
    r>*<around*|(|<frac|<wide|\<varepsilon\><rsub|r>|^>|r>*<frac|\<partial\>|\<partial\>
    r>*<around*|(|r*H<rsub|\<phi\>>|)>|)>*\<b-e\><rsub|z>|]>+<frac|1|r>*<frac|\<partial\>|\<partial\>
    \<phi\>> <around*|(|\<um\><wide|\<varepsilon\><rsub|r>|^><frac|\<partial\>|\<partial\>
    z>H<rsub|\<phi\>>|)>*\<b-e\><rsub|z>|]>*\<mathe\><rsup|j \<omega\>
    t>,>>>>
  </eqnarray*>

  where <math|<wide|\<varepsilon\><rsub|r>|^>=<around*|(|\<varepsilon\><rsub|r>-j*<frac|\<sigma\>|\<omega\>*\<varepsilon\><rsub|0>>|)><rsup|-1>>.
  As <math|<frac|\<partial\>|\<partial\> \<phi\>> H<rsub|\<phi\>>\<equiv\>0>,
  (note this is most easily calculated by Lagrange's formula for del, taking
  <math|\<b-H\>=<around*|(|0,H<rsub|\<phi\>>,0|)>>),

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\>\<times\><wide|\<varepsilon\><rsub|r>|^><around*|(|\<nabla\>\<times\>\<b-H\>|)>>|<cell|=>|<cell|\<um\><around*|[|<frac|\<partial\><rsup|>|\<partial\>
    z<rsup|>><around*|[|<wide|\<varepsilon\><rsub|r>|^><frac|\<partial\>|\<partial\>z>H<rsub|\<phi\>>|]>+<frac|\<partial\>|\<partial\>
    r>*<around*|(|<frac|<wide|\<varepsilon\><rsub|r>|^>|r>*<frac|\<partial\>|\<partial\>
    r>*<around*|(|r*H<rsub|\<phi\>>|)>|)>|]>*\<b-e\><rsub|\<phi\>>*\<mathe\><rsup|j
    \<omega\> t>.<eq-number><label|eqn-curlcurlH>>>>>
  </eqnarray*>

  Hence, our governing equation becomes,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|[|<frac|\<partial\><rsup|>|\<partial\>
    z<rsup|>><around*|[|<wide|\<varepsilon\><rsub|r>|^><frac|\<partial\>|\<partial\>z>H<rsub|\<phi\>>|]>+<frac|\<partial\>|\<partial\>
    r>*<around*|(|<frac|<wide|\<varepsilon\><rsub|r>|^>|r>*<frac|\<partial\>|\<partial\>
    r>*<around*|(|r*H<rsub|\<phi\>>|)>|)>|]>*\<b-e\><rsub|\<phi\>>*\<mathe\><rsup|j
    \<omega\> t>+\<mu\><rsub|r>*k<rsub|0><rsup|2>*H<rsub|\<phi\>>*\<b-e\><rsub|\<phi\>>*\<mathe\><rsup|j*\<omega\>*t>>|<cell|=>|<cell|0.>>>>
  </eqnarray*>

  Removing the oscillating term and taking the <math|\<b-e\><rsub|\<phi\>>>
  component,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|[|<frac|\<partial\><rsup|>|\<partial\>
    z<rsup|>><around*|[|<wide|\<varepsilon\><rsub|r>|^><frac|\<partial\>|\<partial\>z>H<rsub|\<phi\>>|]>+<frac|\<partial\>|\<partial\>
    r>*<around*|(|<frac|<wide|\<varepsilon\><rsub|r>|^>|r>*<frac|\<partial\>|\<partial\>
    r>*<around*|(|r*H<rsub|\<phi\>>|)>|)>|]>+\<mu\><rsub|r>*k<rsub|0><rsup|2>*H<rsub|\<phi\>>>|<cell|=>|<cell|0.<eq-number><label|eqn-scalar-form>>>>>
  </eqnarray*>

  This is the scalar form used in the code.

  <subsubsection|Scalar Variational Form>

  Starting from (<reference|eqn-vecvar>), we may introduce the axisymmetric,
  time-harmonic approximation, using a test function
  <math|\<b-T\>=<around*|(|0,T,0|)>>,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|int><rsub|\<Omega\>><around*|[|<wide|\<varepsilon\>|^><rsub|r><rsup|-1>*<around*|(|<frac|\<partial\>|\<partial\>
    z>*T\<cdot\><frac|\<partial\>|\<partial\>
    z>H<rsub|\<phi\>>+<frac|1|r<rsup|2>>*<frac|\<partial\>|\<partial\>
    r>*<around*|(|r*T|)>\<cdot\><frac|\<partial\>|\<partial\>
    r>*<around*|(|r*H<rsub|\<phi\>>|)>|)>|]> d
    \<Omega\>>|<cell|>|<cell|>>|<row|<cell|-<big|int><rsub|\<Omega\>>k<rsub|0><rsup|2>*\<mu\><rsub|r>*T\<cdot\>H<rsub|\<phi\>>
    d \<Omega\>>|<cell|=>|<cell|R*<big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>>T\<cdot\><around*|(|2*H<rsub|\<phi\>
    0>-H<rsub|\<phi\>>|)> d \<Gamma\>,>>>>
  </eqnarray*>

  with the Dirichlet condition <math|H<rsub|\<phi\>>\<equiv\>0> on
  <math|\<Gamma\><rsub|3>>. Rearranging for linear-bilinear forms,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|int><rsub|\<Omega\>><around*|[|<wide|\<varepsilon\>|^><rsub|r><rsup|-1>*<around*|(|<frac|\<partial\>|\<partial\>
    z>*T\<cdot\><frac|\<partial\>|\<partial\>
    z>H<rsub|\<phi\>>+<frac|1|r<rsup|2>>*<frac|\<partial\>|\<partial\>
    r>*<around*|(|r*T|)>\<cdot\><frac|\<partial\>|\<partial\>
    r>*<around*|(|r*H<rsub|\<phi\>>|)>|)>|]> d
    \<Omega\>>|<cell|>|<cell|>>|<row|<cell|-k<rsub|0><rsup|2>*\<mu\><rsub|r><big|int><rsub|\<Omega\>>*T\<cdot\>H<rsub|\<phi\>>
    d \<Omega\>+R*<big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>>T\<cdot\>H<rsub|\<phi\>>
    d \<Gamma\>>|<cell|=>|<cell|2*R*<big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>>T\<cdot\>H<rsub|\<phi\>
    0> d \<Gamma\>.>>>>
  </eqnarray*>

  However, note that <math|T> and <math|H<rsub|\<phi\>>> are complex, and so,
  for the benefit of a real variable solver, we consider the functions as
  maps <math|\<bbb-R\><rsup|2>\<rightarrow\>\<bbb-R\><rsup|2>> to be
  inner-producted. This, with the Dirichlet condition
  <math|H<rsub|\<phi\>><rsup|im>=H<rsub|\<phi\>><rsup|re>=0> on
  <math|\<Gamma\><rsub|3>>, gives us our final equation.

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|int><rsub|\<Omega\>><around*|[|<wide|\<varepsilon\>|^><rsub|r><rsup|-1>*<around*|(|<frac|\<partial\>|\<partial\>
    z>*T\<cdot\><frac|\<partial\>|\<partial\>
    z>H<rsub|\<phi\>>+<frac|1|r<rsup|2>>*<frac|\<partial\>|\<partial\>
    r>*<around*|(|r*T|)>\<cdot\><frac|\<partial\>|\<partial\>
    r>*<around*|(|r*H<rsub|\<phi\>>|)>|)>|]> d
    \<Omega\>>|<cell|>|<cell|>>|<row|<cell|+k<rsub|0><rsup|2>*\<mu\><rsub|r><big|int><rsub|\<Omega\>>*T\<cdot\>H<rsub|\<phi\>>
    d \<Omega\>+R*<big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>>T\<cdot\>H<rsub|\<phi\>>
    d \<Gamma\>>|<cell|=>|<cell|2*R*<big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>>T\<cdot\>H<rsub|\<phi\>
    0> d \<Gamma\>.>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|n\<times\><sqrt|\<varepsilon\>>*E-<sqrt|\<mu\>>H>|<cell|=>|<cell|-2<sqrt|\<mu\>>H<rsub|\<phi\>0>>>|<row|<cell|n\<times\><sqrt|\<varepsilon\><rsub|0>\<varepsilon\><rsub|r>>*E-<sqrt|\<mu\><rsub|0>*\<mu\><rsub|r>>H>|<cell|=>|<cell|-2<sqrt|\<mu\><rsub|0>*\<mu\><rsub|r>>H<rsub|\<phi\>0>>>|<row|<cell|\<mathi\>*\<omega\><sqrt|\<varepsilon\><rsub|0>\<mu\><rsub|0>\<varepsilon\><rsub|r>\<mu\><rsub|r>>
    n\<times\>E-\<mathi\>*\<omega\>*\<mu\><rsub|r>\<mu\><rsub|0>H>|<cell|=>|<cell|>>|<row|<cell|\<mathi\>*k*n\<times\>E+<around*|(|\<nabla\>\<times\>E|)>>|<cell|=>|<cell|\<ldots\>>>>>
  </eqnarray*>

  <subsubsection|Expanded integrals>

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|int><rsub|\<Omega\>><around*|[|<wide|\<varepsilon\>|^><rsub|r><rsup|-1>*<around*|(|<frac|\<partial\>|\<partial\>
    z>*T\<cdot\><frac|\<partial\>|\<partial\>
    z>H<rsub|\<phi\>>+<frac|1|r<rsup|2>>*<frac|\<partial\>|\<partial\>
    r>*<around*|(|r*T|)>\<cdot\><frac|\<partial\>|\<partial\>
    r>*<around*|(|r*H<rsub|\<phi\>>|)>|)>|]> d
    \<Omega\>>|<cell|>|<cell|>>|<row|<cell|+k<rsub|0><rsup|2>*\<mu\><rsub|r><big|int><rsub|\<Omega\>>*T\<cdot\>H<rsub|\<phi\>>
    d \<Omega\>+R*<big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>>T\<cdot\>H<rsub|\<phi\>>
    d \<Gamma\>>|<cell|=>|<cell|2*R*<big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>>T\<cdot\>H<rsub|\<phi\>
    0> d \<Gamma\>.>>>>
  </eqnarray*>

  However, note that <math|T> and <math|H<rsub|\<phi\>>> are complex, and so,
  for the benefit of a real variable solver, we consider the functions as
  maps <math|\<bbb-R\><rsup|2>\<rightarrow\>\<bbb-R\><rsup|2>>, with the
  Dirichlet condition <math|H<rsub|\<phi\>><rsup|im>=H<rsub|\<phi\>><rsup|re>=0>
  on <math|\<Gamma\><rsub|3>> giving us our final equation.

  <subsubsection|Expanding complex-valued functions>

  The assembly of the discretised problem will result in a complex coupling
  matrix and source term. This is due to the property
  <math|H<around*|(|x|)>=a<around*|(|x|)>+j
  b<around*|(|x|)>=<big|sum><rsub|<rsub|n=1>><rsup|\<infty\>>f<rsub|n><around*|(|x|)>a<rsub|n>+j<big|sum><rsub|<rsub|n=1>><rsup|\<infty\>>f<rsub|n><around*|(|x|)>b<rsub|n>=<big|sum><rsub|<rsub|n=1>><rsup|\<infty\>>f<rsub|n><around*|(|x|)><around*|[|a<rsub|n>+j
  b<rsub|n>|]>=<big|sum><rsub|<rsub|n=1>><rsup|\<infty\>>f<rsub|n><around*|(|x|)>h<rsub|n>>
  where <math|f<rsub|n><around*|(|x|)>\<in\>\<bbb-R\>> and
  <math|h<rsub|n>\<in\>\<bbb-C\>>. Solving complex system of equations in
  fenics requires a hack as it is not supported by default <cite|Marais2011>.

  <subsubsection|Vector Variational Form With Inhomogeneous Material
  Properties>

  It may be necessary to develop a full 3D model in the future and this
  section is included for reference, and to show the derivation of the
  boundary conditions. Start by multiplying the strong form in Eq.
  (<reference|strong-form>) by a test function
  <math|\<b-T\>:\<bbb-R\><rsup|3>\<rightarrow\>\<bbb-C\>> without loss of
  generality the test function can be <math|\<b-T\>:\<bbb-R\><rsup|3>\<rightarrow\>\<bbb-R\>>
  according to Section <reference|overall>. The starting point is thus,

  <\equation*>
    <big|int><rsub|\<Omega\>><with|math-font-series|bold|T>\<cdot\><around*|[|<with|math-font-series|bold|\<nabla\>>\<times\><wide|\<epsilon\><rsub|r><with|math-font-series|bold|>|^><around*|(|<with|math-font-series|bold|\<nabla\>>\<times\><with|math-font-series|bold|H>|)>-\<mu\><rsub|r>
    k<rsub|0><rsup|2> <with|math-font-series|bold|H>|]>d\<Omega\>=0.
  </equation*>

  Using the scalar triple product identity for the divergence operator,

  <\equation*>
    <with|math-font-series|bold|\<nabla\>>\<cdot\><around*|(|<with|math-font-series|bold|A>\<times\><with|math-font-series|bold|B>|)>=<with|math-font-series|bold|B>\<cdot\><around*|(|<with|math-font-series|bold|\<nabla\>>\<times\><with|math-font-series|bold|A>|)>-<with|math-font-series|bold|A>\<cdot\><around*|(|<with|math-font-series|bold|\<nabla\>>\<times\><with|math-font-series|bold|B>|)>,
  </equation*>

  one obtains,

  <\equation*>
    <big|int><rsub|\<Omega\>><around*|[|<around*|(|<with|math-font-series|bold|\<nabla\>>\<times\><with|math-font-series|bold|T>|)>\<cdot\><wide|\<epsilon\><rsub|r><with|math-font-series|bold|>|^><around*|(|
    <with|math-font-series|bold|\<nabla\>>\<times\><with|math-font-series|bold|H>|)>-\<mu\><rsub|r>
    k<rsub|0><rsup|2> <with|math-font-series|bold|H>|]>d\<Omega\>=<big|int><rsub|\<Omega\>><with|math-font-series|bold|\<nabla\>>\<cdot\><around*|(|<with|math-font-series|bold|T>\<times\><wide|\<epsilon\><rsub|r><with|math-font-series|bold|>|^><around*|(|
    <with|math-font-series|bold|\<nabla\>>\<times\><with|math-font-series|bold|H>|)>|)>d\<Omega\>.
  </equation*>

  Next apply Gauss' Theorem,

  <\equation*>
    <big|int><rsub|\<Omega\>><with|math-font-series|bold|\<nabla\>>\<cdot\><with|math-font-series|bold|A>
    d\<Omega\> = <big|oint><rsub|\<Gamma\>><with|math-font-series|bold|A>\<cdot\><with|math-font-series|bold|n>
    d\<Gamma\>,
  </equation*>

  to yield a boundary integral (which also appears to be the variational form
  although derived as weak form),

  <\equation*>
    <big|int><rsub|\<Omega\>><around*|[|<around*|(|<with|math-font-series|bold|\<nabla\>>\<times\><with|math-font-series|bold|T>|)>\<cdot\><wide|\<epsilon\><rsub|r><with|math-font-series|bold|>|^><around*|(|
    <with|math-font-series|bold|\<nabla\>>\<times\><with|math-font-series|bold|H>|)>-\<mu\><rsub|r>
    k<rsub|0><rsup|2> <with|math-font-series|bold|H>|]>d\<Omega\>=<big|oint><rsub|\<Gamma\>><around*|(|<with|math-font-series|bold|T>\<times\><wide|\<epsilon\><rsub|r><with|math-font-series|bold|>|^><around*|(|<with|math-font-series|bold|\<nabla\>>\<times\><with|math-font-series|bold|H>|)>|)>\<cdot\><with|math-font-series|bold|n>
    d\<Gamma\>.
  </equation*>

  The surface <math|\<Gamma\>> enclosing the volume <math|\<Omega\>> is
  formed from the union of four types of boundary condition in this case
  <math|\<Gamma\>=\<Gamma\><rsub|1><rsub|>\<cup\>\<Gamma\><rsub|2>\<cup\>\<Gamma\><rsub|3>\<cup\>\<Gamma\><rsub|4>>.
  Let <math|\<Gamma\><rsub|1>> be the scattering boundary,
  <math|\<Gamma\><rsub|2>> be the conductor boundary,
  <math|\<Gamma\><rsub|3>> be the axisymmetric boundary and
  <math|\<Gamma\><rsub|4>> be the feed port.

  <subsubsection|First-Order Absorbing Boundary Condition (Scattering)>

  In order to solve the computational electromagnetic problem in a finite
  domain an approximation of the Sommerfeld radiation condition is required
  on the external boundary. A first-order absorbing boundary condition should
  be sufficient for this application, but if further accuracy is required a
  perfectly matched layer could be used <cite|Jin2009>. The distance at which
  to apply any boundary condition can be determined using a sensitivity
  analysis. The decay of the microwaves as they propagate through tissue is
  significant, which may relax the requirements on the accuracy of the
  boundary condition.

  According to Jin <cite|Jin2009> the first-order absorbing boundary
  condition is given by:

  <\equation*>
    n\<times\>\<nabla\>\<times\>\<b-H\>+j k<rsub|0>
    n\<times\>n\<times\>\<b-H\>\<approx\>0,<space|1em>\<b-r\> \<epsilon\>
    \<Gamma\><rsub|1>.
  </equation*>

  This can be substituted into the boundary integral and simplified,

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|int><rsub|\<Gamma\><rsub|1>>\<b-T\>\<times\><wide|\<epsilon\><rsub|r>|^><around*|(|\<nabla\>\<times\>\<b-H\>|)>\<cdot\>n
    d\<Gamma\>>|<cell|=>|<cell|<big|int><rsub|\<Gamma\><rsub|1>>\<b-T\>\<times\><wide|\<epsilon\><rsub|r>|^><around*|(|-j
    k<rsub|0><around*|(|n\<times\>\<b-H\>|)>|)>\<cdot\>n
    d\<Gamma\>,>>|<row|<cell|>|<cell|=>|<cell|<big|int><rsub|\<Gamma\><rsub|1>><wide|\<epsilon\><rsub|r>|^>
    j k<rsub|0><around*|(|n\<times\>\<b-H\>|)>\<times\>\<b-T\>\<cdot\>n
    d\<Gamma\>,>>|<row|<cell|>|<cell|=>|<cell|<big|int><rsub|\<Gamma\><rsub|1>><wide|\<epsilon\><rsub|r>|^>
    j k<rsub|0><around*|(|n\<times\>\<b-H\>|)>\<cdot\><around*|(|\<b-T\>\<times\>n|)>
    d\<Gamma\>,>>|<row|<cell|>|<cell|=>|<cell|-<big|int><rsub|\<Gamma\><rsub|1>><wide|\<epsilon\><rsub|r>|^>
    j k<rsub|0> <around*|(|n\<times\>\<b-H\>|)>\<cdot\><around*|(|n\<times\>\<b-T\>|)>
    d\<Gamma\>,>>>>
  </eqnarray*>

  using the vector product identities,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-A\>\<cdot\><around*|(|\<b-B\>\<times\>\<b-C\>|)>>|<cell|=>|<cell|\<b-B\>\<cdot\><around*|(|\<b-C\>\<times\>\<b-A\>|)>>>|<row|<cell|>|<cell|<with|font-family|rm|and>>|<cell|>>|<row|<cell|\<b-A\>\<times\>\<b-B\>>|<cell|=>|<cell|-\<b-B\>\<times\>\<b-A\>,>>>>
  </eqnarray*>

  which agrees with Jin <cite|Jin2009>.

  <subsubsection|Waveguide Port Boundary Condition>

  When modelling antennas it is sometimes necessary to model the waveguide to
  obtain accurate solutions <cite|Jin2009>. It is not necessary to model the
  whole waveguide, just the region in the vicinity of the antenna to capture
  the modified modal distribution. This truncation of the feeding waveguide
  to form a waveguide port requires a special boundary condition. This port
  boundary condition is required to absor<em|>b the field reflected back into
  the waveguide and accurately represent the modal distribution in the
  waveguide region forming the antenna.

  The following WPBC can be used <cite|Gas2012>

  <\equation*>
    <sqrt|\<epsilon\>-j<frac|\<sigma\>|\<omega\>>>
    <with|math-font-series|bold|n>\<times\><with|math-font-series|bold|E>-<sqrt|\<mu\>><with|math-font-series|bold|H<rsub|\<phi\>>>=-2<sqrt|\<mu\>><with|math-font-series|bold|H<rsub|\<phi\>0>>,
  </equation*>

  where <math|<with|math-font-series|bold|H<rsub|\<phi\>0>>=<around*|(|0,H<rsub|\<phi\>0>,0|)><rsup|T>>
  and <math|<with|math-font-series|bold|H<rsub|\<phi\>>>=<around*|(|0,H<rsub|\<phi\>>,0|)><rsup|T>>.
  Further, note that

  <\equation*>
    j<frac|\<omega\>\<epsilon\><rsub|o>|<wide|\<epsilon\><rsub|r>|^>><with|math-font-series|bold|E>=<with|math-font-series|bold|\<nabla\>>\<times\><with|math-font-series|bold|H>.
  </equation*>

  Considering the boundary term in the variational form we can write

  <\eqnarray*>
    <tformat|<table|<row|<cell|<big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>>\<b-T\>\<times\><wide|\<epsilon\><rsub|r>|^><around*|(|\<nabla\>\<times\>\<b-H\>|)>\<cdot\><with|math-font-series|bold|n>
    d\<Gamma\>>|<cell|=>|<cell|<big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>>\<b-T\>\<cdot\><around*|(|<wide|\<epsilon\><rsub|r>|^><around*|(|\<nabla\>\<times\>\<b-H\>|)>\<times\><with|math-font-series|bold|n>|)>
    d\<Gamma\>,>>|<row|<cell|>|<cell|=>|<cell|-<big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>>\<b-T\>\<cdot\><around*|(|<with|math-font-series|bold|n>\<times\><wide|\<epsilon\><rsub|r>|^><around*|(|\<nabla\>\<times\>\<b-H\>|)>|)>
    d\<Gamma\>,>>|<row|<cell|\<omega\>\<epsilon\><rsub|o><space|1em><with|font-family|rm|assumed
    constant and factors out of curl>>|<cell|=>|<cell|-j\<omega\>\<epsilon\><rsub|o><big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>>\<b-T\>\<cdot\><around*|(|<with|math-font-series|bold|n>\<times\><with|math-font-series|bold|E>|)>
    d\<Gamma\>.>>>>
  </eqnarray*>

  Now the boundary condition can be rewritten slightly as

  <\equation*>
    <frac|<sqrt|\<epsilon\><rsub|o>>|<sqrt|<wide|\<epsilon\><rsub|r>|^>>>
    <with|math-font-series|bold|n>\<times\><with|math-font-series|bold|E>-<sqrt|\<mu\>><with|math-font-series|bold|H<rsub|\<phi\>>>=-2<sqrt|\<mu\>><with|math-font-series|bold|H<rsub|\<phi\>0>>.
  </equation*>

  Substituting this into the boundary term

  <\eqnarray*>
    <tformat|<table|<row|<cell|-j\<omega\>\<epsilon\><rsub|o><big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>>\<b-T\>\<cdot\><around*|(|<with|math-font-series|bold|n>\<times\><with|math-font-series|bold|E>|)>
    d\<Gamma\>>|<cell|=>|<cell|-j\<omega\><sqrt|\<mu\><rsub|o>><sqrt|\<epsilon\><rsub|o>><big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>><sqrt|\<mu\><rsub|r>><sqrt|<wide|\<epsilon\><rsub|r>|^>>\<b-T\>\<cdot\><around*|(|<with|math-font-series|bold|H<rsub|\<phi\>>>-2<with|math-font-series|bold|H<rsub|\<phi\>0>>|)>
    d\<Gamma\>>>|<row|<cell|>|<cell|=>|<cell|-j
    k<rsub|o><big|int><rsub|\<Gamma\><rsub|1>\<cup\>\<Gamma\><rsub|4>><sqrt|\<mu\><rsub|r>><sqrt|<wide|\<epsilon\><rsub|r>|^>>\<b-T\>\<cdot\><around*|(|<with|math-font-series|bold|H<rsub|\<phi\>>>-2<with|math-font-series|bold|H<rsub|\<phi\>0>>|)>
    d\<Gamma\>>>>>
  </eqnarray*>

  <subsubsection|Verification Cases (Axisymmetric)>

  In this section the theory behind some of the MWA test cases is described.

  <paragraph|Coaxial cable analytic solution>

  Let <math|H<rsub|\<phi\>>=<frac|C|r*Z>*\<mathe\><rsup|-\<mathi\>*k*z>>.
  Then,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\>\<times\><around*|(|H<rsub|\<phi\>>*\<b-e\><rsub|\<phi\>>|)>>|<cell|=>|<cell|<around*|(|-<frac|\<partial\>
    H<rsub|\<phi\>>|\<partial\> z>*\<b-e\><rsub|r>+<frac|1|r>*<frac|\<partial\>|\<partial\>
    r><around*|(|r*H<rsub|\<phi\>>|)>*\<b-e\><rsub|z>|)>>>|<row|<cell|>|<cell|=>|<cell|<frac|C|Z>*<around*|(|\<mathi\>*k*<frac|1|r>\<mathe\><rsup|\<um\>\<mathi\>*k*z>*\<b-e\><rsub|r>|)>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\>\<times\><around*|(|\<nabla\>\<times\><around*|(|H<rsub|\<phi\>>*\<b-e\><rsub|\<phi\>>|)>|)>>|<cell|=>|<cell|<frac|C|Z>*<around*|(|k<rsup|2>*<frac|1|r>*\<mathe\><rsup|\<um\>*\<mathi\>*k*z>|)>*\<b-e\><rsub|\<phi\>>>>|<row|<cell|>|<cell|=>|<cell|k<rsup|2>*H<rsub|\<phi\>>*\<b-e\><rsub|\<phi\>>>>|<row|<cell|>|<cell|=>|<cell|k<rsub|0><rsup|2>*\<varepsilon\><rsub|r>*\<mu\><rsub|r>*H<rsub|\<phi\>>*\<b-e\><rsub|\<phi\>>>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\>\<times\><around*|(|\<varepsilon\><rsub|r><rsup|-1>*\<nabla\>\<times\><around*|(|H<rsub|\<phi\>>*\<b-e\><rsub|\<phi\>>|)>|)>>|<cell|=>|<cell|k<rsub|0><rsup|2>*\<mu\><rsub|r>*H<rsub|\<phi\>>*\<b-e\><rsub|\<phi\>>.>>>>
  </eqnarray*>

  Trying the boundary condition,

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\>\<times\>H<rsub|\<phi\>>>|<cell|=>|<cell|<frac|C|Z>*<around*|(|\<mathi\>*k*<frac|1|r>*\<mathe\><rsup|\<um\>\<mathi\>*k*z>*\<b-e\><rsub|r>|)>>>|<row|<cell|>|<cell|=>|<cell|\<omega\>*<around*|(|\<mathi\>*\<varepsilon\>*<frac|1|r>*\<mathe\><rsup|\<um\>\<mathi\>*k*z>*\<b-e\><rsub|r>|)>>>|<row|<cell|\<Rightarrow\>>|<cell|>|<cell|>>|<row|<cell|E>|<cell|=>|<cell|<frac|1|\<mathi\>*\<omega\>*\<varepsilon\>>*\<nabla\>\<times\>H<rsub|\<phi\>>>>|<row|<cell|>|<cell|=>|<cell|C*<around*|(|<frac|1|r>\<mathe\><rsup|\<um\>\<mathi\>*k*z>*\<b-e\><rsub|r>|)>>>|<row|<cell|n\<times\>E>|<cell|=>|<cell|<around*|(|n<rsub|r>*\<b-e\><rsub|r>+n<rsub|z>*\<b-e\><rsub|z>|)>\<times\>E>>|<row|<cell|>|<cell|=>|<cell|n<rsub|z>*<frac|C|r>*\<mathe\><rsup|\<um\>\<mathi\>*k*z>*\<b-e\><rsub|\<phi\>>>>|<row|<cell|>|<cell|=>|<cell|n<rsub|z>*<frac|<sqrt|\<mu\>|>|<sqrt|\<varepsilon\>>>*<around*|(|<frac|C|r*Z>*\<mathe\><rsup|-\<mathi\>*k*z>*\<b-e\><rsub|\<phi\>>|)>>>|<row|<cell|>|<cell|=>|<cell|n<rsub|z>*<frac|<sqrt|\<mu\>|>|<sqrt|\<varepsilon\>>>*H<rsub|\<phi\>>*\<b-e\><rsub|\<phi\>>>>|<row|<cell|\<Rightarrow\>>|<cell|>|<cell|>>|<row|<cell|<sqrt|\<varepsilon\>>*n\<times\>E>|<cell|=>|<cell|n<rsub|z>*<sqrt|\<mu\>|>*H<rsub|\<phi\>>*\<b-e\><rsub|\<phi\>>.>>>>
  </eqnarray*>

  To relate this to the WPBC already implemented realise from the COMSOL
  manual that

  <\equation*>
    P<rsub|in>=<frac|\<pi\>C<rsup|2>|Z>ln<around*|(|r<rsub|1>/r<rsub|2>|)>\<rightarrow\>C=<sqrt|<frac|P<rsub|in>Z|\<pi\>ln<around*|(|r<rsub|1>/r<rsub|2>|)>>>,
  </equation*>

  this allows us to test the WPBC we want to use. It appears that the
  first-order absorbing boundary condition is not sufficient for this
  problem, producing around 10% errors. If however, the WPBC is used with
  <math|H<rsub|\<phi\>0>\<equiv\>0> then the error becomes negligible. This
  makes sense as the WPBC is setting the modes of the waveguide to <math|0>
  whereas the first-order absorbing boundary condition is approximate.
  Conversely using the WPBC on the conc-cylinder test gives worse solutions.

  <paragraph|Concentric Cylinders>

  This case consists of two concentric cylinders of biological tissues, which
  can have different electrical properties and sizes, with an imposed
  z-directional electric field. In this simple geometry an analytic solution
  can be derived, which is accepted without confirmation. The solution is

  <\equation*>
    E<rsub|z><around*|(|r|)>=A<rsub|i>J<rsub|o><around*|(|k<rsub|i>r|)>+B<rsub|i>Y<rsub|o><around*|(|k<rsub|i>r|)>
  </equation*>

  where <math|J<rsub|o>> is the Bessel function of the first kind,
  <math|Y<rsub|o>> is the Bessel function of the second kind,
  <math|k=\<omega\><sqrt|\<mu\>\<epsilon\><rsub|o><rsup|>/<wide|\<epsilon\><rsub|r>|^><rsup|*>>>,
  <math|A> and <math|B> are determined from the boundary conditions and the
  subscript <math|i> corresponds to the material properties within a
  particular cylinder. To obtain <math|A> and <math|B> the following system
  of equations must be solved

  <paragraph|Annular Slot Antenna Test Case>

  This test case is based on the paper by Nevels <cite|Nevels1985>. It is a
  coaxially fed annular aperture antenna propagating into a lossy biological
  medium. The dimensions are given by products of <math|k> and the radii.
  This is a test of the full 2D axisymmetric domain, not a pseudo 2D case.
  Further it is also a test of: the WPBC, first order absorbing boundary
  condition, perfect electrical conductor boundary condition, symmetry
  condition and complex dielectric material properties. Note that the imposed
  TEM mode in the coaxial cable is defined as

  <\equation*>
    E<rsub|r>=<frac|V|ln<around*|(|b/a|)>><space|1spc><frac|1|r><space|1spc>e<rsup|-jkz>,
  </equation*>

  where <math|k=\<omega\><sqrt|\<mu\>\<epsilon\>>>, <math|V> is the applied
  voltage, <math|a> and <math|b> are the inner and outer radii of the coaxial
  respectively. In order to apply an appropriate WPBC, \ and ensure that the
  solutions correspond with those in the paper, the WPBC must be position at
  <math|z=n2\<pi\>/k>, where <math|n> is a posotive integer. This places the
  boundary condition at multiples of the wavelength away from the interface
  with the external lossy domain. Converting this to <math|H<rsub|\<phi\>>>
  gives us

  <\equation*>
    <with|math-font-series|bold|<wide|\<phi\>|^>><space|1spc>H<rsub|\<phi\>>=-<space|1spc><frac|1|j\<omega\>\<mu\>><around*|(|<with|math-font-series|bold|\<nabla\>>\<times\><with|math-font-series|bold|<wide|r|^>>|)><space|1spc>E<rsub|r>=<frac|k
    V|\<omega\> \<mu\> ln<around*|(|b/a|)>><space|1spc><frac|1|r><space|1spc>e<rsup|-jkz>.
  </equation*>

  Relating terms to the <math|H<rsub|\<phi\>0>> incident field in the Gas
  paper <cite|Gas2012> and the implemented WPBC we can show that

  <\equation*>
    H<rsub|\<phi\>>=<frac|V|Z ln<around*|(|b/a|)>><space|1spc><frac|1|r><space|1spc>e<rsup|-jkz>\<rightarrow\>C=<frac|V|ln<around*|(|b/a|)>>
  </equation*>

  <subsubsection|Dependencies of Dielectric Tissue Properties>

  In the general case we assume that the dielectric properties are given by

  <\equation*>
    \<epsilon\><rsub|r><around*|(|F,DI,T|)><text| and
    >\<sigma\><around*|(|F,DI,T|)>,
  </equation*>

  where <math|F> is the water content, <math|DI> is the heat damage integral,
  <math|T> is the temperature. The actual values of these parameters, with
  varying temperature, have been measured by Lapresto et al. & Ji and Brace.
  These measurements will implicitly capture <math|DI> and <math|F> as they
  were taken during MWA. The functional form used to describe the temperature
  dependence is

  <\equation*>
    \<eta\><around*|(|T|)>=s<rsub|1><around*|[|1-<frac|1|1+exp<around*|(|s<rsub|2>-s<rsub|3>T|)>>|]>+s<rsub|4>,
  </equation*>

  where <math|\<eta\>\<in\><around*|{|\<epsilon\><rsub|r>,\<sigma\>|}>> and
  <math|s<rsub|n>> are parameters. It should be noted that
  <math|s<rsub|4>\<in\><around*|{|0,1|}>> is not a true parameter, but in
  fact represents a different functional form for <math|\<epsilon\><rsub|r>>
  used by Ji and Brace. For this reason only parameter values for
  <math|\<sigma\>> are probably directly comparable. The quoted values are
  given in the table

  <tabular|<tformat|<table|<row|<cell|>|<cell|<math|s<rsub|1>>>|<cell|<math|s<rsub|2>>>|<cell|<math|s<rsub|3>>>|<cell|<math|s<rsub|4>>>>|<row|<cell|<math|\<epsilon\><rsub|r>>
  - Ji & Brace>|<cell|48.391>|<cell|6.286>|<cell|0.0764>|<cell|1>>|<row|<cell|<math|\<epsilon\><rsub|r>>
  - Lapresto>|<cell|45.57>|<cell|5.223>|<cell|0.05243>|<cell|0>>|<row|<cell|<math|\<sigma\>>
  - Ji & Brace>|<cell|2.173>|<cell|5.951>|<cell|0.0697>|<cell|0>>|<row|<cell|<math|\<sigma\>>
  - Lapresto>|<cell|1.803>|<cell|6.538>|<cell|0.05984>|<cell|0>>>>>

  <subsubsection|Probes>

  The shape of the SAR is highly dependent on the structure of the applicator
  (microwave antenna) used. For this reason it is necessary to have a library
  of applicator models, to match those commercially available, even in
  axisymmetric coordinates. Currently these are implemented in test cases and
  includes:

  <\enumerate-numeric>
    <item>Annular slot antenna (Gas),

    <item>Coaxial monopole (Ji and Brace),

    <item>AMICA (prototype 1),

    <item>AMICA (prototype 2).
  </enumerate-numeric>

  <subsubsection|Code Parameters>

  The relevant code variables for the equations derived in this section are
  given in Tab. <reference|tab:mwa>. There are examples of setting up a
  static computation of the SAR for MWA in tosoratti-2003.py and
  tosoratti-2011.py.

  <small-table|<tabular|<tformat|<table|<row|<cell|Symbol>|<cell|Variable
  Name>|<cell|Description>>|<row|<cell|<math|\<omega\>>>|<cell|om>|<cell|angular
  frequency>>|<row|<cell|<math|P<rsub|in>>>|<cell|Pin>|<cell|applied
  applicator power>>|<row|<cell|<math|r<rsub|1>>>|<cell|r_1>|<cell|inner
  radius coax feed line>>|<row|<cell|<math|r<rsub|2>>>|<cell|r_2>|<cell|outer
  radius coax feed line>>|<row|<cell|Re(<math|H<rsub|\<phi\>0>>)>|<cell|H_phi_0_re>|<cell|real
  component of the coax fundamental mode>>|<row|<cell|Im(<math|H<rsub|\<phi\>0>>)>|<cell|H_phi_0_im>|<cell|imaginary
  component of the coax fundamental mode>>|<row|<cell|<math|\<epsilon\><rsub|r>>>|<cell|eps_r_by_subdomain>|<cell|relative
  permittivity>>|<row|<cell|<math|\<sigma\>>>|<cell|sigma_by_subdomain>|<cell|effective
  conductivity>>|<row|<cell|<math|\<mu\><rsub|r>>>|<cell|mu_r_by_subdomain>|<cell|relative
  permeability>>|<row|<cell|<math|C>>|<cell|C_dielectric>|<cell|parameter of
  waveguide port boundary condition>>|<row|<cell|<math|Z>>|<cell|Z_dielectric>|<cell|parameter
  of waveguide port boundary condition>>|<row|<cell|<math|s<rsub|1>>>|<cell|es1
  / ss1>|<cell|parameter of sigmoidal form of <math|\<epsilon\><rsub|r>> /
  <math|\<sigma\>>>>|<row|<cell|<math|s<rsub|2>>>|<cell|es2 /
  ss2>|<cell|parameter of sigmoidal form of <math|\<epsilon\><rsub|r>> /
  <math|\<sigma\>>>>|<row|<cell|<math|s<rsub|3>>>|<cell|es3 /
  ss3>|<cell|parameter of sigmoidal form of <math|\<epsilon\><rsub|r>> /
  <math|\<sigma\>>>>|<row|<cell|<math|s<rsub|4>>>|<cell|es4 /
  ss4>|<cell|parameter of sigmoidal form of <math|\<epsilon\><rsub|r>> /
  <math|\<sigma\>>>>>>>|<label|tab:mwa>>

  <subsection|Cell Death Model (RFA & MWA)>

  The cell death model that has been implemented is that of O'Neill et al.
  The parameters are set, by default, to the optimal values given in Tab. 2
  of the paper. The cell death model is solved for every nodal point in the
  finite element mesh, over the current time step, to produce a scalar field
  representing the cell viability vs time. To turn on the cell death model a
  flag must be set cda_update = True.

  In order to stop the perfusion during a calculation, as a result of damage
  to the vasculature, the parameter perf_model = ``stop'' must be set. This
  corresponds to a viability <math|\<less\> 0.8>.

  <section|Test Cases><label|sec:test_case>

  In order to perform code verification, a collection of test cases has been
  assembled that will eventually cover all implemented features. This is a
  collection of problems with analytic solutions, manufactured solutions and
  results from the literature.

  <subsection|Solidification due to line heat sink>

  File: 1d-solidification-test.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: this test case sets up a line heat sink at
  <math|r=0> in a domain initially at <math|310>K and then models the phase
  change as a function of time. An analytic solution is available for this
  problem, which is computed and stored in the file analytic.pvd.

  \;

  <no-indent>Expected Result: the predicted temperature profile in
  enthalpy.pvd should agree closely with that in analytic.pvd. Further,
  numerical predictions of the interface location (defined as
  Tu\<gtr\>T\<gtr\>Tl) \ should coincide with the analytic solution at all
  times.

  <subsection|Empirical SAR>

  File: ai-2012.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: The experiments performed by Haiming Ai et al.
  <cite|Ai2012> are partially implemented to provide a comparison against
  experimental results.

  \;

  <no-indent>Expected Results: Reproduces the results in the publication.

  <subsection|Simple Coaxial Antenna>

  File: annular-slot-antenna.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: This is a comparison against a different method of
  computation. A coaxial cable is cut perpendicular to the direction of wave
  propagation and the end brought into contact with a lossy biological tissue
  <cite|Nevels1985>.\ 

  \;

  <no-indent>Expected Results: Reproduces the results in the publication.

  <subsection|Cell death test>

  File: cell-death-test.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: this test case runs a transient calculation with a
  fixed homogeneous temperature. This is to allow comparison of the
  implemented cell death ODE solver with a standard adaptive ODE solver.

  \;

  <no-indent>Expected Result: V should agree with the solutions generated by
  cell_death_test.m in the matlab directory. The matlab results are stored in
  matlab.pvd for easy comparison in paraview.

  <subsection|Coaxial cable>

  File: coaxial-test.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: this test case computes the <math|E> and <math|H>
  field in a piece of coaxial cable for which analytic solutions are
  available.

  \;

  <no-indent>Expected Result: should broadly agree with the analytic
  solutions. The agreement is not the best because the absorbing boundary
  condition is not accurate enought for this application. A much more
  accurate boundary condition is commented out in the code, but of no use for
  modelling realistic MITA cases.

  <subsection|Concentric cylinder test>

  File: conc-cyl-test.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: this test case imposes an electric field on the
  outer surface of two concentric cylinders of tissue. An analytic solution
  exists for this situation and the comparison is performed and results
  displayed graphically and on the command line.

  \;

  <no-indent>Expected Result: The L2 norm error should be low and the second
  plot should be close to zero in the middle of the domain (not throughout
  due to boundary conditions at cylinder ends).

  <subsection|Concentric cylinder analytic solution>

  File: conc-cyl-analy-sol.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: implementation of the analytic solution in FEniCS.

  \;

  <no-indent>Expected Result: produces the analytic solution on computational
  mesh.

  <subsection|Method of Manufactured solutions steady-state Pennes equation>

  File: mms-heat-only-test.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: this test case uses the method of manufactured
  solutions to test that steady state solutions to the Pennes equation can be
  successfully computed.

  \;

  <no-indent>Expected Results: The L2 norm error should be low, the error is
  plotted in T-error.pvd.

  <subsection|Method of Manufactured solutions transient nonlinear Pennes
  equation>

  File: mms-heat-only-nonlinear.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: this test case uses the method of manufactured
  solutions to test the transient nonlinear capabilities of fenics_mwa.

  \;

  <no-indent>Expected Results: The L2 norm error should be low, the
  temperature profiles should agree at all times.

  <subsection|Transient Laser Heating Benchmark>

  File: vyas-1992.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: this test case involves a laser source in tissue to
  provide heat that is only active for an initial period and then switched
  off <cite|Vyas1992>. The source is provided as an SAR, and therefore this
  tests the transient linear Pennes equations solution. There is an analytic
  solution and this is used for comparison.

  \;

  <no-indent>Expected Results: The publication uses a cross section of the
  temperature at a specific time point to evaluate the solution accuracy.

  <subsection|Single Tine RFA>

  File: single-tine-sym-rfa.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: this test case is going to be the base of the RFA
  sensitivity study. For this reason it has all of the nonlinearities present
  and the impedance based control system switched on. It is based on the test
  case used by Trujillo et al. cite{Trujillo2012} which has some example
  solutions and experimental data for comparison.

  \;

  <no-indent>Expected Results: Currently expect code to converge to a
  solution with physically reasonable temeprature and cell death predictions.

  <subsection|Single Tine RFA convergence tests>

  File: single-tine-sym-rfa-conv.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: this test runs through a variety of meshes of
  increasing resolution to determine the mesh sensitivity. The total mesh
  volume is also increased to see if there is a dependence on this.

  \;

  <no-indent>Expected Results: array of convergence of max and min field
  quantities.

  <subsection|Annular Slot Antenna>

  File: gas-probe-sensitivity.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: this test case implements the MWA applicator
  specified in Gas et al.

  \;

  <no-indent>Expected Results: -

  <subsection|Temperature Dependent Electrical Properties (MWA)>

  File: ji-brace-2011.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: this test case evaluates temperature dependent
  electrical properties specified in Ji & Brace cite{Ji2011}. This problem is
  nonlinear and therefore needs small timesteps. Very clear difference can be
  seen between a constant value and temperature dependence.

  \;

  <no-indent>Expected Results: Should agree with results in paper.

  <subsection|AMICA Prototype 2003>

  File: tosoratti-2003.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: computes the SAR for the AMICA prototype.

  \;

  <no-indent>Expected Results: Should agree with results in paper.

  <subsection|AMICA Prototype 2003 Comparison>

  File: tosoratti-2003-unchoked.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: computes the SAR for the comparison unchoked
  applicator to show the difference.

  \;

  <no-indent>Expected Results: Should agree with results in paper.

  <subsection|AMICA Prototype 2011>

  File: tosoratti-probe-2011.py

  \;

  <no-indent>Status: working

  \;

  <no-indent>Description: computes the SAR for the AMICA prototype.

  \;

  <no-indent>Expected Results: Should agree with results in paper.

  <section|Implementation Notes>

  Various implementation issues inevitably arise when developing a code of
  this type. These are detailed here for reference.

  <subsection|Parameter Dependencies>

  There are a variety of functional forms for the dependencies of the model
  parameters proposed in the literature. Many of these vary spatially through
  their dependence on temperature or cell viability and can pose a
  computational issue as a result of the polynomial basis in FEM. The
  dramatic drop in conductivity as the temperature approaches 373K can occur
  in the middle of an element leading to negative values of the interpolating
  polynomials.

  For this reason the state variable e.g. temperature is projected onto a
  discontinuous galerkin (DG) order 0 function space before evaluating the
  dependent parameter e.g. <math|\<sigma\>>. This ensures that no negative
  values are erroneously generated. It also ensures that the parameter values
  are piecewise constant, which is compatible with FEM. The fine mesh in the
  areas where these effects are observed along with the size of the inherent
  uncertainties in the parameter values mean that this approach should not
  impact the solution accuracy greatly.

  <subsection|Nedelec Elements in MWA>

  Lagrange finite elements are not appropriate for the solution of Maxwell's
  Equations due to the interference as a result of spurious modes. To address
  this Nedelec elements need to be used in their place. These are a standard
  option in FEniCS.

  <subsection|Strang Splitting>

  The final system of equations for RFA and MWA includes PDEs describing the
  transfer of heat and the electric field and a system of ODEs describing the
  cell death. In order to solve this system of nonlinear equations an
  <with|font-shape|italic|ad hoc> linearisation has been performed. It should
  be possible to improve this in the future using operator splitting, which
  can generate iterative schemes for this type of problem. This type of
  method is commonly employed to solve the bidomain equations in cardiac
  modelling.

  <subsection|git Repository Details>

  The branching model used is available here:
  http://nvie.com/posts/a-successful-git-branching-model/ (19/09/14). The
  develop branch will always be kept on a local machine and will be the most
  complete version containing files not available on the master branch. The
  master branch will be published online on github to allow users to access
  working copies and to start making their own modifications. If the need
  arises, the develop branch can be posted online for collaboration. The url
  for the master branch on github is:

  <\pseudo-code>
    https://github.com/sheldonkhall/MITA-model.git
  </pseudo-code>

  <no-indent>The release branch will normally contain a subset of the develop
  branch and therefore care must be taken when performing a merge.
  Essentially, when merging a bugfix performed in the feature branch into the
  develop branch the ours merge strategy must be selected. The command is:

  <\pseudo-code>
    git checkout develop

    git merge --no-f -s ours release
  </pseudo-code>

  <no-indent>This should stop files removed for the release being removed
  from the develop branch.

  When preparing a release the command:

  <\pseudo-code>
    git clean
  </pseudo-code>

  <no-indent>is useful for removing untracked (automatically generated)
  files. This will NOT remove .pyc files however and these may cause problems
  when not doing a fresh install.

  In order to list all deleted files (useful when merging changes from hotfix
  and release back to develop use:

  <\pseudo-code>
    git log --diff-filter=D --summary
  </pseudo-code>

  In order to keep the online repository clean: a new branch was created.
  This is called trunk and will contain the released publically available
  code. This is stripped of the history during creating by using the command:

  <\pseudo-code>
    git checkout --orphan trunk
  </pseudo-code>

  <no-indent>then trunk is merged with release. This branch can then be
  pushed to github and lacks the unneeded develop history.

  When updating for release 1.1 and 2.0 the commit history was removed using:

  <\pseudo-code>
    git merge --squash
  </pseudo-code>

  <no-indent>After release 2.0 develop will be copied to develop-old and
  develop will become a fork of trunk.

  <\bibliography|bib|tm-plain|fenics-mwa>
    <\bib-list|8>
      <bibitem*|1><label|bib-Ai2012>Haiming<nbsp>Ai, Shuicai<nbsp>Wu,
      Hongjian<nbsp>Gao, Lei<nbsp>Zhao, Chunlan<nbsp>Yang<localize| and
      >Yi<nbsp>Zeng.<newblock> Temperature distribution analysis of tissue
      water vaporization during microwave ablation: experiments and
      simulations.<newblock> <with|font-shape|italic|International journal of
      hyperthermia : the official journal of European Society for
      Hyperthermic Oncology, North American Hyperthermia Group>,
      28(7):674--85, jan 2012.<newblock>

      <bibitem*|2><label|bib-Gas2012>Piotr<nbsp>Gas.<newblock> Temperature
      Distribution of Human Tissue in Interstitial Microwave
      Hyperthermia.<newblock> <with|font-shape|italic|PRZEGLĄD
      ELEKTROTECHNICZNY>, 88(7a):144--146, 2012.<newblock>

      <bibitem*|3><label|bib-Haemmerich2003a>Dieter<nbsp>Haemmerich,
      Louay<nbsp>Chachati, Andrew S<nbsp>Wright, David M<nbsp>Mahvi, Fred
      T<nbsp>Lee<localize| and >John G<nbsp>Webster.<newblock> Hepatic
      radiofrequency ablation with internally cooled probes: effect of
      coolant temperature on lesion size.<newblock>
      <with|font-shape|italic|IEEE transactions on bio-medical engineering>,
      50(4):493--500, apr 2003.<newblock>

      <bibitem*|4><label|bib-Jin2009>Jian-Ming<nbsp>Jin<localize| and
      >Douglas J.<nbsp>Riley.<newblock> <with|font-shape|italic|Finite
      Element Analysis of Antennas and Arrays>.<newblock> John Wiley & Sons,
      2009.<newblock>

      <bibitem*|5><label|bib-Marais2011>Neilen<nbsp>Marais.<newblock>
      Electrical Dipole and First Order Absorbing Boundary Condition (ABC)
      with Dolfin.<newblock> 2011.<newblock>

      <bibitem*|6><label|bib-Nevels1985>R.D.<nbsp>Nevels,
      C.M.<nbsp>Butler<localize| and >W.<nbsp>Yablon.<newblock> The Annular
      Slot Antenna in a Lossy Biological Medium (Short Papers).<newblock>
      <with|font-shape|italic|IEEE Transactions on Microwave Theory and
      Techniques>, 33(4):314--319, apr 1985.<newblock>

      <bibitem*|7><label|bib-Trujillo2013>Macarena<nbsp>Trujillo<localize|
      and >Enrique<nbsp>Berjano.<newblock> Review of the mathematical
      functions used to model the temperature dependence of electrical and
      thermal conductivities of biological tissue in radiofrequency
      ablation.<newblock> <with|font-shape|italic|International journal of
      hyperthermia : the official journal of European Society for
      Hyperthermic Oncology, North American Hyperthermia Group>, , jul
      2013.<newblock>

      <bibitem*|8><label|bib-Vyas1992>Reeta<nbsp>Vyas.<newblock> Green’s
      function solution to the tissue bioheat equation.<newblock>
      <with|font-shape|italic|Medical Physics>, 19(5):1319, sep
      1992.<newblock>
    </bib-list>
  </bibliography>
</body>

<\initial>
  <\collection>
    <associate|language|american>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|4>>
    <associate|auto-10|<tuple|2|6>>
    <associate|auto-11|<tuple|3.2|7>>
    <associate|auto-12|<tuple|3.2.1|7>>
    <associate|auto-13|<tuple|3.2.2|7>>
    <associate|auto-14|<tuple|3.2.2.1|7>>
    <associate|auto-15|<tuple|3.2.3|7>>
    <associate|auto-16|<tuple|3.3|8>>
    <associate|auto-17|<tuple|3.3.1|8>>
    <associate|auto-18|<tuple|3.3.2|9>>
    <associate|auto-19|<tuple|3.3.3|9>>
    <associate|auto-2|<tuple|1.1|4>>
    <associate|auto-20|<tuple|3.3.4|9>>
    <associate|auto-21|<tuple|3.3.5|9>>
    <associate|auto-22|<tuple|3.3.6|10>>
    <associate|auto-23|<tuple|3.3.7|10>>
    <associate|auto-24|<tuple|3.3.8|11>>
    <associate|auto-25|<tuple|3.3.8.1|11>>
    <associate|auto-26|<tuple|3.3.8.2|12>>
    <associate|auto-27|<tuple|3.3.8.3|12>>
    <associate|auto-28|<tuple|3.3.9|12>>
    <associate|auto-29|<tuple|3.3.10|13>>
    <associate|auto-3|<tuple|2|5>>
    <associate|auto-30|<tuple|3.3.11|13>>
    <associate|auto-31|<tuple|3|13>>
    <associate|auto-32|<tuple|3.4|13>>
    <associate|auto-33|<tuple|4|13>>
    <associate|auto-34|<tuple|4.1|14>>
    <associate|auto-35|<tuple|4.2|14>>
    <associate|auto-36|<tuple|4.3|14>>
    <associate|auto-37|<tuple|4.4|14>>
    <associate|auto-38|<tuple|4.5|14>>
    <associate|auto-39|<tuple|4.6|15>>
    <associate|auto-4|<tuple|1|5>>
    <associate|auto-40|<tuple|4.7|15>>
    <associate|auto-41|<tuple|4.8|15>>
    <associate|auto-42|<tuple|4.9|15>>
    <associate|auto-43|<tuple|4.10|15>>
    <associate|auto-44|<tuple|4.11|16>>
    <associate|auto-45|<tuple|4.12|16>>
    <associate|auto-46|<tuple|4.13|16>>
    <associate|auto-47|<tuple|4.14|16>>
    <associate|auto-48|<tuple|4.15|17>>
    <associate|auto-49|<tuple|4.16|17>>
    <associate|auto-5|<tuple|3|5>>
    <associate|auto-50|<tuple|4.17|17>>
    <associate|auto-51|<tuple|5|17>>
    <associate|auto-52|<tuple|5.1|17>>
    <associate|auto-53|<tuple|5.2|17>>
    <associate|auto-54|<tuple|5.3|18>>
    <associate|auto-55|<tuple|5.4|18>>
    <associate|auto-56|<tuple|5.4|18>>
    <associate|auto-57|<tuple|5.4|?>>
    <associate|auto-6|<tuple|3.1|5>>
    <associate|auto-7|<tuple|1|6>>
    <associate|auto-8|<tuple|3.1.1|6>>
    <associate|auto-9|<tuple|3.1.2|6>>
    <associate|bib-Ai2012|<tuple|1|18>>
    <associate|bib-Gas2012|<tuple|2|18>>
    <associate|bib-Haemmerich2003a|<tuple|3|19>>
    <associate|bib-Jin2009|<tuple|4|19>>
    <associate|bib-Marais2011|<tuple|5|19>>
    <associate|bib-Nevels1985|<tuple|6|19>>
    <associate|bib-Trujillo2013|<tuple|7|19>>
    <associate|bib-Vyas1992|<tuple|8|19>>
    <associate|eqn-curlH|<tuple|2|8>>
    <associate|eqn-curlcurlH|<tuple|3|8>>
    <associate|eqn-scalar-form|<tuple|4|8>>
    <associate|eqn-vecvar|<tuple|6|7>>
    <associate|eqn-weak|<tuple|5|7>>
    <associate|fig:code-schem|<tuple|1|5>>
    <associate|overall|<tuple|3.3.11|11>>
    <associate|ref-dispersive-material|<tuple|4|16>>
    <associate|ref-magvar|<tuple|3|16>>
    <associate|ref-nummeth|<tuple|5|16>>
    <associate|ref-stellenbosch|<tuple|2|16>>
    <associate|ref-tempdist|<tuple|1|16>>
    <associate|sec:test_case|<tuple|4|13>>
    <associate|strong-form|<tuple|1|8>>
    <associate|tab:mwa|<tuple|3|13>>
    <associate|tav:pennes|<tuple|1|6>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      Trujillo2013

      Haemmerich2003a

      Trujillo2013,Trujillo2012

      Gas2012

      Gas2012

      Marais2011

      Jin2009

      Jin2009

      Jin2009

      Jin2009

      Gas2012

      Nevels1985

      Gas2012

      Ai2012

      Nevels1985

      Vyas1992
    </associate>
    <\associate|figure>
      <tuple|normal|<label|fig:code-schem>Schematic of the code highlighting
      the importance of the Pennes equation solver. The information being
      passed between the components is given explicitly:
      <with|mode|<quote|math>|T> is the temperature,
      <with|mode|<quote|math>|V> is the cell viability and SAR the specific
      absorption rate.|<pageref|auto-4>>
    </associate>
    <\associate|table>
      <tuple|normal|<label|tav:pennes>|<pageref|auto-7>>

      <tuple|normal||<pageref|auto-10>>

      <tuple|normal|<label|tab:mwa>|<pageref|auto-31>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Getting
      Started> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|1.1<space|2spc>Mesh Generation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Code
      Structure> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Code
      Components> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|3.1<space|2spc>Apparent heat capacity form
      of the Pennes equation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|3fn>|3.1.1<space|2spc>Thermal conductivity
      <with|mode|<quote|math>|k> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|3fn>|3.1.2<space|2spc>Steady-state temperature
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|<quote|1.5fn>|3.2<space|2spc>Electric Potential Solver
      (RFA) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|3fn>|3.2.1<space|2spc>Electrical conductivity
      <with|mode|<quote|math>|\<sigma\>> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <with|par-left|<quote|3fn>|3.2.2<space|2spc>Control System
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>

      <with|par-left|<quote|6fn>|Impedance
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14><vspace|0.15fn>>

      <with|par-left|<quote|3fn>|3.2.3<space|2spc>Computation of useful
      quantities. <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|1.5fn>|3.3<space|2spc>Vector Helmholtz Solver
      (MWA) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|<quote|3fn>|3.3.1<space|2spc>Scalar Form
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17>>

      <with|par-left|<quote|3fn>|3.3.2<space|2spc>Scalar Variational Form
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18>>

      <with|par-left|<quote|3fn>|3.3.3<space|2spc>Expanded integrals
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19>>

      <with|par-left|<quote|3fn>|3.3.4<space|2spc>Expanding complex-valued
      functions <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20>>

      <with|par-left|<quote|3fn>|3.3.5<space|2spc>Vector Variational Form
      With Inhomogeneous Material Properties
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21>>

      <with|par-left|<quote|3fn>|3.3.6<space|2spc>First-Order Absorbing
      Boundary Condition (Scattering) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-22>>

      <with|par-left|<quote|3fn>|3.3.7<space|2spc>Waveguide Port Boundary
      Condition <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-23>>

      <with|par-left|<quote|3fn>|3.3.8<space|2spc>Verification Cases
      (Axisymmetric) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-24>>

      <with|par-left|<quote|6fn>|Coaxial cable analytic solution
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-25><vspace|0.15fn>>

      <with|par-left|<quote|6fn>|Concentric Cylinders
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-26><vspace|0.15fn>>

      <with|par-left|<quote|6fn>|Annular Slot Antenna Test Case
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-27><vspace|0.15fn>>

      <with|par-left|<quote|3fn>|3.3.9<space|2spc>Dependencies of Dielectric
      Tissue Properties <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-28>>

      <with|par-left|<quote|3fn>|3.3.10<space|2spc>Probes
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-29>>

      <with|par-left|<quote|3fn>|3.3.11<space|2spc>Code Parameters
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-30>>

      <with|par-left|<quote|1.5fn>|3.4<space|2spc>Cell Death Model (RFA &
      MWA) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-32>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Test
      Cases> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-33><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|4.1<space|2spc>Solidification due to line
      heat sink <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-34>>

      <with|par-left|<quote|1.5fn>|4.2<space|2spc>Empirical SAR
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-35>>

      <with|par-left|<quote|1.5fn>|4.3<space|2spc>Simple Coaxial Antenna
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-36>>

      <with|par-left|<quote|1.5fn>|4.4<space|2spc>Cell death test
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-37>>

      <with|par-left|<quote|1.5fn>|4.5<space|2spc>Coaxial cable
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-38>>

      <with|par-left|<quote|1.5fn>|4.6<space|2spc>Concentric cylinder test
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-39>>

      <with|par-left|<quote|1.5fn>|4.7<space|2spc>Concentric cylinder
      analytic solution <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-40>>

      <with|par-left|<quote|1.5fn>|4.8<space|2spc>Method of Manufactured
      solutions steady-state Pennes equation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-41>>

      <with|par-left|<quote|1.5fn>|4.9<space|2spc>Method of Manufactured
      solutions transient nonlinear Pennes equation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-42>>

      <with|par-left|<quote|1.5fn>|4.10<space|2spc>Transient Laser Heating
      Benchmark <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-43>>

      <with|par-left|<quote|1.5fn>|4.11<space|2spc>Single Tine RFA
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-44>>

      <with|par-left|<quote|1.5fn>|4.12<space|2spc>Single Tine RFA
      convergence tests <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-45>>

      <with|par-left|<quote|1.5fn>|4.13<space|2spc>Annular Slot Antenna
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-46>>

      <with|par-left|<quote|1.5fn>|4.14<space|2spc>Temperature Dependent
      Electrical Properties (MWA) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-47>>

      <with|par-left|<quote|1.5fn>|4.15<space|2spc>AMICA Prototype 2003
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-48>>

      <with|par-left|<quote|1.5fn>|4.16<space|2spc>AMICA Prototype 2003
      Comparison <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-49>>

      <with|par-left|<quote|1.5fn>|4.17<space|2spc>AMICA Prototype 2011
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-50>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>Implementation
      Notes> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-51><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|5.1<space|2spc>Parameter Dependencies
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-52>>

      <with|par-left|<quote|1.5fn>|5.2<space|2spc>Nedelec Elements in MWA
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-53>>

      <with|par-left|<quote|1.5fn>|5.3<space|2spc>Strang Splitting
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-54>>

      <with|par-left|<quote|1.5fn>|5.4<space|2spc>git Repository Details
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-55>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-56><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>