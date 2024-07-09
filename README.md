# NX-Teracheminterface
This Perl code enables Surface Hopping dynamics on Newton-X with TeraChem.

This code requires Newton-X installed on your computer along with TeraChem. The interface was created by Vaibhav Singh*, Tolga Karsilli**, Todd Martinez***, Spiridoula Matsika*
* Temple University, Philadelphia, USA
** Stanford university, USA
** University of Loisiana, Lafayette, USA

Instructions to use the code:

Note - The use of code requires you to be familiar with Newton-X software (https://newtonx.org/documentation-tutorials/) as well as TeraChem (http://www.petachem.com/products.html).



Described below are conditions required in order to run Terachem with Newton-X

A. Running JOB_AD (adiabatic) with CASSCF

1) The geometry file has to be saved as geom.xyz
2) The output file has to be saved as geom.out (If Terachem job is run with geom.xyz, the output file is automatically geom.out) 
3) The gradient file grad.xyz should be saved in scr.geom (If Terachem job is run with geom.xyz, grad.xyz is automatically saved in scr.geom) 

B. Running JOB_NAD (non-adiabatic) with CASSCF

1) make 3 directories in JOB_NAD as described below
	i)   NAD_N1N0 : This directory should contain the non-adiabatic couplings in between the current state and the state below. Note-If the current state is the ground state, just copy all the files from NAD_N1N2 here.
	ii)  NAD_N1N2 : This directory should contain the non-adiabatic couplings in between the current state and the state above. In the case of the ground state, the coupling is in between state 1 and state 2. Note: if the current state is the highest excited state, just copy all the files from NAD_N1N0 here. 
	iii) GRAD     : This directory just calculates the gradient of the current state.

2) In addition to the 1st point above, all the points in 'A' should also be satisfied. Note that NACs in Terachem are saved as grad.xyz, here under the directories JOB_NAD/NAD_N1Ni/

3) Input files required:
   a) sh.inp is a important input file which should be present for running dynamics. If not provided vdoth keyword is by default=2 which is not reading non-adiabatic couplings from Terachem output. A sample of sh.out is provided below:

 &shinp
        vdoth        = 0
        integrator   = 5
        phase        = 1
        nohop        = 0
        forcesurf    = 1
        getphase     = 1
        nrelax       = 0
        tully        = 1
        ms           = 20
        mom          = 1
        adjmom       = -1
        probmin      = 0
        popdev       = 0.05
        decohmod     = 1
        decay        = 0.1
        run_complex  = 0
        seed         = 1
/

   b) control.dyn sample: Note that the program key (prog) for Terachem given below is 21.0. Other keywords are the same as keywords generated from NX/nxinp for COLUMBUS. I will suggest generating control.dyn by running $NX/nxinp for COLUMBUS and then change the prog keyword to 21.0, if you are not sure about other keywords. 

 &input
        nat        = 6
        nstat      = 3
        nstatdyn   = 2
        dt         = 0.5
        tmax       = 24
        prog       = 21.0
        thres      = 100
        killstat   = 1
        timekill   = 0
        ndamp      = 0
        lvprt      = 1
        kt         = 1
        etot_jump  = 0.5
        etot_drift = 0.5
        nxrestart  = 0
/

   c) geom and veloc of course should be provided containing geometry and the velocity.


  



