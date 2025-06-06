import os
import sys
import numpy as np
import matplotlib as mb
import matplotlib.pyplot as plt
import datetime
import numpy as np
from pathlib import Path

oldGetelecPath = str(Path(__file__).parents[1].absolute()) + '/oldGETELEC/python'
sys.path.insert(0,oldGetelecPath)

newGetelecPath = str(Path(__file__).parents[1].absolute()) + '/src'
sys.path.insert(0,newGetelecPath)

import getelec as gt
import getelec_mod as getelec_old

def _Load_Semiconductor_Data():
    try:
        field = np.load("tests/20221024_semiconductor_test_data/semi_field.npy")
        radius = np.load("tests/20221024_semiconductor_test_data/semi_radius.npy")
        gamma = np.load("tests/20221024_semiconductor_test_data/semi_gamma.npy")
        ec = np.load("tests/20221024_semiconductor_test_data/semi_ec.npy")
        ef = np.load("tests/20221024_semiconductor_test_data/semi_ef.npy")
        eg = np.load("tests/20221024_semiconductor_test_data/semi_eg.npy")
        temp = np.load("tests/20221024_semiconductor_test_data/semi_temp.npy")
        current = np.load("tests/20221024_semiconductor_test_data/semi_current.npy")
        heat = np.load("tests/20221024_semiconductor_test_data/semi_heat.npy")
        return field, radius, gamma, ec, ef, eg, temp, current, heat

    except(IOError):
        print("Semiconductor test data not found")
        return False
   
def _Load_Metal_Data():
    try:
        field = np.load("tests/20221024_metal_test_data/metal_field.npy")
        radius = np.load("tests/20221024_metal_test_data/metal_radius.npy")
        gamma = np.load("tests/20221024_metal_test_data/metal_gamma.npy")
        ef = np.load("tests/20221024_metal_test_data/metal_ef.npy")
        temp = np.load("tests/20221024_metal_test_data/metal_temp.npy")
        current = np.load("tests/20221024_metal_test_data/metal_current.npy")
        heat = np.load("tests/20221024_metal_test_data/metal_heat.npy")
        return field, radius, gamma, ef, temp, current, heat

    except(IOError):
        print("Metal test data not found")
        return False
   
def _Save_Table_to_Files():
    path_to_file = "/home/salva/Documents/getelec_priv/python/spectra"
    data = np.loadtxt(path_to_file)

    field = data[:,0]
    radius = data[:,1]
    gamma = data[:,2]
    ec = data[:,3]
    ef = data[:,4]
    eg = data[:,5]
    temp = data[:,6]

    try:
        os.mkdir("tests/20221024_metal_test_data")
    except:
        pass

    np.save("tests/20221024_metal_test_data/metal_field", field)
    np.save("tests/20221024_metal_test_data/metal_radius", radius)
    np.save("tests/20221024_metal_test_data/metal_gamma", gamma)
    #np.save("tests/20221024_semiconductor_test_data/semi_ec", ec)
    np.save("tests/20221024_metal_test_data/metal_ef", ef)
    #np.save("tests/20221024_semiconductor_test_data/semi_eg", eg)
    np.save("tests/20221024_metal_test_data/metal_temp", temp)
    
    current=gt.current_metal_emitter(field,radius,gamma,ef,temp)
    heat=gt.heat_metal_emitter(field,radius,gamma,ef,temp)

    np.save("tests/20221024_metal_test_data/metal_current",current)
    np.save("tests/20221024_metal_test_data/metal_heat",heat)

    return True

def _Getelec_Installation_Semiconductor_Test():
    
    ref_field, ref_radius, ref_gamma, ref_ec, ref_ef, ref_eg, ref_temp, J_ref_getelec, Pn_ref_getelec = _Load_Semiconductor_Data()

    J_new_getelec = gt.current_semiconductor_emitter(ref_field, ref_radius, ref_gamma, ref_ec, ref_ef, ref_eg, ref_temp)
    Pn_new_getelec = gt.heat_semiconductor_emitter(ref_field, ref_radius, ref_gamma, ref_ec, ref_ef, ref_eg, ref_temp)

    J_absolute_error = abs(J_new_getelec - J_ref_getelec)

    if J_absolute_error.any() != 0:
        J_relative_error = J_absolute_error / J_ref_getelec
        J_bad = np.where(np.logical_and(J_relative_error > 0.5, J_absolute_error > 1.e-25))[0]
        J_rms_error = np.sqrt(np.mean(J_relative_error[J_absolute_error > 1.e-25]**2))
    else:
        J_relative_error = np.zeros(len(J_absolute_error))
        J_bad = np.where(np.logical_and(J_relative_error > 0.5, J_absolute_error > 1.e-25))[0]
        J_rms_error = np.zeros(len(J_absolute_error))
    
    Pn_absolute_error = abs(Pn_new_getelec - Pn_ref_getelec)
    if Pn_absolute_error.any() !=0:
        Pn_relative_error = Pn_absolute_error / Pn_ref_getelec
        Pn_bad = np.where(np.logical_and(Pn_relative_error > 0.5, Pn_absolute_error > 1.e-25))[0]
        Pn_rms_error = np.sqrt(np.mean(Pn_relative_error[Pn_absolute_error > 1.e-25]**2))
    else:
        Pn_relative_error = np.zeros(len(Pn_absolute_error))
        Pn_bad = np.where(np.logical_and(Pn_relative_error > 0.5, Pn_absolute_error > 1.e-25))[0]
        Pn_rms_error = np.zeros(len(Pn_absolute_error))
    
    

    if J_rms_error.any() > 0.1:
        print("\nSemiconductor current test: NOT PASSED\nDebug your installation or seek assistance.")
        print("\nrms error in J = ", J_rms_error)
        print("J bad = ", J_bad)
        for i in J_bad:
            print("J_ref_getelec_bad", J_ref_getelec[i]," J_new_getelec_bad", J_new_getelec[i])
            print("     Parametres:",ref_field[i], ref_radius[i], ref_gamma[i], ref_ec[i], ref_ef[i], ref_eg[i], ref_temp[i])
            #em.set(Field[i], Radius[i], Gamma[i])
            #em.interpolate()
            #em.get_lims(Workfunction[i], kT[i])
            #em.integrate_quad(Workfunction[i], kT[i])
            #em.integrate_quad_Nottingham(Workfunction[i], kT[i])
    else:
        print("\nSemiconductor current test: PASSED")

    if Pn_rms_error.any() > .1:
        print("\n\nSemiconductor Nottighamm heat test: NOT PASSED\nDebug your installation or seek assistance.")
        print("\nrms error in Pn = ", Pn_rms_error)
        print("Pn bad = ", Pn_bad)
        for i in Pn_bad:
            print("Pn_ref_getelec_bad", Pn_ref_getelec[i]," Pn_new_getelec_bad", Pn_new_getelec[i])
            print("     Parametres:",ref_field[i], ref_radius[i], ref_gamma[i], ref_ec[i], ref_ef[i], ref_eg[i], ref_temp[i])
            #em.set(Field[i], Radius[i], Gamma[i])
            #em.interpolate()
            #em.get_lims(Workfunction[i], kT[i])
            #em.integrate_quad(Workfunction[i], kT[i])
            #em.integrate_quad_Nottingham(Workfunction[i], kT[i])
    else:
        print("\nSemiconductor Nottigham heat test: PASSED")        

    if len(J_bad)==0 and len(Pn_bad)==0 and Pn_rms_error.all()==0 and J_rms_error.all()==0:
        print("\nSemiconductor installation successful") 
    else:
        print("Semiconductor intstallation unsuccessful. Debug yout installation or seek assistance")

        font = 20
        x = 15
        y = x

        mb.rcParams["font.size"] = font
        mb.rcParams["axes.labelsize"] = font
        mb.rcParams["xtick.labelsize"] = font
        mb.rcParams["ytick.labelsize"] = font
        #mb.rcParams["legend.fontsize"] = font
        mb.rcParams["lines.linewidth"] = 2

        fig1, ax = plt.subplots(figsize=(x,y)) 
        #ax.ticklabel_format(scilimits=[-1,1])
        plt.loglog(J_ref_getelec, J_new_getelec, '.' ,color = "steelblue")
        plt.loglog(abs(J_ref_getelec), abs(J_new_getelec), '.',color = "steelblue")
        plt.loglog([1.e-80, 1.], [1.e-80, 1.], color = "orange", label = "1to1 line")
        plt.grid("True")
        plt.legend()
        plt.xlim([10E-25,10E0])
        plt.ylim([10E-25,10E0])
        plt.xlabel("Ref GETELEC current densities (A$nm^{-2}$)")
        plt.ylabel("New GETELEC current densities (A$nm^{-2}$)")
        plt.title("Comparing current from new and ref GETELEC")
        plt.savefig("Installation_test: J.png")
        #plt.show()

        fig2, ax = plt.subplots(figsize=(x,y))
        plt.loglog(Pn_ref_getelec, Pn_new_getelec, '.',color = "steelblue")
        plt.loglog(abs(Pn_ref_getelec), abs(Pn_new_getelec), '.',color = "steelblue")
        plt.loglog([1.e-80, 1.], [1.e-80, 1.],color = "orange", label = "1to1 line")
        plt.grid("True")
        plt.legend()
        plt.xlim([10E-25,10E0])
        plt.ylim([10E-25,10E0])
        plt.xlabel("Ref GETELEC Nottigham heat (W$nm^{-2}$)")
        plt.ylabel("New GETELEC Nottigham heat (W$nm^{-2}$)")
        plt.title("Comparing heat from new and ref GETELEC")
        plt.savefig("Installation_test_Pn.png")
        #plt.show()

    return True

def _Getelec_Installation_Metal_Test():
    
    ref_field, ref_radius, ref_gamma, ref_ef, ref_temp, J_ref_getelec, Pn_ref_getelec = _Load_Metal_Data()

    J_new_getelec = gt.current_metal_emitter(ref_field, ref_radius, ref_gamma, ref_ef, ref_temp)
    Pn_new_getelec = gt.heat_metal_emitter(ref_field, ref_radius, ref_gamma, ref_ef, ref_temp)

    J_absolute_error = abs(J_new_getelec - J_ref_getelec)

    if J_absolute_error.any() != 0:
        J_relative_error = J_absolute_error / J_ref_getelec
        J_bad = np.where(np.logical_and(J_relative_error > 0.5, J_absolute_error > 1.e-25))[0]
        J_rms_error = np.sqrt(np.mean(J_relative_error[J_absolute_error > 1.e-25]**2))
    else:
        J_relative_error = np.zeros(len(J_absolute_error))
        J_bad = np.where(np.logical_and(J_relative_error > 0.5, J_absolute_error > 1.e-25))[0]
        J_rms_error = np.zeros(len(J_absolute_error))
    
    Pn_absolute_error = abs(Pn_new_getelec - Pn_ref_getelec)
    if Pn_absolute_error.any() !=0:
        Pn_relative_error = Pn_absolute_error / Pn_ref_getelec
        Pn_bad = np.where(np.logical_and(Pn_relative_error > 0.5, Pn_absolute_error > 1.e-25))[0]
        Pn_rms_error = np.sqrt(np.mean(Pn_relative_error[Pn_absolute_error > 1.e-25]**2))
    else:
        Pn_relative_error = np.zeros(len(Pn_absolute_error))
        Pn_bad = np.where(np.logical_and(Pn_relative_error > 0.5, Pn_absolute_error > 1.e-25))[0]
        Pn_rms_error = np.zeros(len(Pn_absolute_error))
    
    

    if J_rms_error.any() > .1:
        print("\nMetal current test: NOT PASSED\nDebug your installation or seek assistance.")
        print("\nrms error in J = ", J_rms_error)
        print("J bad = ", J_bad)
        for i in J_bad:
            print("J_ref_getelec_bad", J_ref_getelec[i]," J_new_getelec_bad", J_new_getelec[i])
            print("     Parametres:",ref_field[i], ref_radius[i], ref_gamma[i], ref_ef[i], ref_temp[i])
            #em.set(Field[i], Radius[i], Gamma[i])
            #em.interpolate()
            #em.get_lims(Workfunction[i], kT[i])
            #em.integrate_quad(Workfunction[i], kT[i])
            #em.integrate_quad_Nottingham(Workfunction[i], kT[i])
    else:
        print("\nMetal current test: PASSED")

    if Pn_rms_error.any() > .1:
        print("\n\nMetal Nottighamm heat test: NOT PASSED\nDebug your installation or seek assistance.")
        print("\nrms error in Pn = ", Pn_rms_error)
        print("Pn bad = ", Pn_bad)
        for i in Pn_bad:
            print("Pn_ref_getelec_bad", Pn_ref_getelec[i]," Pn_new_getelec_bad", Pn_new_getelec[i])
            print("     Parametres:",ref_field[i], ref_radius[i], ref_gamma[i], ref_ef[i], ref_temp[i])
            #em.set(Field[i], Radius[i], Gamma[i])
            #em.interpolate()
            #em.get_lims(Workfunction[i], kT[i])
            #em.integrate_quad(Workfunction[i], kT[i])
            #em.integrate_quad_Nottingham(Workfunction[i], kT[i])
    else:
        print("\nMetal Nottigham heat test: PASSED")        

    if len(J_bad)==0 and len(Pn_bad)==0 and Pn_rms_error.all()==0 and J_rms_error.all()==0:
        print("\nMetal installation successful") 
    else:
        print("Metal intstallation unsuccessful. Debug yout installation or seek assistance")

        font = 20
        x = 15
        y = x

        mb.rcParams["font.size"] = font
        mb.rcParams["axes.labelsize"] = font
        mb.rcParams["xtick.labelsize"] = font
        mb.rcParams["ytick.labelsize"] = font
        #mb.rcParams["legend.fontsize"] = font
        mb.rcParams["lines.linewidth"] = 2

        fig1, ax = plt.subplots(figsize=(x,y)) 
        #ax.ticklabel_format(scilimits=[-1,1])
        plt.loglog(J_ref_getelec, J_new_getelec, '.' ,color = "steelblue")
        plt.loglog(abs(J_ref_getelec), abs(J_new_getelec), '.',color = "steelblue")
        plt.loglog([1.e-80, 1.], [1.e-80, 1.], color = "orange", label = "1to1 line")
        plt.grid("True")
        plt.legend()
        plt.xlim([10E-25,10E0])
        plt.ylim([10E-25,10E0])
        plt.xlabel("Ref GETELEC current densities (A$nm^{-2}$)")
        plt.ylabel("New GETELEC current densities (A$nm^{-2}$)")
        plt.title("Comparing current from new and ref GETELEC")
        plt.savefig("Installation_test_J.png")
        #plt.show()

        fig2, ax = plt.subplots(figsize=(x,y))
        plt.loglog(Pn_ref_getelec, Pn_new_getelec, '.',color = "steelblue")
        plt.loglog(abs(Pn_ref_getelec), abs(Pn_new_getelec), '.',color = "steelblue")
        plt.loglog([1.e-80, 1.], [1.e-80, 1.],color = "orange", label = "1to1 line")
        plt.grid("True")
        plt.legend()
        plt.xlim([10E-25,10E0])
        plt.ylim([10E-25,10E0])
        plt.xlabel("Ref GETELEC Nottigham heat (W$nm^{-2}$)")
        plt.ylabel("New GETELEC Nottigham heat (W$nm^{-2}$)")
        plt.title("Comparing heat from new and ref GETELEC")
        plt.savefig("Installation_test_Pn.png")
        #plt.show()

    return True

def Randomised_Tabulator_Test():

    tab = gt.Interpolator()

    Fmax = 1/tab.Finv[0]
    Fmin = 1/tab.Finv[-1]
    Rmax = 1/tab.Rinv[0]
    Rmin = 1/tab.Rinv[-1]
    gammax = 1/tab.gaminv[0]
    gammin = 1/tab.gaminv[-1]

    Np = 8096

    Field = np.random.rand(Np) * (Fmax - Fmin) + Fmin
    Radius = np.random.rand(Np) * (Rmax - Rmin) + Rmin
    Gamma = np.random.rand(Np) * (gammax - gammin) + gammin
    Workfunction = np.random.rand(Np) * (7.5 - 2.5) + 2.5
    Temperature = np.random.rand(Np) * (3000 - 100) + 100
    kT = Temperature * 8.6173324e-5 

    J_new_getelec = np.copy(Field)
    Pn_new_getelec = np.copy(Field)
    J_ref_getelec = np.copy(Field)
    Pn_ref_getelec = np.copy(Field)

    metal_emitter = gt.Metal_Emitter(tab)

    print("\nCalculating from NEW GETELEC")
    new_getelec_start = datetime.datetime.now()
    for i in range(len(Field)):
        metal_emitter._emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        metal_emitter._emitter.Interpolate_Gammow()
        metal_emitter.Define_Metal_Emitter_Parameters(Workfunction[i], kT[i])
        J_new_getelec[i] = metal_emitter.currentDensity()
        Pn_new_getelec[i] = metal_emitter.Nottingham_Heat_from_Metals()
    new_getelec_end = datetime.datetime.now()

    print("\nCalculating from REF GETELEC\n")
    ref_getelec_start = datetime.datetime.now()
    em = getelec_old.emission_create(approx=2)
    for i in range(len(Field)):   
        em.F = Field[i]
        em.W = Workfunction[i]
        em.Temp = Temperature[i]
        em.gamma = Gamma[i]
        em.R = Radius[i]
        em.cur_dens()
        J_ref_getelec[i] = em.Jem
        Pn_ref_getelec[i] = em.heat
    ref_getelec_end = datetime.datetime.now()


    J_absolute_error = abs(J_new_getelec - J_ref_getelec)
    J_relative_error = J_absolute_error / J_ref_getelec
    J_bad = np.where(np.logical_and(J_relative_error > 0.5, J_absolute_error > 1.e-25))[0]

    Pn_absolute_error = abs(Pn_new_getelec - Pn_ref_getelec)
    Pn_relative_error = Pn_absolute_error / Pn_ref_getelec
    Pn_bad = np.where(np.logical_and(Pn_relative_error > 0.5, Pn_absolute_error > 1.e-25))[0]

    J_rms_error = np.sqrt(np.mean(J_relative_error[J_absolute_error > 1.e-25]**2))
    Pn_rms_error = np.sqrt(np.mean(Pn_relative_error[Pn_absolute_error > 1.e-25]**2))

    if J_rms_error > 0.01:
        print("\nTabulator current test: NOT PASSED")
        if len(J_bad)!=0:
            print("WARNING!\nSome values do not complie with the rms error tolerance.\nDebug your code or seek assistance.\n")
            print("\nrms error in J = ", J_rms_error)
            print("J bad = ", J_bad)
            for i in J_bad:
               print("J_old_getelec_bad", J_ref_getelec[i]," J_new_getelec_bad", J_new_getelec[i])
               print("     Parametres:",Field[i], Radius[i], Gamma[i], Workfunction[i], Temperature[i])
               #em.set(Field[i], Radius[i], Gamma[i])
               #em.interpolate()
               #em.get_lims(Workfunction[i], kT[i])
               #em.integrate_quad(Workfunction[i], kT[i])
               #em.integrate_quad_Nottingham(Workfunction[i], kT[i])
        else:
            pass
    else:
        print("\nTabulator current test: PASS")

    if Pn_rms_error > 0.01:
        print("\nTabulator Nottigham heat test: NOT PASSED")
        if len(Pn_bad) !=0:
            print("WARNING!\nSome values do not complie with the rms error tolerance.\nDebug your code or seek assistance.\n")
            print("\nrms error in Pn = ", Pn_rms_error)
            print("Pn bad = ", Pn_bad)
            for i in Pn_bad:
               print("Pn_old_getelec_bad", Pn_ref_getelec[i]," Pn_new_getelec_bad", Pn_new_getelec[i])
               print("     Parametres:",Field[i], Radius[i], Gamma[i], Workfunction[i], Temperature[i])
               #em.set(Field[i], Radius[i], Gamma[i])
               #em.interpolate()
               #em.get_lims(Workfunction[i], kT[i])
               #em.integrate_quad(Workfunction[i], kT[i])
               #em.integrate_quad_Nottingham(Workfunction[i], kT[i])
        else:
            pass
    else:
        print("\nTabulator Nottigham heat test: PASSED")

    if len(J_bad)==0 and len(Pn_bad)==0 and Pn_rms_error<=0.01 and J_rms_error<=0.01:
        print("\nTabulator WITHIN the tolerance")
    else:
        print("\nTabulator OUT of tolerance") 

    print("\nRef GETELEC running time =", ref_getelec_end-ref_getelec_start)
    print("New GETELEC running time =", new_getelec_end-new_getelec_start)

    font = 20
    x = 15
    y = x

    mb.rcParams["font.size"] = font
    mb.rcParams["axes.labelsize"] = font
    mb.rcParams["xtick.labelsize"] = font
    mb.rcParams["ytick.labelsize"] = font
    #mb.rcParams["legend.fontsize"] = font
    mb.rcParams["lines.linewidth"] = 2

    fig1, ax = plt.subplots(figsize=(x,y)) 
    #ax.ticklabel_format(scilimits=[-1,1])
    plt.loglog(J_ref_getelec, J_new_getelec, '.' ,color = "steelblue")
    plt.loglog(abs(J_ref_getelec), abs(J_new_getelec), '.',color = "steelblue")
    plt.loglog([1.e-80, 1.], [1.e-80, 1.], color = "orange", label = "1to1 line")
    plt.grid("True")
    plt.legend()
    plt.xlim([10E-90,10E0])
    plt.ylim([10E-90,10E0])
    plt.xlabel("Old GETELEC current densities (A$nm^{-2}$)")
    plt.ylabel("New GETELEC current densities (A$nm^{-2}$)")
    plt.title("Comparing current from new and old GETELEC")
    plt.savefig("Validation_test_J.png")
    #plt.show()

    fig2, ax = plt.subplots(figsize=(x,y))
    plt.loglog(Pn_ref_getelec, Pn_new_getelec, '.',color = "steelblue")
    plt.loglog(abs(Pn_ref_getelec), abs(Pn_new_getelec), '.',color = "steelblue")
    plt.loglog([1.e-80, 1.], [1.e-80, 1.],color = "orange", label = "1to1 line")
    plt.grid("True")
    plt.legend()
    plt.xlim([10E-90,10E0])
    plt.ylim([10E-90,10E0])
    plt.xlabel("Old GETELEC Nottigham heat (W$nm^{-2}$)")
    plt.ylabel("New GETELEC Nottigham heat (W$nm^{-2}$)")
    plt.title("Comparing heat from new and old GETELEC")
    plt.savefig("Validation_test_Pn.png")
    #plt.show()
    return True

def Referenced_Tabulator_Test():
    ref_field, ref_radius, ref_gamma, ref_ef, ref_temp, J_ref_getelec, Pn_ref_getelec = _Load_Metal_Data()
    tab = gt.Interpolator()

    Field = ref_field
    Radius = ref_radius
    Gamma = ref_gamma
    Workfunction = ref_ef
    Temperature = ref_temp
    kT = Temperature * 8.6173324e-5 

    J_new_getelec = np.copy(Field)
    Pn_new_getelec = np.copy(Field)
    J_ref_getelec = np.copy(Field)
    Pn_ref_getelec = np.copy(Field)

    metal_emitter = gt.Metal_Emitter(barrier=gt.Barrier(tabulationFolder='tabulated'), supply=gt.Supply())

    print("\nCalculating from NEW GETELEC")
    new_getelec_start = datetime.datetime.now()
    for i in range(len(Field)):
        metal_emitter._emitter.Define_Barrier_Parameters(Field[i], Radius[i], Gamma[i])
        metal_emitter._emitter.Interpolate_Gammow()
        metal_emitter.Define_Metal_Emitter_Parameters(Workfunction[i], kT[i])
        J_new_getelec[i] = metal_emitter.currentDensity()
        Pn_new_getelec[i] = metal_emitter.Nottingham_Heat_from_Metals()
    new_getelec_end = datetime.datetime.now()

    print("\nCalculating from REF GETELEC\n")
    ref_getelec_start = datetime.datetime.now()
    em = getelec_old.emission_create(approx=2)
    for i in range(len(Field)):   
        em.F = Field[i]
        em.W = Workfunction[i]
        em.Temp = Temperature[i]
        em.gamma = Gamma[i]
        em.R = Radius[i]
        em.cur_dens()
        J_ref_getelec[i] = em.Jem
        Pn_ref_getelec[i] = em.heat
    ref_getelec_end = datetime.datetime.now()


    J_absolute_error = abs(J_new_getelec - J_ref_getelec)
    J_relative_error = J_absolute_error / J_ref_getelec
    J_bad = np.where(np.logical_and(J_relative_error > 0.5, J_absolute_error > 1.e-25))[0]

    Pn_absolute_error = abs(Pn_new_getelec - Pn_ref_getelec)
    Pn_relative_error = Pn_absolute_error / Pn_ref_getelec
    Pn_bad = np.where(np.logical_and(Pn_relative_error > 0.5, Pn_absolute_error > 1.e-25))[0]

    J_rms_error = np.sqrt(np.mean(J_relative_error[J_absolute_error > 1.e-25]**2))
    Pn_rms_error = np.sqrt(np.mean(Pn_relative_error[Pn_absolute_error > 1.e-25]**2))

    if J_rms_error > 0.01:
        print("\nTabulator current test: NOT PASSED")
        print("WARNING!\nSome values do not complie with the rms error tolerance.\nDebug your code or seek assistance.\n")
        print("\nrms error in J = ", J_rms_error)
        print("J bad = ", J_bad)
        for i in J_bad:
            print("J_old_getelec_bad", J_ref_getelec[i]," J_new_getelec_bad", J_new_getelec[i])
            print("     Parametres:",Field[i], Radius[i], Gamma[i], Workfunction[i], Temperature[i])
            #em.set(Field[i], Radius[i], Gamma[i])
            #em.interpolate()
            #em.get_lims(Workfunction[i], kT[i])
            #em.integrate_quad(Workfunction[i], kT[i])
            #em.integrate_quad_Nottingham(Workfunction[i], kT[i])
        else:
            pass
    else:
        print("\nTabulator current test: PASSED")

    if Pn_rms_error > 0.01:
        print("\nTabulator Nottigham heat test: NOT PASSED")
        print("WARNING!\nSome values do not complie with the rms error tolerance.\nDebug your code or seek assistance.\n")
        print("\nrms error in Pn = ", Pn_rms_error)
        print("Pn bad = ", Pn_bad)
        for i in Pn_bad:
            print("Pn_old_getelec_bad", Pn_ref_getelec[i]," Pn_new_getelec_bad", Pn_new_getelec[i])
            print("     Parametres:",Field[i], Radius[i], Gamma[i], Workfunction[i], Temperature[i])
            #em.set(Field[i], Radius[i], Gamma[i])
            #em.interpolate()
            #em.get_lims(Workfunction[i], kT[i])
            #em.integrate_quad(Workfunction[i], kT[i])
            #em.integrate_quad_Nottingham(Workfunction[i], kT[i])
        else:
            pass
    else:
        print("\nTabulator Nottigham heat test: PASSED")

    if len(J_bad)==0 and len(Pn_bad)==0 and Pn_rms_error<=0.01 and J_rms_error<=0.01:
        print("\nTabulator within the tolerance")
    else:
        print("\nTabulator out of tolerance") 

    print("\nRef GETELEC running time =", ref_getelec_end-ref_getelec_start)
    print("New GETELEC running time =", new_getelec_end-new_getelec_start)

    font = 20
    x = 15
    y = x

    mb.rcParams["font.size"] = font
    mb.rcParams["axes.labelsize"] = font
    mb.rcParams["xtick.labelsize"] = font
    mb.rcParams["ytick.labelsize"] = font
    #mb.rcParams["legend.fontsize"] = font
    mb.rcParams["lines.linewidth"] = 2

    fig1, ax = plt.subplots(figsize=(x,y)) 
    #ax.ticklabel_format(scilimits=[-1,1])
    plt.loglog(J_ref_getelec, J_new_getelec, '.' ,color = "steelblue")
    plt.loglog(abs(J_ref_getelec), abs(J_new_getelec), '.',color = "steelblue")
    plt.loglog([1.e-80, 1.], [1.e-80, 1.], color = "orange", label = "1to1 line")
    plt.grid("True")
    plt.legend()
    plt.xlim([10E-90,10E0])
    plt.ylim([10E-90,10E0])
    plt.xlabel("Old GETELEC current densities (A$nm^{-2}$)")
    plt.ylabel("New GETELEC current densities (A$nm^{-2}$)")
    plt.title("Comparing current from new and old GETELEC")
    plt.savefig("Validation_test_J.png")
    #plt.show()

    fig2, ax = plt.subplots(figsize=(x,y))
    plt.loglog(Pn_ref_getelec, Pn_new_getelec, '.',color = "steelblue")
    plt.loglog(abs(Pn_ref_getelec), abs(Pn_new_getelec), '.',color = "steelblue")
    plt.loglog([1.e-80, 1.], [1.e-80, 1.],color = "orange", label = "1to1 line")
    plt.grid("True")
    plt.legend()
    plt.xlim([10E-90,10E0])
    plt.ylim([10E-90,10E0])
    plt.xlabel("Old GETELEC Nottigham heat (W$nm^{-2}$)")
    plt.ylabel("New GETELEC Nottigham heat (W$nm^{-2}$)")
    plt.title("Comparing heat from new and old GETELEC")
    plt.savefig("Validation_test_Pn.png")
    #plt.show()
    return True

def Getelec_Installation_Test():
    _Getelec_Installation_Metal_Test()
    # _Getelec_Installation_Semiconductor_Test()
    return True

Getelec_Installation_Test()

#Randomised_Tabulator_Test()

# region One data point calculation routine
    # This routine calculates the current density from semiconductors (two methods) and metals, as well as the plotting of the energy distributions
"""
Npoly = 5
NGi = 512
zs = 1.6183e-4
kBoltz = 8.6173324e-5 

tab = gtab.Tabulator()

Workfunction = 4.7
Ef = Workfunction

Eg = 0.661
Ec = 4.5
m = 9.1093837015e-31 
me = 1.59*m # effective electron mass @300K, https://doi.org/10.1142/S0219749909005833
mp = 0.33*m # effective hole mass @300k, https://doi.org/10.1142/S0219749909005833
Temp = 300.
kT = kBoltz * Temp

metal_emitter = gtab.Metal_Emitter(tab)

metal_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
metal_emitter.emitter.Interpolate_Gammow()

metal_emitter.Define_Emitter_Parameters(Workfunction, kT)

j_metal = metal_emitter.Current_Density()
pn_metal = metal_emitter.Nottingham_Heat()
energy_space_metal, distribution_metal = metal_emitter.Energy_Distribution()

semiconductor_emitter = gtab.Semiconductor_Emitter(tab)

semiconductor_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
semiconductor_emitter.emitter.Interpolate_Gammow()

semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec, Ef, Eg, kT, m, me, mp)

j_c, j_v, j_total = semiconductor_emitter.Current_Density_from_Semiconductors()

pn_c, pn_v, pn_total = semiconductor_emitter.Nottingham_Heat_from_Semiconductors()

energy_c, distribution_c, energy_v, distribution_v = semiconductor_emitter.Energy_Distribution_from_Semiconductors()

j_c2, j_v2, j_total2 = semiconductor_emitter.Distributed_Current_Density_from_Semiconductors()

c_abs_error = abs(j_c2-j_c)
c_rel_error = c_abs_error/j_c

v_abs_error = abs(j_v2-j_v)
v_rel_error = v_abs_error/j_v

total_abs_error = abs(j_total2-j_total)
total_rel_error = total_abs_error/j_total

print("current_c", j_c, j_c2, "relative error:", c_rel_error)
print("current_v", j_v, j_v2, "relative error:", v_rel_error)
print("current_semi", j_total, j_total2, "relative error:", total_rel_error)
print("current_metal", j_metal)


font = 120
x = 52
y = 40


mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font+5
mb.rcParams["legend.fontsize"] = font/1.5
mb.rcParams["lines.linewidth"] = 7

fig, ax = plt.subplots(figsize=(x,y))
ax.ticklabel_format(scilimits=[-1,1])
#plt.plot(energy_space_metal1-4.5, distribution_metal1, color = "steelblue", label = "electrons from metal WF=4.5") 
#plt.plot(energy_space_metal2-4.5, distribution_metal2, color = "orange", label = "electrons from metal WF=2") 
plt.plot(energy_space_metal-Workfunction, distribution_metal/max(distribution_metal), color = "steelblue", label = "electrons from metal") #
plt.plot(energy_c-Ef, distribution_c/max(distribution_c), color = "orange", label = "electrons from conduction band") 
plt.plot(energy_v-Ef, distribution_v/max(distribution_v), color = "green", label = "electrons from valence band") 
plt.legend()
plt.grid("True")
plt.xlabel("Energy (eV)")
plt.ylabel("$J_{eN}$(E) (A$nm^{-2}$$eV^{-1}$)")
plt.title("Normalised electron energy distribution")
plt.savefig("Normalised electron energy distribution.svg")
plt.savefig("Normalised electron energy distribution.png")

fig, ax = plt.subplots(figsize=(x,y))
ax.ticklabel_format(scilimits=[-1,1])
plt.plot(energy_space_metal-Workfunction, distribution_metal, color = "steelblue", label = "electrons from metal") 
plt.legend()
plt.grid("True")
plt.xlabel("Energy (eV)")
plt.ylabel("$J_e$(E) (A$nm^{-2}$$eV^{-1}$)")
plt.title("Electron energy distribution from a metal")
plt.savefig("Electron energy distribution from a metal.svg")
plt.savefig("Electron energy distribution from a metal.png")

fig, ax = plt.subplots(figsize=(x,y))
ax.ticklabel_format(scilimits=[-1,1])
plt.plot(energy_c-Ef+Eg/2, distribution_c, color = "orange", label = "electrons from conduction band") 
plt.legend()
plt.grid("True")
plt.xlabel("Energy (eV)")
plt.ylabel("$J_e$(E) (A$nm^{-2}$$eV^{-1}$)")
plt.title("Electron energy distribution from conduction band")
plt.savefig("Electron energy distribution from conduction band.svg")
plt.savefig("Electron energy distribution from conduction band.png")

fig, ax = plt.subplots(figsize=(x,y))
ax.ticklabel_format(scilimits=[-1,1])
plt.plot(energy_v-Ef+Eg/2, distribution_v, color = "green", label = "electrons from valence band") 
plt.legend()
plt.grid("True")
plt.xlabel("Energy (eV)")
plt.ylabel("$J_e$(E) (A$nm^{-2}$$eV^{-1}$)")
plt.title("Electron energy distribution from valence band")
plt.savefig("Electron energy distribution from valence band.svg")
plt.savefig("Electron energy distribution from valence band.png")
"""
# endregion

# region Multiple data point calculation routine - metal
"""
#Npoly = 5
#NGi = 128
#zs = 1.6183e-4
kBoltz = 8.6173324e-5 

tab = gtab.Tabulator()

Fmax = 1/tab.Finv[0]
Fmin = 1/tab.Finv[-1]
Rmax = 1/tab.Rinv[0]
Rmin = 1/tab.Rinv[-1]
gammax = 1/tab.gaminv[0]
gammin = 1/tab.gaminv[-1]

Np = 8096

Fi = np.random.rand(Np) * (Fmax - Fmin) + Fmin
Ri = np.random.rand(Np) * (Rmax - Rmin) + Rmin
gami = np.random.rand(Np) * (gammax - gammin) + gammin
Wi = np.random.rand(Np) * (7.5 - 2.5) + 2.5
Ti = np.random.rand(Np) * (3000 - 100) + 100
kT = Ti * kBoltz
Ji = np.copy(Fi)
Pi = np.copy(Fi)
Jget = np.copy(Ji)
Pget = np.copy(Ji)

metal_emitter = gtab.Metal_Emitter(tab)

print("calculating from tabulator")
tab_start = datetime.datetime.now()
for i in range(len(Fi)):
    metal_emitter.emitter.Define_Barrier_Parameters(Fi[i], Ri[i], gami[i])
    metal_emitter.emitter.Interpolate_Gammow()
    metal_emitter.Define_Emitter_Parameters(Wi[i], kT[i])
    Ji[i] = metal_emitter.Current_Density()
    Pi[i] = metal_emitter.Nottingham_Heat()
tab_end = datetime.datetime.now()

print("calculating from getelec")
get_start = datetime.datetime.now()
em = gt.emission_create(approx=2)
for i in range(len(Fi)):   
    em.F = Fi[i]
    em.W = Wi[i]
    em.Temp = Ti[i]
    em.gamma = gami[i]
    em.R = Ri[i]
    em.cur_dens()
    Jget[i] = em.Jem
    Pget[i] = em.heat
get_end = datetime.datetime.now()


abserr = abs(Ji - Jget)
relerr = abserr / Jget
bad = np.where(np.logical_and(relerr > 0.5, abserr > 1.e-25))[0]

Pn_abserr = abs(Pi - Pget)
Pn_relerr = Pn_abserr / Pget
Pbad = np.where(np.logical_and(Pn_relerr > 0.5, Pn_abserr > 1.e-25))[0]

print("bad = ", bad)
print("rms error in J = ", np.sqrt(np.mean(relerr[abserr > 1.e-25]**2)))
print("rms error in Pn = ", np.sqrt(np.mean(Pn_relerr[Pn_abserr > 1.e-25]**2)))
print("getelec running time =", get_end-get_start)
print("tabulat running time =", tab_end-tab_start)


#for i in bad:
#    print("Jget, Ji : ", Jget[i], Ji[i])
#    emit.set(Fi[i], Ri[i], gami[i])
#    emit.interpolate()
#    emit.get_lims(Wi[i], kT[i])
#    emit.integrate_quad(Wi[i], kT[i])
#    emit.integrate_quad_Nottingham(W[i], kT[i])

font = 20
x = 15
y = x

mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font
#mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 2

fig1, ax = plt.subplots(figsize=(x,y))
#ax.ticklabel_format(scilimits=[-1,1])
plt.loglog(Jget, Ji, '.' ,color = "steelblue")
plt.loglog(abs(Jget), abs(Ji), '.',color = "steelblue")
plt.loglog([1.e-80, 1.], [1.e-80, 1.], color = "orange", label = "1to1 line")
plt.grid("True")
plt.legend()
plt.xlim([10E-90,10E0])
plt.ylim([10E-90,10E0])
plt.xlabel("GETELEC1.0 current densities (A$nm^{-2}$)")
plt.ylabel("GETELEC2.0 current densities (A$nm^{-2}$)")
plt.title("Comparing current from GETELEC1.0 & 2.0")
plt.savefig("J comparison.png")
#plt.show()

fig2, ax = plt.subplots(figsize=(x,y))
plt.loglog(Pget, Pi, '.',color = "steelblue")
plt.loglog(abs(Pget), abs(Pi), '.',color = "steelblue")
plt.loglog([1.e-80, 1.], [1.e-80, 1.],color = "orange", label = "1to1 line")
plt.grid("True")
plt.legend()
plt.xlim([10E-90,10E0])
plt.ylim([10E-90,10E0])
plt.xlabel("GETELEC1.0 Nottigham heat (W$nm^{-2}$)")
plt.ylabel("GETELEC2.0 Nottigham heat (W$nm^{-2}$)")
plt.title("Comparing heat from GETELEC1.0 & 2.0")
plt.savefig("Pn comparison.png")
#plt.show()
"""
#endregion

# region Multiple data point calculation routine - semiconductor
"""
Npoly = 5
NGi = 128
zs = 1.6183e-4
kBoltz = 8.6173324e-5 

tab = gtab.Tabulator()

Ef = 4.75
Eg = 0.661
Temp = 300.
m = 9.1093837015e-31 
me = 1.64*m # effective electron mass @300K, https://doi.org/10.1142/S0219749909005833
mp = 0.28*m # effective hole mass @300k, https://doi.org/10.1142/S0219749909005833
kT = kBoltz * Temp

semiconductor_emitter = gtab.Semiconductor_Emitter(tab)

semiconductor_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
semiconductor_emitter.emitter.Interpolate_Gammow()

metal_emitter = gtab.Metal_Emitter(tab)

metal_emitter.emitter.Define_Barrier_Parameters(5., 10., 10.)
metal_emitter.emitter.Interpolate_Gammow()

resolution = 1000

j_c = np.zeros(resolution) 
j_v = np.zeros(resolution) 
j_total = np.zeros(resolution) 
j_c2 = np.zeros(resolution) 
j_v2 = np.zeros(resolution) 
j_total2 = np.zeros(resolution) 
pn_c = np.zeros(resolution) 
pn_v = np.zeros(resolution) 
pn_total = np.zeros(resolution) 
Bottom_Ec = np.zeros(resolution) 
j_metal = np.zeros(resolution) 
pn_metal = np.zeros(resolution)
c_rel_error = np.zeros(resolution)
v_rel_error = np.zeros(resolution)
total_rel_error = np.zeros(resolution)

for i in range(resolution):
    Ec = i*(1.3/resolution)+3.9

    semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec, Ef, Eg, kT, m, me, mp)

    metal_emitter.Define_Emitter_Parameters(Ec, kT)


    j_c[i], j_v[i], j_total[i] = semiconductor_emitter.Current_Density_from_Semiconductors()
    
    j_c2[i], j_v2[i], j_total2[i] = semiconductor_emitter.Distributed_Current_Density_from_Semiconductors()
    
    pn_c[i], pn_v[i], pn_total[i] = semiconductor_emitter.Nottingham_Heat_from_Semiconductors()
    
    energy_c, distribution_c, energy_v, distribution_v = semiconductor_emitter.Energy_Distribution_from_Semiconductors()
    
    
    j_metal[i] = metal_emitter.Current_Density()
    
    pn_metal[i] = metal_emitter.Nottingham_Heat()
    
    energy_space_metal, distribution_metal = metal_emitter.Energy_Distribution()
    
    
    Bottom_Ec[i] = -Ec
    
    
    c_abs_error = abs(j_c2[i]-j_c[i])
    c_rel_error[i] = c_abs_error/j_c[i]

    v_abs_error = abs(j_v2[i]-j_v[i])
    v_rel_error[i] = v_abs_error/j_v[i]

    total_abs_error = abs(j_total2[i]-j_total[i])
    total_rel_error[i] = total_abs_error/j_total[i]
    
    
    title = "Normalised energy distrution, Ec = " + str(Ec)
    
    save_png = "Normalised energy distrution, Ec = " + str(Ec) + ".png"
    
    #fig = plt.figure(figsize=(16,6))
    #plt.plot(energy_space_metal, distribution_metal/max(distribution_metal))
    #plt.plot(energy_c, distribution_c/max(distribution_c))
    #plt.plot(energy_v, distribution_v/max(distribution_v))
    #plt.grid("True")
    #plt.title(title)
    #plt.savefig(save_png)

c_mean = np.mean(c_rel_error)
c_std = np.std(c_rel_error)
c_max = max(c_rel_error)

v_mean = np.mean(v_rel_error)
v_std = np.std(v_rel_error)
v_max = max(v_rel_error)


total_mean = np.mean(total_rel_error)
total_std = np.std(total_rel_error)
total_max = max(total_rel_error)

font = 60
x = 40
y = 17

mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font+5
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 5

fig, ax = plt.subplots(figsize=(x,y))
ax.ticklabel_format(scilimits=[-1,1])
plt.semilogy(Bottom_Ec, j_c, color = "steelblue", label = "electrons from conduction band")
plt.semilogy(Bottom_Ec, j_v, color = "green", label = "electrons from valence band")
plt.semilogy(Bottom_Ec, j_total, color = "orange", label = "total emitted electrons")
#plt.plot(Bottom_Ec, j_metal)
#plt.yscale("symlog", linscale=-1)
plt.legend()
plt.grid("True")
plt.xlabel("Bottom of conduction band (eV)")
plt.ylabel("$J_{e}$ (A$nm^{-2}$)")
plt.title("Current density as function of band bending")
plt.savefig("Detail.png")
plt.savefig("Detail.svg")

file = open("Field_Emission_Data.txt", "w+")

file.write("Simulation temperature = %f\r\n\r\n" % Temp)

for i in range(len(Bottom_Ec)):
    file.write("\r\nSemiconductor Ec = %d" % (Bottom_Ec[i]))
    file.write(" Metal Ef = %d\r\n" % (Bottom_Ec[i]))
    file.write("    Current from semiconductor %e" % (j_total[i]))
    file.write("    Current from Metal = %e\r\n" % (j_metal[i]))
    file.write("        Semiconductor current relative error = %e\r\n" % (total_rel_error[i]))
    file.write("    Heat from semiconductor = %e" % (pn_total[i]))
    file.write("    Heat from metal = %e\r\n" % (pn_metal[i]))
    
file.write("\r\nMax relative error = %e\r\n" % (total_max))
  
file.close()
"""
# endregion

# region
#comsol_data =  np.array([-3.704766376545930000E+00,-3.710674031866700000E+00,-3.720906416070840000E+00,-3.732116601817530000E+00,-3.740652716975520000E+00,-3.751598406358010000E+00,-3.761510201146620000E+00,-3.770604439815580000E+00,-3.780260163707000000E+00,-3.791351630345240000E+00,-3.800263660657250000E+00,-3.811481390616470000E+00,-3.822632072265170000E+00,-3.830607338695320000E+00,-3.840082463019220000E+00,-3.850422700845510000E+00,-3.860524243732050000E+00,-3.870099020687190000E+00,-3.880047324080420000E+00,-3.890040198098470000E+00,-3.900573044943670000E+00,-3.904948391402420000E+00])
"""
Npoly = 5
NGi = 128
zs = 1.6183e-4
kBoltz = 8.6173324e-5 

tab = gtab.Tabulator()
Ec = np.array([5.07E+00,5.05E+00,4.96E+00,4.91E+00,4.89E+00])
Ef = np.array([4.783503100772620000E+00,4.781818198377520000E+00,4.778199674411920000E+00,4.772537667287590000E+00,4.766705287842640000E+00])
Eg = 1.12
Temp = np.array([5E+02,10E+02,4.79E+02,4.61E+02,4.46E+02])
Field = np.array([6.24E+00,5.69E+00,3.63E+00,2.68E+00,2.29E+00])
m = 9.1093837015e-31 
me = 1.64*m # effective electron mass @300K, https://doi.org/10.1142/S0219749909005833
mp = 0.28*m # effective hole mass @300k, https://doi.org/10.1142/S0219749909005833
kT = kBoltz * Temp

resolution = len(Field)

j_c = np.zeros(resolution) 
j_v = np.zeros(resolution) 
j_total = np.zeros(resolution) 
j_c2 = np.zeros(resolution) 
j_v2 = np.zeros(resolution) 
j_total2 = np.zeros(resolution) 
pn_c = np.zeros(resolution) 
pn_v = np.zeros(resolution) 
pn_total = np.zeros(resolution) 
Bottom_Ec = np.zeros(resolution) 
j_metal = np.zeros(resolution) 
pn_metal = np.zeros(resolution)
c_rel_error = np.zeros(resolution)
v_rel_error = np.zeros(resolution)
total_rel_error = np.zeros(resolution)



semiconductor_emitter = gtab.Semiconductor_Emitter(tab)

semiconductor_emitter.emitter.Define_Barrier_Parameters(Field[0], 20., 10.)
semiconductor_emitter.emitter.Interpolate_Gammow()
semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec[0], Ef[0], Eg, kT[0], m, me, mp)
energy_c0, distribution_c0, energy_v0, distribution_v0 = semiconductor_emitter.Energy_Distribution_from_Semiconductors()

semiconductor_emitter.emitter.Define_Barrier_Parameters(Field[2], 20., 10.)
semiconductor_emitter.emitter.Interpolate_Gammow()
semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec[2], Ef[2], Eg, kT[2], m, me, mp)
energy_c2, distribution_c2, energy_v2, distribution_v2 = semiconductor_emitter.Energy_Distribution_from_Semiconductors()

semiconductor_emitter.emitter.Define_Barrier_Parameters(Field[4], 20., 10.)
semiconductor_emitter.emitter.Interpolate_Gammow()
semiconductor_emitter.Define_Semiconductor_Emitter_Parameters(Ec[4], Ef[4], Eg, kT[4], m, me, mp)
energy_c4, distribution_c4, energy_v4, distribution_v4 = semiconductor_emitter.Energy_Distribution_from_Semiconductors()

font = 60
x = 40
y = 17

mb.rcParams["font.size"] = font
mb.rcParams["axes.labelsize"] = font
mb.rcParams["xtick.labelsize"] = font
mb.rcParams["ytick.labelsize"] = font+5
mb.rcParams["legend.fontsize"] = font
mb.rcParams["lines.linewidth"] = 5


fig = plt.figure(figsize=(x,y))

plt.plot(energy_c0-Ef[0]/2, distribution_c0/max(distribution_c0), color = "steelblue", label = "emission from apex")
#plt.plot(energy_v0, distribution_v0/max(distribution_v0), color = "steelblue")
plt.plot(energy_c2-Ef[2]/2, distribution_c2/max(distribution_c2), color = "green", label = "emission from 1/4 arc")
#plt.plot(energy_v2, distribution_v2/max(distribution_v2), color = "green")
plt.plot(energy_c4-Ef[4]/2, distribution_c4/max(distribution_c4), color = "orange", label = "emission from half arc")
#plt.plot(energy_v4, distribution_v4/max(distribution_v4), color = "orange")

plt.legend()
plt.grid("True")
plt.xlabel("Energy (eV)")
plt.ylabel("$J_e$(E) (A$nm^{-2}$$eV^{-1}$)")
plt.title("Normalised electron energy distribution")
plt.savefig("Energy_Distribution.png")
plt.savefig("Energy_Distribution.svg")

#fig, ax = plt.subplots(figsize=(x,y))
#ax.ticklabel_format(scilimits=[-1,1])
#plt.semilogy(Bottom_Ec, j_c, color = "steelblue", label = "electrons from conduction band")
#plt.semilogy(Bottom_Ec, j_v, color = "green", label = "electrons from valence band")
#plt.semilogy(Bottom_Ec, j_total, color = "orange", label = "total emitted electrons")
##plt.plot(Bottom_Ec, j_metal)
#plt.yscale("symlog", linscale=-1)
#plt.legend()
#plt.grid("True")
#plt.xlabel("Bottom of conduction band (eV)")
#plt.ylabel("$J_{e}$ (A$nm^{-2}$)")
#plt.title("Current density as function of band bending")
#plt.savefig("Detail.png")
#plt.savefig("Detail.svg")

"""
# endregion


