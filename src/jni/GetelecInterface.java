public class GetelecInterface {
    // Load the shared library (Ensure libgetelec.so or getelec.dll is available)
    static {
        System.loadLibrary("getelec");
    }

    // Native methods matching the C++ functions
    private native long Getelec_new_with_config(String configPath, String barrierType);
    private native void Getelec_delete(long obj);
    private native void Getelec_setRadius(long obj, double[] radii);
    private native void Getelec_setField(long obj, double[] fields);
    private native void Getelec_setGamma(long obj, double[] gammas);
    private native void Getelec_setkT(long obj, double[] kTs);
    private native void Getelec_setWorkFunction(long obj, double[] workFunctions);
    private native void Getelec_setBandDepth(long obj, double[] bandDepths);
    private native void Getelec_setEffectiveMass(long obj, double[] effectiveMasses);
    private native int Getelec_run(long obj, boolean calculateSpectra);
    private native double[] Getelec_getCurrentDensities(long obj);
    private native double[] Getelec_getNottinghamHeats(long obj);
    private native double Getelec_calculateTransmissionForEnergy(long obj, double energy, int paramsIndex);
    private native double[] Getelec_calculateTransmissionForEnergies(long obj, double[] energies, int paramsIndex);
    private native double[] Getelec_calculateTransmissionForManyEnergies(long obj, double[] energies, int paramsIndex);
    private native double[] Getelec_getSpectraEnergies(long obj, int index);
    private native double[] Getelec_getSpectraValues(long obj, int index);
    private native double[] Getelec_getSpectraDerivatives(long obj, int index);
    private native void Getelec_getBarrierValues(long obj, double[] x, double[] potential, int paramsIndex);
    private native void Getelec_getBarrierIntegrationLimits(long obj, double[] xInitial, double[] xFinal, int paramsIndex);
    
    private long getelecObjPtr;

    public GetelecInterface(String configPath, String barrierType) {
        getelecObjPtr = Getelec_new_with_config(configPath, barrierType);
    }

    public void setRadius(double[] radii) { Getelec_setRadius(getelecObjPtr, radii); }
    public void setField(double[] fields) { Getelec_setField(getelecObjPtr, fields); }
    public void setGamma(double[] gammas) { Getelec_setGamma(getelecObjPtr, gammas); }
    public void setkTCopy(double[] kTs) { Getelec_setkT(getelecObjPtr, kTs); }
    public void setWorkFunction(double[] workFunctions) { Getelec_setWorkFunction(getelecObjPtr, workFunctions); }
    public void setBandDepth(double[] bandDepths) { Getelec_setBandDepth(getelecObjPtr, bandDepths); }
    public void setEffectiveMass(double[] effectiveMasses) { Getelec_setEffectiveMass(getelecObjPtr, effectiveMasses); }
    
    public int run(boolean calculateSpectra) { return Getelec_run(getelecObjPtr, calculateSpectra); }
    public double[] getCurrentDensities() { return Getelec_getCurrentDensities(getelecObjPtr); }
    public double[] getNottinghamHeats() { return Getelec_getNottinghamHeats(getelecObjPtr); }
    public double calculateTransmissionForEnergy(double energy, int paramsIndex) { return Getelec_calculateTransmissionForEnergy(getelecObjPtr, energy, paramsIndex); }
    public double[] calculateTransmissionForEnergies(double[] energies, int paramsIndex) { return Getelec_calculateTransmissionForEnergies(getelecObjPtr, energies, paramsIndex); }
    public double[] calculateTransmissionForManyEnergies(double[] energies, int paramsIndex) { return Getelec_calculateTransmissionForManyEnergies(getelecObjPtr, energies, paramsIndex); }
    public double[] getSpectraEnergies(int index) { return Getelec_getSpectraEnergies(getelecObjPtr, index); }
    public double[] getSpectraValues(int index) { return Getelec_getSpectraValues(getelecObjPtr, index); }
    public double[] getSpectraDerivatives(int index) { return Getelec_getSpectraDerivatives(getelecObjPtr, index); }
    public void getBarrierValues(double[] x, double[] potential, int paramsIndex) { Getelec_getBarrierValues(getelecObjPtr, x, potential, paramsIndex); }
    public void getBarrierIntegrationLimits(double[] xInitial, double[] xFinal, int paramsIndex) { Getelec_getBarrierIntegrationLimits(getelecObjPtr, xInitial, xFinal, paramsIndex); }
    
    public void close() { Getelec_delete(getelecObjPtr); }
}
