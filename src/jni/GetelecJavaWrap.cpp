#include "Getelec.h"
#include "GetelecInterface.h"
#include <vector>

extern "C" {

JNIEXPORT jlong JNICALL Java_GetelecInterface_Getelec_1new_1with_1config
(JNIEnv* env, jobject obj, jstring configFileName, jstring barrierType) {
    const char* config = env->GetStringUTFChars(configFileName, nullptr);
    const char* barrier = env->GetStringUTFChars(barrierType, nullptr);
    getelec::Getelec* instance = new getelec::Getelec(config, barrier);
    env->ReleaseStringUTFChars(configFileName, config);
    env->ReleaseStringUTFChars(barrierType, barrier);
    return (jlong)instance;
}

// Destructor
JNIEXPORT void JNICALL Java_GetelecInterface_Getelec_1delete
(JNIEnv* env, jobject obj, jlong instancePtr) {
    delete (getelec::Getelec*)instancePtr;
}

JNIEXPORT jint JNICALL Java_GetelecInterface_Getelec_1run
(JNIEnv* env, jobject obj, jlong instancePtr, jboolean calculateSpectra) {
    return ((getelec::Getelec*)instancePtr)->run(getelec::CalculationFlags::CurrentDensity | getelec::CalculationFlags::TotalEnergyDistribution | getelec::CalculationFlags::NottinghamHeat);
}

JNIEXPORT void JNICALL Java_GetelecInterface_Getelec_1setField
(JNIEnv* env, jobject obj, jlong instancePtr, jdoubleArray fieldArray) {
    
    // Get array length
    jsize length = env->GetArrayLength(fieldArray);

    // Get direct pointer to Java array (may pin memory, avoids copy)
    jboolean isCopy;
    void* rawData = env->GetPrimitiveArrayCritical(fieldArray, &isCopy);
    jdouble* data = static_cast<jdouble*>(rawData);  // ✅ Explicit cast

    if (!data) throw std::runtime_error("Cannot get primitive array critical");

    ((Getelec*)instancePtr)->setField(data, length);

    // Release the array
    env->ReleasePrimitiveArrayCritical(fieldArray, rawData, JNI_ABORT);
}

JNIEXPORT void JNICALL Java_GetelecInterface_Getelec_1setRadius
(JNIEnv* env, jobject obj, jlong instancePtr, jdoubleArray radiusArray) {
    
    // Get array length
    jsize length = env->GetArrayLength(radiusArray);

    // Get direct pointer to Java array (may pin memory, avoids copy)
    jboolean isCopy;
    void* rawData = env->GetPrimitiveArrayCritical(radiusArray, &isCopy);
    jdouble* data = static_cast<jdouble*>(rawData);  // ✅ Explicit cast

    if (!data) throw std::runtime_error("Cannot get primitive array critical");

    ((getelec::Getelec*)instancePtr)->setRadius(data, length);

    // Release the array
    env->ReleasePrimitiveArrayCritical(radiusArray, rawData, JNI_ABORT);
}

JNIEXPORT void JNICALL Java_GetelecInterface_Getelec_1setGamma
(JNIEnv* env, jobject obj, jlong instancePtr, jdoubleArray gammaArray) {
    
    // Get array length
    jsize length = env->GetArrayLength(gammaArray);

    // Get direct pointer to Java array (may pin memory, avoids copy)
    jboolean isCopy;
    void* rawData = env->GetPrimitiveArrayCritical(gammaArray, &isCopy);
    jdouble* data = static_cast<jdouble*>(rawData);  // ✅ Explicit cast

    if (!data) throw std::runtime_error("Cannot get primitive array critical");

    ((Getelec*)instancePtr)->setGamma(data, length);

    // Release the array
    env->ReleasePrimitiveArrayCritical(gammaArray, rawData, JNI_ABORT);
}

JNIEXPORT void JNICALL Java_GetelecInterface_Getelec_1setkT
(JNIEnv* env, jobject obj, jlong instancePtr, jdoubleArray kTArray) {
    
    // Get array length
    jsize length = env->GetArrayLength(kTArray);

    // Get direct pointer to Java array (may pin memory, avoids copy)
    jboolean isCopy;
    void* rawData = env->GetPrimitiveArrayCritical(kTArray, &isCopy);
    jdouble* data = static_cast<jdouble*>(rawData);  // ✅ Explicit cast

    if (!data) throw std::runtime_error("Cannot get primitive array critical");

    ((Getelec*)instancePtr)->setkT(data, length);

    // Release the array
    env->ReleasePrimitiveArrayCritical(kTArray, rawData, JNI_ABORT);
}


JNIEXPORT void JNICALL Java_GetelecInterface_Getelec_1setWorkFunction
(JNIEnv* env, jobject obj, jlong instancePtr, jdoubleArray workFunctionArray) {
    
    // Get array length
    jsize length = env->GetArrayLength(workFunctionArray);

    // Get direct pointer to Java array (may pin memory, avoids copy)
    jboolean isCopy;
    void* rawData = env->GetPrimitiveArrayCritical(workFunctionArray, &isCopy);
    jdouble* data = static_cast<jdouble*>(rawData);  // ✅ Explicit cast

    if (!data) throw std::runtime_error("Cannot get primitive array critical");

    ((Getelec*)instancePtr)->setWorkFunction(data, length);

    // Release the array
    env->ReleasePrimitiveArrayCritical(workFunctionArray, rawData, JNI_ABORT);
}

JNIEXPORT void JNICALL Java_GetelecInterface_Getelec_1setBandDepth
(JNIEnv* env, jobject obj, jlong instancePtr, jdoubleArray bandDepthArray) {
    
    // Get array length
    jsize length = env->GetArrayLength(bandDepthArray);

    // Get direct pointer to Java array (may pin memory, avoids copy)
    jboolean isCopy;
    void* rawData = env->GetPrimitiveArrayCritical(bandDepthArray, &isCopy);
    jdouble* data = static_cast<jdouble*>(rawData);  // ✅ Explicit cast

    if (!data) throw std::runtime_error("Cannot get primitive array critical");

    ((Getelec*)instancePtr)->setBandDepth(data, length);

    // Release the array
    env->ReleasePrimitiveArrayCritical(bandDepthArray, rawData, JNI_ABORT);
}

JNIEXPORT void JNICALL Java_GetelecInterface_Getelec_1setEffectiveMass
(JNIEnv* env, jobject obj, jlong instancePtr, jdoubleArray effectiveMassArray) {
    
    // Get array length
    jsize length = env->GetArrayLength(effectiveMassArray);

    // Get direct pointer to Java array (may pin memory, avoids copy)
    jboolean isCopy;
    void* rawData = env->GetPrimitiveArrayCritical(effectiveMassArray, &isCopy);
    jdouble* data = static_cast<jdouble*>(rawData);  // ✅ Explicit cast

    if (!data) throw std::runtime_error("Cannot get primitive array critical");

    ((getelec::Getelec*)instancePtr)->setEffectiveMass(data, length);

    // Release the array
    env->ReleasePrimitiveArrayCritical(effectiveMassArray, rawData, JNI_ABORT);
}



JNIEXPORT jdoubleArray JNICALL Java_GetelecInterface_Getelec_1getCurrentDensities
(JNIEnv* env, jobject obj, jlong instancePtr) {
    size_t size;
    const double* densities = ((getelec::Getelec*)instancePtr)->getCurrentDensities(&size);

    jdoubleArray result = env->NewDoubleArray(size);
    env->SetDoubleArrayRegion(result, 0, size, densities);

    return result;
}

JNIEXPORT jdoubleArray JNICALL Java_GetelecInterface_Getelec_1getNottinghamHeats
(JNIEnv* env, jobject obj, jlong instancePtr) {
    size_t size;
    const double* heats = ((getelec::Getelec*)instancePtr)->getNottinghamHeats(&size);

    jdoubleArray result = env->NewDoubleArray(size);
    env->SetDoubleArrayRegion(result, 0, size, heats);

    return result;
}

JNIEXPORT jdouble JNICALL Java_GetelecInterface_Getelec_1calculateTransmissionForEnergy
(JNIEnv* env, jobject obj, jlong instancePtr, jdouble energy, jint paramsIndex) {
    return ((getelec::Getelec*)instancePtr)->calculateTransmissionCoefficientForEnergy(energy, paramsIndex);
}

JNIEXPORT jdoubleArray JNICALL Java_GetelecInterface_Getelec_1calculateTransmissionForEnergies
(JNIEnv* env, jobject obj, jlong instancePtr, jdoubleArray energyArray, jint paramsIndex) {
    
    // Get array length
    jsize length = env->GetArrayLength(energyArray);

    // Get direct pointer to Java array (may pin memory, avoids copy)
    jboolean isCopy;
    void* rawData = env->GetPrimitiveArrayCritical(energyArray, &isCopy);
    jdouble* data = static_cast<jdouble*>(rawData);  // ✅ Explicit cast

    if (!data) throw std::runtime_error("Cannot get primitive array critical");

    std::vector<double> energies(data, data + length);
    std::vector<double> result = ((getelec::Getelec*)instancePtr)->calculateTransmissionCoefficientForEnergies(energies, paramsIndex);

    // Release the array
    env->ReleasePrimitiveArrayCritical(energyArray, rawData, JNI_ABORT);

    jdoubleArray resultArray = env->NewDoubleArray(result.size());
    env->SetDoubleArrayRegion(resultArray, 0, result.size(), result.data());

    return resultArray;
}

JNIEXPORT jdoubleArray JNICALL Java_GetelecInterface_Getelec_1calculateTransmissionForManyEnergies
(JNIEnv* env, jobject obj, jlong instancePtr, jdoubleArray energyArray, jint paramsIndex) {
    
    // Get array length
    jsize length = env->GetArrayLength(energyArray);

    // Get direct pointer to Java array (may pin memory, avoids copy)
    jboolean isCopy;
    void* rawData = env->GetPrimitiveArrayCritical(energyArray, &isCopy);
    jdouble* data = static_cast<jdouble*>(rawData);  // ✅ Explicit cast

    if (!data) throw std::runtime_error("Cannot get primitive array critical");

    std::vector<double> energies(data, data + length);
    std::vector<double> result = ((getelec::Getelec*)instancePtr)->calculateTransmissionCoefficientForManyEnergies(energies, paramsIndex);

    // Release the array
    env->ReleasePrimitiveArrayCritical(energyArray, rawData, JNI_ABORT);

    jdoubleArray resultArray = env->NewDoubleArray(result.size());
    env->SetDoubleArrayRegion(resultArray, 0, result.size(), result.data());

    return resultArray;
}

JNIEXPORT jdoubleArray JNICALL Java_GetelecInterface_Getelec_1getSpectraEnergies
(JNIEnv* env, jobject obj, jlong instancePtr, jint index) {
    size_t length;
    const double* energies = ((getelec::Getelec*)instancePtr)->getSpectraEnergies(index, &length, 'T');

    jdoubleArray result = env->NewDoubleArray(length);
    env->SetDoubleArrayRegion(result, 0, length, energies);

    return result;
}

JNIEXPORT jdoubleArray JNICALL Java_GetelecInterface_Getelec_1getSpectraValues
(JNIEnv* env, jobject obj, jlong instancePtr, jint index) {
    size_t length;
    const double* values = ((getelec::Getelec*)instancePtr)->getSpectraValues(index, &length, 'T');

    jdoubleArray result = env->NewDoubleArray(length);
    env->SetDoubleArrayRegion(result, 0, length, values);

    return result;
}

JNIEXPORT jdoubleArray JNICALL Java_GetelecInterface_Getelec_1getSpectraDerivatives
(JNIEnv* env, jobject obj, jlong instancePtr, jint index) {
    size_t length;
    const double* derivatives = ((getelec::Getelec*)instancePtr)->getSpectraValues(index, &length, 'D');

    jdoubleArray result = env->NewDoubleArray(length);
    env->SetDoubleArrayRegion(result, 0, length, derivatives);

    return result;
}

JNIEXPORT void JNICALL Java_GetelecInterface_Getelec_1getBarrierValues
(JNIEnv* env, jobject obj, jlong instancePtr, jdoubleArray xArray, jdoubleArray potentialArray, jint paramsIndex) {
    
    // Get array length
    jsize length = env->GetArrayLength(xArray);

    // Get direct pointer to Java array (may pin memory, avoids copy)
    jboolean isCopyX, isCopyPotential;
    void* rawDataX = env->GetPrimitiveArrayCritical(xArray, &isCopyX);
    void* rawDataPotential = env->GetPrimitiveArrayCritical(potentialArray, &isCopyPotential);
    jdouble* xData = static_cast<jdouble*>(rawDataX);  // ✅ Explicit cast
    jdouble* potentialData = static_cast<jdouble*>(rawDataPotential);  // ✅ Explicit cast

    if (!xData || !potentialData) throw std::runtime_error("Cannot get primitive array critical");

    ((getelec::Getelec*)instancePtr)->getBarrierValues(xData, potentialData, length, paramsIndex);

    // Release the arrays
    env->ReleasePrimitiveArrayCritical(xArray, rawDataX, JNI_ABORT);
    env->ReleasePrimitiveArrayCritical(potentialArray, rawDataPotential, 0);
}

JNIEXPORT void JNICALL Java_GetelecInterface_Getelec_1getBarrierIntegrationLimits
(JNIEnv* env, jobject obj, jlong instancePtr, jdoubleArray xInitial, jdoubleArray xFinal, jint paramsIndex) {
    
    // Get direct pointer to Java array (may pin memory, avoids copy)
    auto [xI, xF] = ((getelec::Getelec*)instancePtr)->getBarrierIntegrationLimits(paramsIndex);
    env->SetDoubleArrayRegion(xInitial, 0, 1, &xI);
    env->SetDoubleArrayRegion(xFinal, 0, 1, &xF);
}


  

} // extern "C"
