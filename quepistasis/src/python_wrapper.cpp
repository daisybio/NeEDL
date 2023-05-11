#include "../header/python_wrapper.h"

PythonWrapper::PythonWrapper() { }

PythonWrapper& PythonWrapper::get_instance() {
    static PythonWrapper* instance;

    if(instance == nullptr) {
        instance = new PythonWrapper();
        instance->init();
    }

    return *instance;
}

int PythonWrapper::init() {
    Py_Initialize();

    // allow file to be seen from Python
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");
    PyRun_SimpleString("sys.path.append(\"./quepistasis\")");

    // load modules and functions
    std::cout << "Load QA module" << std::endl << std::flush;
    this->pModuleQuantumAnnealer = get_module("dwave_device");
    if(this->pModuleQuantumAnnealer == NULL) { return -1; }

    std::cout << "Load AZURE module" << std::endl << std::flush;
    this->pModuleAzureOptimizer = get_module("azure_qio_device");
    if(this->pModuleAzureOptimizer == NULL) { return -1; }

    std::cout << "Load QC module" << std::endl << std::flush;
    this->pModuleQAOA = get_module("azure_device");
    if(this->pModuleQAOA == NULL) { return -1; }

    std::cout << "Load QA function" << std::endl << std::flush;
    this->pFuncQuantumAnnealer = get_function("run_quantum_annealer", pModuleQuantumAnnealer);
    if(this->pFuncQuantumAnnealer == NULL) { return -1; }

    std::cout << "Load AZURE function" << std::endl << std::flush;
    this->pFuncAzureOptimizer = get_function("run_azure_qio", pModuleAzureOptimizer);
    if(this->pFuncAzureOptimizer == NULL) { return -1; }

    std::cout << "Load QC function" << std::endl << std::flush;
    this->pFuncQAOA = get_function("run_qaoa", pModuleQAOA);
    if(this->pFuncQAOA == NULL) { return -1; }

    std::cout << "Loaded" << std::endl << std::flush;

    return 0;
}

void PythonWrapper::close() {
    Py_DECREF(this->pFuncQAOA);
    Py_DECREF(this->pFuncAzureOptimizer);
    Py_DECREF(this->pFuncQuantumAnnealer);
    Py_DECREF(this->pModuleQAOA);
    Py_DECREF(this->pModuleAzureOptimizer);
    Py_DECREF(this->pModuleQuantumAnnealer);
    std::cout << "Finalizing..." << std::endl << std::flush;
    if(Py_FinalizeEx() != 0) { // Python >= 3.6
        PyErr_Print();
        fprintf(stderr, "Py_FinalizeEx failes\n");
    };
}

PyObject* PythonWrapper::get_module(const char* filepath) {

    PyObject *pName = PyUnicode_FromString(filepath); // represent name string in python format
    PyObject *pModule = PyImport_Import(pName); // load module

    if(pModule == NULL) {
        PyErr_Print();
        fprintf(stderr, "loading module have failed\n");
        return NULL;
    }

    Py_DECREF(pName); // deallocate name variable
    return pModule;
}

PyObject* PythonWrapper::get_function(const char* functionName, PyObject *pModule) {

    if(pModule == NULL) {
        fprintf(stderr, "Call load_module first\n");
        return NULL;
    }

    PyObject* pFunc = PyObject_GetAttrString(pModule, functionName);

    if(pFunc == NULL) {
        PyErr_Print();
        fprintf(stderr, "pFunc failes\n");
        return NULL;
    }

    if(!PyCallable_Check(pFunc)) {
        fprintf(stderr, "pFunc is not callable\n");
        return NULL;
    }

    return pFunc;
}

PyObject* PythonWrapper::call_function(PyObject *pFunc, PyObject *pArgs) {

    if(!PyCallable_Check(pFunc)) {
        fprintf(stderr, "pFunc is not callable\n");
        return NULL;
    }

    if(!PyTuple_Check(pArgs)) {
        fprintf(stderr, "pArgs is not a function argument\n");
        return NULL;
    }

    PyObject *pValue = PyObject_CallObject(pFunc, pArgs);

    if(pValue == NULL) {
        PyErr_Print();
        fprintf(stderr, "pValue failes\n");
        return NULL;
    }

    return pValue;
}

PyObject* PythonWrapper::create_h(std::vector<double> h, std::vector<PyObject*>& garbage) {
    PyObject *pHList = PyList_New(h.size());
    for(size_t i = 0; i < h.size(); i++) {
        PyObject* pHCoeff = PyFloat_FromDouble(h[i]);
        PyList_SetItem(pHList, i, pHCoeff);
        garbage.push_back(pHCoeff);
    }
    return pHList;
}

PyObject* PythonWrapper::create_j(std::vector<double> J, std::vector<int> start, std::vector<int> end, std::vector<PyObject*>& garbage) {
    PyObject *pJDict = PyDict_New();
    for(size_t i = 0; i < J.size(); i++) {
        PyObject* pValue = PyFloat_FromDouble(J[i]);
        PyObject* pStart = PyLong_FromLong(start[i]);
        PyObject* pEnd = PyLong_FromLong(end[i]);
        PyObject* pKey = PyTuple_New(2);
        PyTuple_SetItem(pKey, 0, pStart);
        PyTuple_SetItem(pKey, 1, pEnd);
        PyDict_SetItem(pJDict, pKey, pValue);
        garbage.push_back(pValue);
        garbage.push_back(pStart);
        garbage.push_back(pEnd);
        garbage.push_back(pKey);
    }
    return pJDict;
}

std::vector<int> PythonWrapper::unpack_result(PyObject* pResult) {

    std::vector<int> spins;

    if(NULL != pResult && !PyList_Check(pResult)) {
        fprintf(stderr, "pResult is not a list\n");
        return spins;
    }

    size_t size = PyList_Size(pResult);

    if(size == 0) {
        // has an error while trying to find an embedding for the problem
        return spins;
    }

    for(size_t i = 0; i < size; i++) {
        PyObject *item = PyList_GetItem(pResult, i);
        long the_spin = PyLong_AsLong(item);
        // Py_DECREF(item);
        spins.push_back(the_spin);
    }

    return spins;    
}

std::vector<int> PythonWrapper::run_qaoa(
        const char* vendor,
        const char* azure_subscription_id, const char* azure_resource_group, const char* azure_name, const char* azure_location,
        const char* azure_backend,
        std::vector<double> h, std::vector<double> J, std::vector<int> start, std::vector<int> end,
        const char* optimizer, int maxiter, int reps, int n_shots, int is_recursive_qaoa, const char* save_path) {

    std::vector<PyObject*> garbage;

    // create parameters
    PyObject* pVendor = PyUnicode_FromString(vendor);
    PyObject* pAzureSubscriptionId = PyUnicode_FromString(azure_subscription_id);
    PyObject* pAzureResourceGroup = PyUnicode_FromString(azure_resource_group);
    PyObject* pAzureName = PyUnicode_FromString(azure_name);
    PyObject* pAzureLocation = PyUnicode_FromString(azure_location);
    PyObject* pAzureBackend = PyUnicode_FromString(azure_backend);
    PyObject* pH = this->create_h(h, garbage);
    PyObject* pJ = this->create_j(J, start, end, garbage);
    PyObject* pOptimizer = PyUnicode_FromString(optimizer);
    PyObject* pMaxIters = PyLong_FromLong(maxiter);
    PyObject* pReps = PyLong_FromLong(reps);
    PyObject* pNShots = PyLong_FromLong(n_shots);
    PyObject* pIsRecursiveQaoa = PyLong_FromLong(is_recursive_qaoa);
    PyObject* pSavePath = PyUnicode_FromString(save_path);

    PyObject *pArgs = PyTuple_New(14);
    PyTuple_SetItem(pArgs, 0, pVendor);
    PyTuple_SetItem(pArgs, 1, pAzureSubscriptionId);
    PyTuple_SetItem(pArgs, 2, pAzureResourceGroup);
    PyTuple_SetItem(pArgs, 3, pAzureName);
    PyTuple_SetItem(pArgs, 4, pAzureLocation);
    PyTuple_SetItem(pArgs, 5, pAzureBackend);
    PyTuple_SetItem(pArgs, 6, pH);
    PyTuple_SetItem(pArgs, 7, pJ);
    PyTuple_SetItem(pArgs, 8, pOptimizer);
    PyTuple_SetItem(pArgs, 9, pMaxIters);
    PyTuple_SetItem(pArgs,10, pReps);
    PyTuple_SetItem(pArgs,11, pNShots);
    PyTuple_SetItem(pArgs,12, pIsRecursiveQaoa);
    PyTuple_SetItem(pArgs,13, pSavePath);

    // call Python
    PyObject *pResult = this->call_function(this->pFuncQAOA, pArgs);

    // convert Python result in C format
    std::vector<int> result = this->unpack_result(pResult);

    // clean parameters and results object
    // for(PyObject* pObj : garbage) {
    //     Py_DECREF(pObj);
    // }
    // Py_DECREF(pVendor);
    // Py_DECREF(pAzureSubscriptionId);
    // Py_DECREF(pAzureName);
    // Py_DECREF(pAzureResourceGroup);
    // Py_DECREF(pAzureLocation);
    // Py_DECREF(pAzureBackend);
    // Py_DECREF(pH);
    // Py_DECREF(pJ);
    // Py_DECREF(pOptimizer);
    // Py_DECREF(pMaxIters);
    // Py_DECREF(pReps);
    // Py_DECREF(pNShots);
    // Py_DECREF(pIsRecursiveQaoa);
    // Py_DECREF(pSavePath);
    // Py_DECREF(pArgs);
    // Py_DECREF(pResult);

    return result;
}


std::vector<int> PythonWrapper::run_azure_optimizers(
    std::vector<double> h, std::vector<double> J, std::vector<int> start, std::vector<int> end,
    const char* azure_subscription_id, const char* azure_resource_group, const char* azure_name, const char* azure_location,
    int timeout_seconds, int seed, const char* save_path) {

    std::vector<PyObject*> garbage;

    // create parameters
    PyObject* pH = this->create_h(h, garbage);
    PyObject* pJ = this->create_j(J, start, end, garbage);
    PyObject* pAzureSubscriptionId = PyUnicode_FromString(azure_subscription_id);
    PyObject* pAzureResourceGroup = PyUnicode_FromString(azure_resource_group);
    PyObject* pAzureName = PyUnicode_FromString(azure_name);
    PyObject* pAzureLocation = PyUnicode_FromString(azure_location);
    PyObject* pTimeout = PyLong_FromLong(timeout_seconds);
    PyObject* pSeed = PyLong_FromLong(seed);
    PyObject* pSavePath = PyUnicode_FromString(save_path);

    PyObject *pArgs = PyTuple_New(9);
    PyTuple_SetItem(pArgs, 0, pH);
    PyTuple_SetItem(pArgs, 1, pJ);
    PyTuple_SetItem(pArgs, 2, pAzureSubscriptionId);
    PyTuple_SetItem(pArgs, 3, pAzureResourceGroup);
    PyTuple_SetItem(pArgs, 4, pAzureName);
    PyTuple_SetItem(pArgs, 5, pAzureLocation);
    PyTuple_SetItem(pArgs, 6, pTimeout);
    PyTuple_SetItem(pArgs, 7, pSeed);
    PyTuple_SetItem(pArgs, 8, pSavePath);

    // call Python
    PyObject *pResult = this->call_function(this->pFuncAzureOptimizer, pArgs);

    // convert Python result in C format
    std::vector<int> result = this->unpack_result(pResult);

    // clean parameters and results object -> Cleaning vars may lead to SEG FAULT...
    // Py_XDECREF(pArgs);
    // Py_XDECREF(pSavePath);
    // Py_XDECREF(pSeed);
    // Py_XDECREF(pTimeout);
    // Py_XDECREF(pAzureLocation);
    // Py_XDECREF(pAzureName);
    // Py_XDECREF(pAzureResourceGroup);
    // Py_XDECREF(pAzureSubscriptionId);
    // Py_XDECREF(pJ);
    // Py_XDECREF(pH);
    // for(PyObject* pObj : garbage) {
    //     Py_XDECREF(pObj);
    // }

    return result;
}

std::vector<int> PythonWrapper::run_quantum_annealer(
    const char* token, std::vector<double> h, std::vector<double> J, std::vector<int> start, std::vector<int> end,
    int num_reads, int solver_idx, double fw_annealing_ramp_time, double fw_annealing_pause_time,
    double rev_annealing_ramp_time, double rev_annealing_pause_time, double rev_annealing_s_target,
    const char* save_path) {

    std::vector<PyObject*> garbage;

    // create parameters
    PyObject* pH = this->create_h(h, garbage);
    PyObject* pJ = this->create_j(J, start, end, garbage);
    PyObject* pToken = PyUnicode_FromString(token);
    PyObject* pNumReads = PyLong_FromLong(num_reads);
    PyObject* pSolverIdx = PyLong_FromLong(solver_idx);
    PyObject* pFwAnnealingRampTime = PyFloat_FromDouble(fw_annealing_ramp_time);
    PyObject* pFwAnnealingPauseTime = PyFloat_FromDouble(fw_annealing_pause_time);
    PyObject* pRevAnnealingRampTime = PyFloat_FromDouble(rev_annealing_ramp_time);
    PyObject* pRevAnnealingPauseTime = PyFloat_FromDouble(rev_annealing_pause_time);
    PyObject* pRevAnnealingSTarget = PyFloat_FromDouble(rev_annealing_s_target);
    PyObject* pSavePath = PyUnicode_FromString(save_path);

    PyObject *pArgs = PyTuple_New(11);
    PyTuple_SetItem(pArgs, 0, pToken);
    PyTuple_SetItem(pArgs, 1, pH);
    PyTuple_SetItem(pArgs, 2, pJ);
    PyTuple_SetItem(pArgs, 3, pNumReads);
    PyTuple_SetItem(pArgs, 4, pSolverIdx);
    PyTuple_SetItem(pArgs, 5, pFwAnnealingRampTime);
    PyTuple_SetItem(pArgs, 6, pFwAnnealingPauseTime);
    PyTuple_SetItem(pArgs, 7, pRevAnnealingRampTime);
    PyTuple_SetItem(pArgs, 8, pRevAnnealingPauseTime);
    PyTuple_SetItem(pArgs, 9, pRevAnnealingSTarget);
    PyTuple_SetItem(pArgs, 10, pSavePath);

    // call Python
    PyObject *pResult = this->call_function(this->pFuncQuantumAnnealer, pArgs);

    // convert Python result in C format
    std::vector<int> result = this->unpack_result(pResult);

    // clean parameters and results object
    // for(PyObject* pObj : garbage) {
    //     Py_DECREF(pObj);
    // }
    // Py_DECREF(pToken);
    // Py_DECREF(pH);
    // Py_DECREF(pJ);
    // Py_DECREF(pNumReads);
    // Py_DECREF(pSolverIdx);
    // Py_DECREF(pFwAnnealingRampTime);
    // Py_DECREF(pFwAnnealingPauseTime);
    // Py_DECREF(pRevAnnealingRampTime);
    // Py_DECREF(pRevAnnealingPauseTime);
    // Py_DECREF(pRevAnnealingSTarget);
    // Py_DECREF(pSavePath);
    // Py_DECREF(pArgs);
    // Py_DECREF(pResult);

    return result;
}

std::vector<int> PythonWrapper::run_sqa(
    std::vector<double> h, std::vector<double> J, std::vector<int> start, std::vector<int> end,
    int clique_size, int shots, const char* save_path) {

    return std::vector<int>();
}