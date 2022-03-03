// Copyright (c) 2022 Ruben Vorderman
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdint.h>

#define COFFSET_MAX 0xFFFFFFFFFFFFULL
#define UOFFSET_MAX 0xFFFFU

typedef struct {
    PyObject_HEAD;
    uint64_t voffset;
} VirtualFileOffset;

static void
VirtualFileOffset_dealloc(VirtualFileOffset * self){
    Py_TYPE(self)->tp_free(self);
}

static int 
VirtualFileOffset__init__(VirtualFileOffset * self, PyObject *args, 
                          PyObject *kwargs) {
    PyObject * coffset_o = NULL;
    PyObject * uoffset_o = NULL; 
    char * keywords[] = {"coffset", "uoffset", NULL};
    const char *format = "O!|O!:VirtualFileOffset.__init__";
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, format, keywords, 
        coffset_o, &PyLong_Type,
        uoffset_o, &PyLong_Type)) {
        return -1;
    }
    uint64_t coffset = PyLong_AsUnsignedLongLong(coffset_o);
    uint64_t uoffset = 0;
    if (uoffset_o != NULL) {
        uoffset = PyLong_AsUnsignedLongLong(uoffset_o);
    }
    if (coffset > COFFSET_MAX) {
        PyErr_Format(PyExc_OverflowError, 
            "%ld is larger than maximum allowed coffset value %ld",
            coffset, COFFSET_MAX);
        return -1;
    }
    if (uoffset > UOFFSET_MAX) {
        PyErr_Format(PyExc_OverflowError, 
            "%ld is larger than maximum allowed uoffset value %ld",
            coffset, COFFSET_MAX);
        return -1;
    }
    uint64_t voffset = (coffset << 16) & uoffset;
    self->voffset = voffset;
    return 0;
}

static PyObject * 
VirtualFileOffset_coffset_get(VirtualFileOffset *self, void *closure) {
    return PyLong_FromUnsignedLongLong(self->voffset >> 16);
}

static PyObject * 
VirtualFileOffset_uoffset_get(VirtualFileOffset * self, void *closure) {
    return PyLong_FromUnsignedLongLong(self->voffset & UOFFSET_MAX);
}

static PyObject *
VirtualFileOffset__voffset_get(VirtualFileOffset *self, void *closure) {
    return PyLong_FromUnsignedLongLong(self->voffset);
}

static PyGetSetDef VirtualFileOffset_properties[] = {
    {"coffset", (getter)VirtualFileOffset_coffset_get, NULL,
     "Offset to the beginning of a BGZF block", NULL},
    {"uoffset", (getter)VirtualFileOffset_uoffset_get, NULL,
     "Offset inside the BGZF block.", NULL},
    {"_voffset", (getter)VirtualFileOffset__voffset_get, NULL,
     "The internal virtual file offset integer.", NULL},
    {NULL}
};

// METHODS
PyDoc_STRVAR(VirtualFileOffset_from_bytes_doc,
"Create a VirtualFileOffset from a bytes object");

#define VIRTUAL_FILEOFFSET_FROM_BYTES_METHODDEF    \
    {"from_bytes", (PyCFunction)(void(*)(void))VirtualFileOffset_from_bytes, \
     METH_O | METH_CLASS, VirtualFileOffset_from_bytes_doc}

static PyObject * 
VirtualFileOffset_from_bytes(PyTypeObject *cls, PyObject* b) {
    if (!PyBytes_CheckExact(b)) {
        PyErr_Format(PyExc_TypeError, "Expected bytes, got %s", 
                     Py_TYPE(b)->tp_name);
        return NULL;
    }
    if (!(PyBytes_GET_SIZE(b) == sizeof(uint64_t))) {
        PyErr_Format(PyExc_ValueError, "b must have a length of exactly %d", 
                     sizeof(uint64_t));
    }
    VirtualFileOffset *vfs = PyObject_NEW(VirtualFileOffset, cls);
    vfs->voffset = *(uint64_t *)(PyBytes_AS_STRING(b));
    return (PyObject *)vfs;
}

static PyMethodDef VirtualFileOffset_methods[] = {
     VIRTUAL_FILEOFFSET_FROM_BYTES_METHODDEF,
    {NULL}
};

static PyTypeObject VirtualFileOffset_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_bgzf.VirtualFileOffset",
    .tp_basicsize = sizeof(VirtualFileOffset),
    .tp_dealloc = (destructor)VirtualFileOffset_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_methods = VirtualFileOffset_methods,
    .tp_getset = VirtualFileOffset_properties,
    .tp_init = (initproc)VirtualFileOffset__init__,
    .tp_new = PyType_GenericNew,
};

static PyObject * VirtualFileOffset_FromUint64(uint64_t i) {
    VirtualFileOffset * vfo = PyObject_NEW(VirtualFileOffset, 
                                           &VirtualFileOffset_Type);
    vfo->voffset = i;
    return (PyObject *)vfo;
}

PyDoc_STRVAR(vfo_list_from_bytes_doc,
"vfo_list_from_bytes($module, data, /)\n"
"--\n"
"\n"
"Creates a list of VirtualFileOffset objects from a bytes object.\n"
"\n"
"  data\n"
"    A bytes object that contains virtual file offset integers\n"
);

static PyObject * 
vfo_list_from_bytes(PyObject *module, PyObject *data) {
    if (!PyBytes_CheckExact(data)) {
        PyErr_Format(
            PyExc_TypeError, "data must be a bytes object, got %s", 
            Py_TYPE(data)->tp_name);
        return NULL;
    }
    Py_ssize_t size = PyBytes_GET_SIZE(data);
    if (size % sizeof(uint64_t)) {
        PyErr_Format(
            PyExc_ValueError, 
            "data must have a length that is a multiple of %ld, got %ld", 
            sizeof(uint64_t), size
        );
        return NULL;
    }
    Py_ssize_t n = size / sizeof(uint64_t);
    uint64_t * vfo = (uint64_t *)PyBytes_AS_STRING(data);
    PyObject * vfo_list = PyList_New(n);
    Py_ssize_t i = 0;
    while (i < n) {
        PyList_SET_ITEM(vfo_list, i, VirtualFileOffset_FromUint64(vfo[i]));
        i += 1;
    }
    return vfo_list;
}

PyDoc_STRVAR(vfo_chunk_list_from_bytes_doc,
"vfo_list_from_bytes($module, data, /)\n"
"--\n"
"\n"
"Creates a list of paired VirtualFileOffset tuples (start, end)\n"
"from a bytes object.\n"
"\n"
"  data\n"
"    A bytes object that contains virtual file offset integers\n"
);

static PyObject * 
vfo_chunk_list_from_bytes(PyObject *module, PyObject *data) {
    if (!PyBytes_CheckExact(data)) {
        PyErr_Format(
            PyExc_TypeError, "data must be a bytes object, got %s", 
            Py_TYPE(data)->tp_name);
        return NULL;
    }
    Py_ssize_t item_size = 2 * sizeof(uint64_t);
    Py_ssize_t size = PyBytes_GET_SIZE(data);
    if (size % item_size) {
        PyErr_Format(
            PyExc_ValueError, 
            "data must have a length that is a multiple of %ld, got %ld", 
            item_size, size
        );
        return NULL;
    }
    Py_ssize_t n = size / item_size;
    uint64_t * vfo = (uint64_t *)PyBytes_AS_STRING(data);
    PyObject * vfo_chunk_list = PyList_New(n);
    Py_ssize_t i = 0;
    PyObject * chunk;
    PyObject * start;
    PyObject * end;
    while (i < n) {
        start = VirtualFileOffset_FromUint64(*vfo);
        vfo += 1;
        end = VirtualFileOffset_FromUint64(*vfo);
        vfo += 1;
        chunk = PyTuple_Pack(2, start, end);
        PyList_SET_ITEM(vfo_chunk_list, i, chunk);
        i += 1;
    }
    return vfo_chunk_list;
}

static PyMethodDef _bgzf_methods[] = {
    {"vfo_list_from_bytes", (PyCFunction)(void(*)(void))vfo_list_from_bytes,
    METH_O, vfo_list_from_bytes_doc},
    {"vfo_chunk_list_from_bytes", (PyCFunction)(void(*)(void))vfo_chunk_list_from_bytes,
    METH_O, vfo_chunk_list_from_bytes_doc},
    {NULL}
};

static struct PyModuleDef _bgzf_module = {
    PyModuleDef_HEAD_INIT,
    "_bgzf",
    NULL,
    -1, 
    _bgzf_methods
};

PyMODINIT_FUNC
PyInit__bgzf(void){
    PyObject *m;
    m = PyModule_Create(&_bgzf_module);
    if (m == NULL){
        return NULL;
    }

    if (PyType_Ready(&VirtualFileOffset_Type) < 0) {
        return NULL;
    }
    PyObject * VirtualFileOffsetType = (PyObject *)&VirtualFileOffset_Type;
    Py_INCREF(VirtualFileOffsetType);
    if (PyModule_AddObject(m, "VirtualFileOffset", VirtualFileOffsetType) < 0){
        return NULL;
    }
    return m;
}
