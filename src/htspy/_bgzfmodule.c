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

static PyTypeObject VirtualFileOffset_Type; // Forward declaration

static PyObject * VirtualFileOffset_FromUint64(uint64_t i) {
    VirtualFileOffset * vfo = PyObject_NEW(VirtualFileOffset, 
                                           &VirtualFileOffset_Type);
    vfo->voffset = i;
    return (PyObject *)vfo;
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
    {"coffset", VirtualFileOffset_coffset_get, NULL, 
     "Offset to the beginning of a BGZF block", NULL},
    {"uoffset", VirtualFileOffset_uoffset_get, NULL, 
     "Offset inside the BGZF block.", NULL},
    {"_voffset", VirtualFileOffset__voffset_get, NULL, 
     "The internal virtual file offset integer.", NULL},
    {NULL}
}
