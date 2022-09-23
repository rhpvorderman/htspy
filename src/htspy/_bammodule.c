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
#include "structmember.h"         // PyMemberDef

#include "_conversions.h"

// Py_SET_SIZE, Py_SET_REFCNT and Py_SET_TYPE where all introduced and 
// recommended in Python 3.9
#if (PY_VERSION_HEX > 0x03090000)
    #define PyObject_SET_REFCOUNT(op, count) Py_SET_REFCNT(op, count)
    #define BamCigar_SET_SIZE(op, size) Py_SET_SIZE(op, size)
    #define BamCigar_SET_TYPE(op) Py_SET_TYPE(op, &BamCigar_Type)
#else 
    #define PyObject_SET_REFCOUNT(op, count) Py_REFCNT(op) = count
    #define BamCigar_SET_SIZE(op, size) Py_SIZE(op) = size
    #define BamCigar_SET_TYPE(op) Py_TYPE(op) = &BamCigar_Type
#endif
#define BAM_CIGAR_MAX_COUNT 0xFFFFFFF
#define BAM_CIGAR_MAX_OP 9

#define BGZF_BLOCK_SIZE 0xff00  // From bgzf.h

typedef struct {
    PyObject_VAR_HEAD
    uint32_t cigar[0];
} BamCigar;

#define BamCigar_OBJECT_SIZE sizeof(BamCigar)
#define BamCigar_GET_CIGAR(op) ((BamCigar *)op)->cigar
static PyTypeObject BamCigar_Type;  // Forward declaration

/**
 * @brief Creates a new BamCigar from the given uint32_t pointer.
 *
 * @param cigar_ptr the pointer to the raw cigar array. If NULL the BamCigar 
 *                  is not initialized.
 * @param n_cigar_op number of cigar units.
 * @return PyObject*
 */
static PyObject *
BamCigar_FromPointerAndSize(uint32_t * cigar_ptr, Py_ssize_t n_cigar_op) {
    BamCigar * obj;
    size_t size = n_cigar_op * sizeof(uint32_t);
    if (size > PY_SSIZE_T_MAX) {
        PyErr_SetString(PyExc_OverflowError, "Cigar array too large.");
        return NULL;
    }
    obj = PyObject_Malloc(BamCigar_OBJECT_SIZE + size);
    if (obj == NULL) {
        return PyErr_NoMemory();
    }
    PyObject_SET_REFCOUNT(obj, 1);
    BamCigar_SET_TYPE(obj);
    BamCigar_SET_SIZE(obj, n_cigar_op);
    if (cigar_ptr != NULL) {
        memcpy(obj->cigar, cigar_ptr, size);
    }
    return (PyObject *)obj;
}

static int
_BamCigar_Resize(PyObject **obj, Py_ssize_t new_n_cigar_op) {
    PyObject * orig_obj = *obj;
    size_t new_size = new_n_cigar_op * sizeof(uint32_t);
    if (new_n_cigar_op < 0 || Py_REFCNT(orig_obj) != 1) {
        PyErr_BadInternalCall();
        return -1;
    }
    if (Py_SIZE(orig_obj) == new_n_cigar_op) {
        return 0;
    }
    if (new_n_cigar_op == 0){
        *obj = BamCigar_FromPointerAndSize(NULL, 0);
        Py_DECREF(orig_obj);
        return 0;
    }
    *obj = PyObject_Realloc(orig_obj, BamCigar_OBJECT_SIZE + new_size);
    if (*obj == NULL) {
        PyObject_Free(orig_obj);
        PyErr_NoMemory();
        return -1;
    }
    BamCigar_SET_SIZE(*obj, new_n_cigar_op);
    return 0;
}

static PyObject *
BamCigar_number_of_operations(BamCigar *self, void *closure){
    return PyLong_FromSsize_t(Py_SIZE(self));
}

static PyGetSetDef BamCigar_properties[] = {
    {"number_of_operations", (getter)BamCigar_number_of_operations, NULL,
    "The number of CIGAR operations (n_cigar_op).", NULL},
    {NULL},
};

static inline int
BamCigar_equals(BamCigar *self, BamCigar *other) {
    if (Py_TYPE(other) != &BamCigar_Type) {
        return 0;
    }
    if (Py_SIZE(self) != Py_SIZE(other)) {
        return 0;
    }
    Py_ssize_t size = Py_SIZE(self) * sizeof(uint32_t);
    return (memcmp(self->cigar, other->cigar, size) == 0);
}

static PyObject *
BamCigar_richcompare(BamCigar *self, BamCigar *other, int op) {
    switch (op) {
        case Py_EQ:
           return PyBool_FromLong(BamCigar_equals(self, other));
        case Py_NE:
            return PyBool_FromLong(!BamCigar_equals(self, other));
        default:
            Py_RETURN_NOTIMPLEMENTED;
    }
}

static int
BamCigar_get_buffer(BamCigar *self, Py_buffer *view, int flags) {
    // Code adapted from PyBuffer_FillInfo as it is nearly the same.
    // Python bytes objects do this too.
    if ((flags & PyBUF_WRITABLE) == PyBUF_WRITABLE) {
        PyErr_SetString(PyExc_BufferError,
                        "BamCigar is not writable.");
       return -1;
    }
    Py_INCREF(self);
    view->obj = (PyObject *)self;
    view->buf = (void *)self->cigar;  // Discard const, readonly is set.
    view->len = sizeof(uint32_t) * Py_SIZE(self);
    view->readonly = 1;
    view->itemsize = sizeof(uint32_t);
    view->format = NULL;
    if ((flags & PyBUF_FORMAT) == PyBUF_FORMAT)
        view->format = "I";
    view->ndim = 1;
    view->shape = NULL;
    if ((flags & PyBUF_ND) == PyBUF_ND)
        view->shape = &(Py_SIZE(self));
    view->strides = NULL;
    if ((flags & PyBUF_STRIDES) == PyBUF_STRIDES)
        view->strides = &(view->itemsize);
    view->suboffsets = NULL;
    view->internal = NULL;
    return 0;
}

static PyBufferProcs BamCigar_as_buffer = {
    .bf_getbuffer = (getbufferproc)BamCigar_get_buffer,
    .bf_releasebuffer = NULL, // Buffer does not use resources
};

static PyObject *
BamCigar__str__(BamCigar *self) {
    // Largest cigar op length is 9 digits (268435455). So 9 digits plus 1
    // op character == 10 characters per cigar op. Assigning max_size memory
    // has the disadvantage that we probably assign way too much memory, but
    // at the advantage that sprintf can never overshoot, so there is no need
    // to check, and the memory never has to be resized.
    Py_ssize_t n_cigar_op = Py_SIZE(self);
    uint32_t * cigar = self->cigar;
    size_t max_size = n_cigar_op * 10;
    char * buffer = PyMem_Malloc(max_size);
    if (buffer == NULL) {
        return PyErr_NoMemory();
    }
    uint32_t cigar_int;
    Py_ssize_t string_size = 0;
    Py_ssize_t i = 0;
    while (i < n_cigar_op) {
        cigar_int = cigar[i];
        // No snprintf, because overshoot is impossible.
        string_size += sprintf(buffer + string_size, "%d%c",
                               bam_cigar_oplen(cigar_int),
                               bam_cigar_opchr(cigar_int));
        i += 1;
    }
    // PyUnicode_New is faster than PyUnicode_DecodeASCII, since we do not need
    // to check for ASCII. The above code can never contain non-ASCII characters.
    // This trick was learned from Marcel Martin. Thanks!
    PyObject * retval = PyUnicode_New(string_size, 127);
    if (retval == NULL) {
        PyMem_Free(buffer);
        return PyErr_NoMemory();
    }
    memcpy(PyUnicode_1BYTE_DATA(retval), buffer, string_size);
    PyMem_Free(buffer);
    return retval;
}

static PyObject *
BamCigar__repr__(BamCigar * self){
    PyObject * cigarstring = BamCigar__str__(self);
    if (cigarstring == NULL){
        return NULL;
    }
    return PyUnicode_FromFormat("Cigar(%R)", cigarstring);
}

PyDoc_STRVAR(BamCigar_from_iter__doc__,
"from_iter($cls, cigartuples, /)\n"
"--\n"
"\n"
"Create a new BamCigar from an iterable of (operation, count) tuples.\n"
);
#define BAM_CIGAR_FROM_ITER_METHODDEF    \
    {"from_iter", (PyCFunction)(void(*)(void))BamCigar_from_iter, \
    METH_O | METH_CLASS, BamCigar_from_iter__doc__}

#define BAMCIGAR_FROM_ITER_ERROR_EXIT \
    Py_DECREF(cigartuples);Py_DECREF(cigar_obj); return NULL;

static PyObject *
BamCigar_from_iter(PyTypeObject *type, PyObject *cigartuples_in) {
    PyObject * cigartuples = PySequence_Fast(
        cigartuples_in, "cigartuples must be an iterable");
    if (cigartuples == NULL){
        return NULL;
    }
    Py_ssize_t n_cigar_op = PySequence_Fast_GET_SIZE(cigartuples);
    PyObject * cigar_obj = BamCigar_FromPointerAndSize(NULL, n_cigar_op);
    if (cigar_obj == NULL){
        Py_DECREF(cigartuples);
        return PyErr_NoMemory();
    }
    uint32_t * cigar = BamCigar_GET_CIGAR(cigar_obj);
    Py_ssize_t i = 0;
    PyObject * tup;
    PyObject * operation;
    PyObject * count;
    Py_ssize_t operation_i;
    Py_ssize_t count_i;

    while (i < n_cigar_op) {
        tup = PySequence_Fast_GET_ITEM(cigartuples, i);
        if (!PyTuple_CheckExact(tup)) {
            PyErr_Format(PyExc_TypeError,
                         "List should only consist of tuples got '%s' for "
                         "item: %R",
                         Py_TYPE(tup)->tp_name, tup);
            BAMCIGAR_FROM_ITER_ERROR_EXIT
        }
        if (PyTuple_GET_SIZE(tup) != 2) {
            PyErr_Format(PyExc_ValueError,
                         "Tuples should consist of 2 items got '%ld' for "
                         "item: %r",
                          PyTuple_GET_SIZE(tup), tup);
            BAMCIGAR_FROM_ITER_ERROR_EXIT
        }
        operation = PyTuple_GET_ITEM(tup, 0);
        count = PyTuple_GET_ITEM(tup, 1);
        if (!PyLong_CheckExact(operation)) {
              PyErr_Format(PyExc_TypeError,
                           "Operation should be of type int, got '%s' for "
                           "cigartuple: %R",
                           Py_TYPE(operation)->tp_name, tup);
            BAMCIGAR_FROM_ITER_ERROR_EXIT
        }
        if (!PyLong_CheckExact(count)) {
              PyErr_Format(PyExc_TypeError,
                           "Count should be of type int, got '%s' for "
                           "cigartuple: %R",
                           Py_TYPE(count)->tp_name, tup);
            BAMCIGAR_FROM_ITER_ERROR_EXIT
        }
        operation_i = PyLong_AsSsize_t(operation);
        count_i = PyLong_AsSsize_t(count);
        if ((operation_i > BAM_CIGAR_MAX_OP) || (operation_i < 0)){
            PyErr_Format(
                PyExc_ValueError,
                "Operation should be between 0 and %d. "
                "Got %ld for cigartuple: %R",
                BAM_CIGAR_MAX_OP, operation_i, tup);
            BAMCIGAR_FROM_ITER_ERROR_EXIT
        }
        if ((count_i > BAM_CIGAR_MAX_COUNT) || (count_i < 0)) {
            PyErr_Format(
                PyExc_ValueError,
                "Count should be between 0 and %d. "
                "Got %ld for cigartuple: %R",
                BAM_CIGAR_MAX_COUNT, count_i, tup);
        }
        cigar[i] = bam_cigar_gen(count_i, operation_i);
        i += 1;
    }
    Py_DECREF(cigartuples);
    return cigar_obj;
}

PyDoc_STRVAR(BamCigar_from_buffer__doc__,
"from_buffer($cls, data, /)\n"
"--\n"
"\n"
"Create a new BamCigar from an object that supports the buffer protocol.\n"
"\n"
"Objects that support the buffer protocol include bytes, bytesarray, array.array,\n"
"numpy arrays and others."
);

#define BAM_CIGAR_FROM_BUFFER_METHODDEF    \
    {"from_buffer", (PyCFunction)(void(*)(void))BamCigar_from_buffer, \
    METH_O | METH_CLASS, BamCigar_from_buffer__doc__}

static PyObject *
BamCigar_from_buffer(PyTypeObject *type, PyObject *data) {
    Py_buffer buffer;
    if (PyObject_GetBuffer(data, &buffer, PyBUF_SIMPLE) != 0) {
        return NULL;
    }
    if (buffer.len % 4) {
        PyErr_SetString(PyExc_ValueError,
            "buffer length not a multiple of 4");
        return NULL;
    }
    Py_ssize_t n_cigar_op = buffer.len / 4;
    return BamCigar_FromPointerAndSize((uint32_t *)buffer.buf, n_cigar_op);
}

PyDoc_STRVAR(BamCigar_init__doc__,
"__init__($cls, cigarstring, /)\n"
"--\n"
"\n"
"Create a new BamCigar from a cigarstring.\n"
);

static PyObject *
BamCigar__new___(PyTypeObject *type, PyObject *args, PyObject *kwargs) {
    PyObject * cigarstring = NULL;
    char * keywords[] = {"", NULL};
    const char *format = "O|:BamCigar.__init__";
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, format, keywords, &cigarstring)) {
        return NULL;
    }
    if (!PyUnicode_CheckExact(cigarstring)) {
        PyErr_Format(PyExc_TypeError, "cigarstring must be of type str, got %s",
            Py_TYPE(cigarstring)->tp_name);
        return NULL;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(cigarstring)) {
        PyErr_SetString(PyExc_ValueError,
            "cigarstring must be a valid ascii string");
        return NULL;
    }
    Py_ssize_t string_size = PyUnicode_GET_LENGTH(cigarstring);
    char * cigar_string_ptr = (char *)PyUnicode_1BYTE_DATA(cigarstring);
    char * cigar_string_end = cigar_string_ptr + string_size;
    // uint32_t is 4 bytes. A string needs at least 2 characters to encode a
    // a cigarop + count. So maximum number of cigarops is string_size / 2.
    Py_ssize_t maximum_cigar_op = string_size / 2;
    PyObject * cigar_obj = BamCigar_FromPointerAndSize(NULL, maximum_cigar_op);
    if (cigar_obj == NULL) {
        return PyErr_NoMemory();
    }
    uint32_t * cigar = BamCigar_GET_CIGAR(cigar_obj);
    Py_ssize_t n_cigar_op = 0;
    char * endptr = NULL;
    long int count;
    char operation;
    char * cursor = cigar_string_ptr;
    while (cursor < cigar_string_end){
        count = strtol(cursor, &endptr, 10);
        if ((count < 0)) {
            PyErr_Format(PyExc_ValueError, "Invalid cigarstring: %R",
                         cigarstring);
            Py_DECREF(cigar_obj); return NULL;
        }
        if ((count > BAM_CIGAR_MAX_COUNT)) {
            PyErr_Format(
                PyExc_ValueError, "Maximum count exceeded: %ld > %ld",
                count, BAM_CIGAR_MAX_COUNT);
            Py_DECREF(cigar_obj); return NULL;
        }
        if (endptr >= cigar_string_end) {
            PyErr_Format(
                PyExc_ValueError, "Truncated cigarstring: %R",
                    cigarstring);
            Py_DECREF(cigar_obj); return NULL;
        }
        operation = bam_cigar_table[(uint8_t)endptr[0]];
        if (operation == -1) {
            PyErr_Format(PyExc_ValueError, "Invalid cigar operation: '%c'",
                endptr[0]);
            Py_DECREF(cigar_obj); return NULL;
        }
        cigar[n_cigar_op] = bam_cigar_gen(count, operation);
        n_cigar_op += 1;
        cursor = endptr + 1;
    }
    // Make sure the bytes object is made smaller if necessary.
    if (_BamCigar_Resize(&cigar_obj, n_cigar_op) == -1){
        Py_DECREF(cigar_obj); return NULL;
    }
    return cigar_obj;
}

static PyMethodDef BamCigar_methods[] = {
    BAM_CIGAR_FROM_ITER_METHODDEF,
    BAM_CIGAR_FROM_BUFFER_METHODDEF,
    {NULL}
};


typedef struct {
    PyObject_HEAD
    PyObject * bam_cigar;
    uint32_t * cigar;
    Py_ssize_t n_cigar_op;
    Py_ssize_t pos;
} BamCigarIter;

static void
BamCigarIter_dealloc(BamCigarIter * self) {
    Py_CLEAR(self->bam_cigar);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject*
BamCigarIter__next__(BamCigarIter * self) {
    if (self->pos == self->n_cigar_op) {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    uint32_t c = self->cigar[self->pos];
    Py_ssize_t cigar_op = bam_cigar_op(c);
    Py_ssize_t cigar_oplen = bam_cigar_oplen(c);
    self->pos += 1;
    PyObject * tup = PyTuple_New(2);
    PyTuple_SET_ITEM(tup, 0, PyLong_FromSsize_t(cigar_op));
    PyTuple_SET_ITEM(tup, 1, PyLong_FromSsize_t(cigar_oplen));
    return tup;
}

static PyObject *
BamCigarIter__iter__(BamCigarIter * self) {
    Py_INCREF(self);
    return (PyObject *)self;
}

static PyTypeObject BamCigarIter_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_bam.CigarIter",
    .tp_basicsize = sizeof(BamCigarIter),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_dealloc = (destructor)BamCigarIter_dealloc,
    .tp_iter = (getiterfunc)BamCigarIter__iter__,
    .tp_iternext = (iternextfunc)BamCigarIter__next__,
};

static PyObject *
BamCigar__iter__(BamCigar * self) {
    BamCigarIter *iter = PyObject_NEW(BamCigarIter, &BamCigarIter_Type);
    Py_INCREF(self);
    iter->bam_cigar = (PyObject *)self;
    iter->cigar = self->cigar;
    iter->n_cigar_op = Py_SIZE(self);
    iter->pos = 0;
    return (PyObject *)iter;
}

static PyTypeObject BamCigar_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_bam.Cigar",
    .tp_basicsize = sizeof(BamCigar),
    .tp_itemsize = sizeof(uint32_t),
    .tp_doc = BamCigar_init__doc__,
    .tp_methods = BamCigar_methods,
    .tp_getset = BamCigar_properties,
    .tp_new = BamCigar__new___,
    .tp_iter = (getiterfunc)BamCigar__iter__,
    .tp_richcompare = (richcmpfunc)BamCigar_richcompare,
    .tp_as_buffer = &BamCigar_as_buffer,
    .tp_str = (reprfunc)BamCigar__str__,
    .tp_repr = (reprfunc)BamCigar__repr__,
    .tp_free = PyObject_Del,
};

typedef struct {
    PyObject_HEAD
    uint32_t block_size;
    int32_t refID; 
    int32_t pos; 
    uint8_t l_read_name;
    uint8_t mapq;
    uint16_t bin;
    uint16_t n_cigar_op;
    uint16_t flag; 
    uint32_t l_seq;
    int32_t next_refID;
    int32_t next_pos;
    int32_t tlen;
    // The compiler automatically inserts 4 padding bytes here to properly
    // align the PyObject pointers in memory.
    PyObject * read_name;
    PyObject * bamcigar;
    PyObject * seq;
    PyObject * qual;
    PyObject * tags;
} BamRecord;

# define BAM_PROPERTIES_STRUCT_START offsetof(BamRecord, block_size)
# define BAM_PROPERTIES_STRUCT_SIZE 36

static void
BamRecord_dealloc(BamRecord *self) {
    Py_CLEAR(self->read_name);
    Py_CLEAR(self->bamcigar);
    Py_CLEAR(self->seq);
    Py_CLEAR(self->qual);
    Py_CLEAR(self->tags);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

// Convert an ASCII string to a bytes object if necessary.
static inline PyObject *
convert_to_new_bytes_reference(PyObject *obj, const char * param_name)
{
    if (obj == NULL) {
        return PyBytes_FromStringAndSize("", 0);
    }
    if (PyBytes_CheckExact(obj)) {
        Py_INCREF(obj);
        return obj;
    }
    if (PyUnicode_CheckExact(obj)) {
        return PyUnicode_AsASCIIString(obj);
    }
    PyErr_Format(PyExc_TypeError, 
                 "'%s' expected a bytes or str object got: %s",
                  param_name,
                  obj->ob_type->tp_name);
    return NULL;
}

static int
BamRecord_init(BamRecord *self, PyObject *args, PyObject *kwargs) {
    char * keywords[] = {
        "reference_id", "position", "read_name", "mapping_quality", 
        "flag", "next_reference_id, next_position", NULL};
    const char *format = "|IIObHII:BamRecord.__init__";
    int32_t reference_id = -1; 
    int32_t position = -1;
    PyObject * read_name = NULL;
    uint8_t mapping_quality = 255;
    uint16_t flag = 0;
    int32_t next_reference_id = -1;
    int32_t next_position = -1;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, format, keywords, 
        &reference_id, &position, &read_name, &mapping_quality, &flag, 
        &next_reference_id, &next_position)) {
        return -1; 
    }
    read_name = convert_to_new_bytes_reference(read_name, "read_name");
    if (read_name == NULL) {
            return -1;
    }
    self->refID = reference_id;
    self->pos = position;
    self->l_read_name = PyBytes_GET_SIZE(read_name) + 1;
    self->mapq = mapping_quality;
    self->bin = 0;
    self->n_cigar_op = 0;
    self->flag = flag;
    self->l_seq = 0;
    self->next_refID = next_reference_id;
    self->next_pos = next_position;
    self->tlen = 0;
    self->read_name = read_name;
    self->bamcigar = BamCigar_FromPointerAndSize(NULL, 0);
    self->seq = PyBytes_FromStringAndSize("", 0);
    self->qual = PyBytes_FromStringAndSize("", 0);
    self->tags = PyBytes_FromStringAndSize("", 0);
    self->block_size = (32 + self->l_read_name + (self->l_seq + 1 / 2) + 
                        self->l_seq + self->n_cigar_op * 4 + 
                        PyBytes_GET_SIZE(self->tags));
    return 0;
}

static PyMemberDef BamRecord_members[] = {
    // All the underlying BAM struct members should be accessible read only.
    // The BAM spec names are accessible READONLY by prepending an underscore. 
    // This way we communicate to the user that they are internal and readonly 
    // while still providing full access for power users. We also circumvent 
    // the dilemma BAM spec names vs Pythonic readable names.
    {"_block_size", T_UINT, offsetof(BamRecord, block_size), READONLY},
    {"_refID", T_INT, offsetof(BamRecord, refID), READONLY},
    {"_pos", T_INT, offsetof(BamRecord, pos), READONLY},
    {"_l_read_name", T_UBYTE, offsetof(BamRecord, l_read_name), READONLY},
    {"_mapq", T_UBYTE, offsetof(BamRecord, mapq), READONLY},
    {"_bin", T_USHORT, offsetof(BamRecord, bin), READONLY},
    {"_n_cigar_op", T_USHORT, offsetof(BamRecord, n_cigar_op), READONLY},
    {"_flag", T_USHORT, offsetof(BamRecord, flag), READONLY},
    {"_l_seq", T_UINT, offsetof(BamRecord, l_seq), READONLY},
    {"_next_refID", T_INT, offsetof(BamRecord, next_refID), READONLY},
    {"_next_pos", T_INT, offsetof(BamRecord, next_pos), READONLY},
    {"_tlen", T_INT, offsetof(BamRecord, tlen), READONLY},
    {"_cigar", T_OBJECT_EX, offsetof(BamRecord, bamcigar), READONLY},
    {"_seq", T_OBJECT_EX, offsetof(BamRecord, seq), READONLY},
    {"_qual", T_OBJECT_EX, offsetof(BamRecord, qual), READONLY},
    {"_tags", T_OBJECT_EX, offsetof(BamRecord, tags), READONLY},
    
    // Pythonic naming for convenient access. Everything here should be 
    // READONLY. Values that are not readonly should be set trough properties 
    // to ensure the internal consistency of the BamRecord. (Correct lengths 
    // etc.)
    {"reference_id", T_INT, offsetof(BamRecord, refID), READONLY, 
     "The index number referring to the reference."},
    {"position", T_INT, offsetof(BamRecord, pos), READONLY, 
     "The leftmost position where the template alignment starts (0-based)."},
    {"mapping_quality", T_UBYTE, offsetof(BamRecord, mapq), READONLY, 
     "The quality of the mapping."},
    {"flag", T_USHORT, offsetof(BamRecord, flag), READONLY, 
     "Bitwise flags."},
    {"next_position", T_INT, offsetof(BamRecord, next_pos), READONLY, 
     "next_pos: The leftmost position of the next segment."},
    {"template_length", T_INT, offsetof(BamRecord, tlen), READONLY},
    {"qualities", T_OBJECT_EX, offsetof(BamRecord, qual), READONLY},
    {NULL}
};

// PROPERTIES

PyDoc_STRVAR(BamRecord__read_name_doc,
"The name of the aligned read as an ASCII encoded bytes object.\n\n"
"This is the underlying property in the BamRecord object. read_name converts\n"
"this property upon access to a string.");

static PyObject * 
BamRecord_get__read_name(BamRecord * self, void* closure) {
    Py_INCREF(self->read_name);
    return self->read_name;
}

static int 
BamRecord_set__read_name(BamRecord * self, PyObject * new_read_name, void* closure) 
{
    if (!PyBytes_CheckExact(new_read_name)){
        PyErr_SetString(PyExc_TypeError, "read_name must be a bytes object");
        return -1;
    }
    Py_ssize_t read_name_size = PyBytes_GET_SIZE(new_read_name);
    if (read_name_size > 254) {
        PyErr_SetString(PyExc_ValueError, 
            "read_name may not be larger than 254 characters.");
        return -1;
    }
    PyObject * old_read_name = self->read_name;
    Py_INCREF(new_read_name);
    self->read_name = new_read_name;
    Py_DECREF(old_read_name);
    // Make sure the internal sizes are correct after updating.
    uint8_t old_l_read_name = self->l_read_name;
    self->l_read_name = (uint8_t)read_name_size + 1;
    self->block_size = self->block_size + self->l_read_name - old_l_read_name;
    return 0;
}

PyDoc_STRVAR(BamRecord_read_name_doc,
"The name of the aligned read as a string.\n");

static PyObject * 
BamRecord_get_read_name(BamRecord * self, void* closure) 
{
    return PyUnicode_DecodeASCII(PyBytes_AS_STRING(self->read_name), 
                                 PyBytes_GET_SIZE(self->read_name), 
                                 "strict");
}

static int 
BamRecord_set_read_name(BamRecord * self, PyObject * new_qname, void* closure) 
{
    PyObject * new_read_name = PyUnicode_AsASCIIString(new_qname);
    if (new_read_name == NULL)
        return -1;
    int ret = BamRecord_set__read_name(self, new_read_name, closure);
    Py_DECREF(new_read_name);
    return ret;
}

PyDoc_STRVAR(BamRecord_tags_doc,
"The raw tags as a bytes object.");

static PyObject * 
BamRecord_get_tags(BamRecord * self, void* closure) 
{
    Py_INCREF(self->tags);
    return self->tags;
}

static int 
BamRecord_set_tags(BamRecord * self, PyObject * new_tags, void* closure) 
{
    if (!PyBytes_CheckExact(new_tags)){
        PyErr_SetString(PyExc_TypeError, "tags must be a bytes object");
        return -1;
    }
    Py_ssize_t new_tags_size = PyBytes_GET_SIZE(new_tags);
    Py_ssize_t old_tags_size = PyBytes_GET_SIZE(self->tags);
    PyObject * old_tags = self->tags;
    Py_INCREF(new_tags);
    self->tags = new_tags;
    Py_DECREF(old_tags);
    // Make sure the internal sizes are correct after updating.
    self->block_size = self->block_size + new_tags_size - old_tags_size;
    return 0;
}

PyDoc_STRVAR(BamRecord_cigar_doc, 
"A BamCigar object representing the CIGAR information.");

static PyObject *
BamRecord_get_cigar(BamRecord * self, void * closure) {
    if (self->n_cigar_op == 2) {
        // Initiate CG tag check
        uint32_t * cigar = BamCigar_GET_CIGAR(self->bamcigar);
        if ((bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) && 
            (bam_cigar_oplen(cigar[0]) == self->l_seq)) {
                PyErr_SetString(PyExc_NotImplementedError, 
                    "Support for cigars longer than 65536 has not yet been implemented.");
                return NULL;
            }
    }
    Py_INCREF(self->bamcigar);
    return (PyObject *)self->bamcigar;
}

static int 
BamRecord_set_cigar(BamRecord * self, BamCigar * new_cigar, void * closure) {
    if (Py_TYPE(new_cigar) != &BamCigar_Type) {
        PyErr_Format(PyExc_TypeError, "cigar must be of BamCigar type, got %s.",
            Py_TYPE(new_cigar)->tp_name);
        return -1; 
    }
    if (Py_SIZE(new_cigar) > 65536) {
        PyErr_SetString(PyExc_NotImplementedError, 
            "Support for cigars longer than 65536 has not yet been implemented.");
        return -1;
    }
    PyObject * tmp = self->bamcigar;
    Py_INCREF(new_cigar);
    self->bamcigar = (PyObject *)new_cigar;
    self->n_cigar_op = Py_SIZE(new_cigar);
    Py_DECREF(tmp);
    return 0;
}

// Flags 
#define GET_FLAG_PROP(prop_name, FLAG) \
static PyObject * \
    prop_name(BamRecord *self, void *closure) {\
    return PyBool_FromLong(self->flag & FLAG); \
}

PyDoc_STRVAR(BamRecord_is_paired_doc,
"The read is paired in sequencing, no matter whether it is mapped in a pair");
GET_FLAG_PROP(BamRecord_is_paired, BAM_FPAIRED)

PyDoc_STRVAR(BamRecord_is_proper_pair_doc,
"The read is mapped in a proper pair");
GET_FLAG_PROP(BamRecord_is_proper_pair, BAM_FPROPER_PAIR)

PyDoc_STRVAR(BamRecord_is_unmapped_doc,
"The read itself is unmapped; conflictive with is_proper_pair.");
GET_FLAG_PROP(BamRecord_is_unmapped, BAM_FUNMAP)

PyDoc_STRVAR(BamRecord_mate_is_unmapped_doc,
"The mate is unmapped");
GET_FLAG_PROP(BamRecord_mate_is_unmapped, BAM_FMUNMAP)

PyDoc_STRVAR(BamRecord_is_reverse_doc,
"The read is mapped to the reverse strand");
GET_FLAG_PROP(BamRecord_is_reverse, BAM_FREVERSE)

PyDoc_STRVAR(BamRecord_mate_is_reverse_doc,
"The mate is mapped to the reverse strand.");
GET_FLAG_PROP(BamRecord_mate_is_reverse, BAM_FMREVERSE)

PyDoc_STRVAR(BamRecord_is_read1_doc,
"This is read1");
GET_FLAG_PROP(BamRecord_is_read1, BAM_FREAD1)

PyDoc_STRVAR(BamRecord_is_read2_doc,
"This is read2");
GET_FLAG_PROP(BamRecord_is_read2, BAM_FREAD2)

PyDoc_STRVAR(BamRecord_is_secondary_doc,
"This is not the primary alignment");
GET_FLAG_PROP(BamRecord_is_secondary, BAM_FSECONDARY)

PyDoc_STRVAR(BamRecord_is_qcfail_doc,
"QC failure for this read");
GET_FLAG_PROP(BamRecord_is_qcfail, BAM_FQCFAIL)

PyDoc_STRVAR(BamRecord_is_duplicate_doc,
"Read is an optical or PCR duplicate");
GET_FLAG_PROP(BamRecord_is_duplicate, BAM_FDUP)

PyDoc_STRVAR(BamRecord_is_supplementary_doc,
"This is a supplementary alignment");
GET_FLAG_PROP(BamRecord_is_supplementary, BAM_FSUPPLEMENTARY)

static PyGetSetDef BamRecord_properties[] = {
    {"read_name", (getter)BamRecord_get_read_name, (setter)BamRecord_set_read_name, 
     BamRecord_read_name_doc, NULL},
    {"_read_name", (getter)BamRecord_get__read_name, (setter)BamRecord_set__read_name,
     BamRecord__read_name_doc, NULL},
    {"tags", (getter)BamRecord_get_tags, (setter)BamRecord_set_tags,
     BamRecord_tags_doc, NULL},
    {"cigar", (getter)BamRecord_get_cigar, (setter)BamRecord_set_cigar,
     BamRecord_cigar_doc, NULL},
    {"is_paired", (getter)BamRecord_is_paired, NULL, 
     BamRecord_is_paired_doc, NULL},
    {"is_proper_pair", (getter)BamRecord_is_proper_pair, NULL, 
     BamRecord_is_proper_pair_doc, NULL},
    {"is_unmapped", (getter)BamRecord_is_unmapped, NULL, 
     BamRecord_is_unmapped_doc, NULL},
    {"mate_is_unmapped", (getter)BamRecord_mate_is_unmapped, NULL,
     BamRecord_mate_is_unmapped_doc, NULL},
    {"is_reverse", (getter)BamRecord_is_reverse, NULL,
     BamRecord_is_reverse_doc, NULL},
    {"mate_is_reverse", (getter)BamRecord_mate_is_reverse, NULL, 
     BamRecord_mate_is_reverse_doc, NULL},
    {"is_read1", (getter)BamRecord_is_read1, NULL,
     BamRecord_is_read1_doc, NULL},
    {"is_read2", (getter)BamRecord_is_read2, NULL, 
     BamRecord_is_read2_doc, NULL},
    {"is_secondary", (getter)BamRecord_is_secondary, NULL, 
     BamRecord_is_secondary_doc},
    {"is_qcfail", (getter)BamRecord_is_qcfail, NULL,
     BamRecord_is_qcfail_doc},
    {"is_duplicate", (getter)BamRecord_is_duplicate, NULL,
     BamRecord_is_duplicate_doc, NULL},
    {"is_supplementary", (getter)BamRecord_is_supplementary, NULL,
     BamRecord_is_supplementary_doc, NULL},
    {NULL}
};

// METHODS
PyDoc_STRVAR(BamRecord_get_sequence__doc__,
"Convert the encoded sequence to an ASCII-string");

#define BAMRECORD_GET_SEQUENCE_METHODDEF \
    {"get_sequence", (PyCFunction)(void(*)(void))BamRecord_get_sequence, \
     METH_NOARGS, BamRecord_get_sequence__doc__}

static PyObject *
BamRecord_get_sequence(BamRecord * self, PyObject * Py_UNUSED(ignore)) {
    uint32_t length = self->l_seq;
    PyObject * retval = PyUnicode_New(length, 127);
    if (retval == NULL) {
        return PyErr_NoMemory();
    }
    uint8_t * decoded_sequence = (uint8_t *)PyUnicode_DATA(retval);
    uint16_t * decoded_sequence_pairs = (uint16_t *)decoded_sequence;
    uint8_t * encoded_sequence = (uint8_t *)PyBytes_AS_STRING(self->seq);
    Py_ssize_t encoded_length = PyBytes_GET_SIZE(self->seq);
    Py_ssize_t i = 0;
    uint8_t index;
    // The memory address should be aligned at an even number. Because uint16_t
    // uses two bytes. Since python strings at UCS4 (4 bytes) this should never
    // pose a problem.
    assert(!((size_t)decoded_sequence_pairs & 1));
    while (i < encoded_length) {
        index = encoded_sequence[i];
        decoded_sequence_pairs[i] = number_to_nucleotide_pair_le[index];
        i += 1;
    }
    if (length & 1) {
        // When length is odd, say 7, then 4 encoded bytes are needed to 
        // store the sequence. This leads to 4 decoded pairs and thus 8 bytes.
        // Since we only requested 7 bytes when creating the return string, 
        // the above algorithm overshoots slightly. This is not a problem, 
        // since CPython always allocates 1 byte of extra space for the 
        // terminating NULL byte. This NULL byte will have been overwritten, 
        // so we have to put it back.
        decoded_sequence[length] = 0;
    }
    return retval;
}

PyDoc_STRVAR(BamRecord_set_sequence__doc__,
"Encode an ASCII-string containing IUPAC nucleotides and set it in the BAM record.");

#define BAMRECORD_SET_SEQUENCE_METHODDEF \
    {"set_sequence", (PyCFunction)(void(*)(void))BamRecord_set_sequence, \
     METH_FASTCALL, BamRecord_set_sequence__doc__}

static PyObject *
BamRecord_set_sequence(BamRecord *self, PyObject *const *args, Py_ssize_t nargs) {
    if (nargs != 1 && nargs != 2) {
        PyErr_Format(PyExc_TypeError,
                     "set_sequence expected at least 1 argument "
                     "and at most 2 arguments, got %ld", nargs);
        return NULL;
    }
    PyObject * qualities = NULL;
    PyObject * sequence = args[0];
    if (!PyUnicode_CheckExact(sequence)) {
        PyErr_Format(PyExc_TypeError, "sequence must be of type str, got %s",
                     Py_TYPE(sequence)->tp_name);
        return NULL;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(sequence)) {
        PyErr_SetString(PyExc_ValueError, 
                        "sequence must only contain ASCII characters");
        return NULL;
    }

    Py_ssize_t sequence_length = PyUnicode_GET_LENGTH(sequence);
    PyObject * new_qualities;
    if (nargs == 2) {
        qualities = args[1];
        if (!PyBytes_CheckExact(qualities)) {
            PyErr_Format(PyExc_TypeError, 
                         "qualities must be of type bytes, got %s",
                         Py_TYPE(qualities)->tp_name);
            return NULL;
        } 
        if (PyBytes_GET_SIZE(qualities) != sequence_length) {
            PyErr_SetString(PyExc_ValueError, 
                "sequence and qualities must have the same length");
            return NULL;
        }
        new_qualities = qualities;
        Py_INCREF(new_qualities);
    } else {
        // From the spec: "When base qualities are omitted but the sequence is
        // not, qual is filled with 0xFF bytes"
        new_qualities = PyBytes_FromStringAndSize(NULL, sequence_length);
        memset(PyBytes_AS_STRING(new_qualities), 0xFF, sequence_length);
    }
    Py_ssize_t encoded_length = (sequence_length + 1) / 2;
    PyObject * encoded_sequence = PyBytes_FromStringAndSize(NULL, encoded_length);
    if (encoded_sequence == NULL) {
        Py_DECREF(new_qualities);
        return PyErr_NoMemory();
    }
    Py_ssize_t i = 0;
    Py_ssize_t j = 0;
    uint8_t * sequence_chars = (uint8_t *)PyUnicode_DATA(sequence);
    uint8_t * encoded_sequence_chars = (uint8_t *)PyBytes_AS_STRING(encoded_sequence);
    int8_t iupac_int_first;
    int8_t iupac_int_second;
    while (i < sequence_length) {
        iupac_int_first = nucleotide_to_number[sequence_chars[i]];
        if (iupac_int_first == -1) {
            PyErr_Format(PyExc_ValueError, "Not a IUPAC character: %c", 
                         sequence_chars[i]);
            Py_DECREF(encoded_sequence); 
            Py_DECREF(new_qualities); 
            return NULL;
        }
        i += 1;
        iupac_int_second = nucleotide_to_number[sequence_chars[i]];
        if (iupac_int_second == -1) {
            PyErr_Format(PyExc_ValueError, "Not a IUPAC character: %c", 
                         sequence_chars[i]);
            Py_DECREF(encoded_sequence); 
            Py_DECREF(new_qualities);  
            return NULL;
        }
        i += 1;
        encoded_sequence_chars[j] = ((iupac_int_first << 4) | iupac_int_second);
        j += 1;
    }
    PyObject * old_sequence = self->seq;
    PyObject * old_qualities = self->qual; 
    Py_ssize_t old_l_seq = self->l_seq;
    Py_ssize_t old_encoded_length = (self->l_seq + 1) / 2;
    self->seq = encoded_sequence;
    self->qual = new_qualities;
    self->l_seq = sequence_length;
    self->block_size = self->block_size + sequence_length + encoded_length - 
        (old_l_seq + old_encoded_length);
    Py_DECREF(old_sequence);
    Py_DECREF(old_qualities);
    Py_RETURN_NONE;
}

static inline int
value_type_size(uint8_t value_type) {
    switch(value_type) {
        case 'A': 
        case 'c': 
        case 'C':
            return 1;
        case 's': 
        case 'S':
            return 2;
        case 'f': 
        case 'i': 
        case 'I':
            return 4; 
        case 'd':
            return 8;
        default:
            PyErr_Format(PyExc_ValueError, "Unkown value type: %c", value_type);
            return 0;
    }
}

static inline char *
bam_array_type_to_python_type(uint8_t array_type){
    switch (array_type) {
        case 'c': return "b";
        case 'C': return "B";
        case 's': return "h";
        case 'S': return "H";
        case 'i': return "i";
        case 'I': return "I";
        case 'f': return "f";
        case 'd': return "d";
        default:
            PyErr_Format(PyExc_ValueError, "Unknown array type: %c", array_type);
            return NULL;
    }
}

static inline uint8_t
python_array_type_to_bam_type(char *array_type){
    uint8_t tp = array_type[0];
    switch (tp) {
        case 'b': return 'c';
        case 'B': return 'C';
        case 'h': return 's';
        case 'H': return 'S';
        case 'i': return 'i';
        case 'I': return 'I';
        case 'f': return 'f';
        case 'd': return 'd';
        default:
            PyErr_Format(PyExc_ValueError, "Unknown array type: %s", array_type);
            return 0;
    }
}

static const uint8_t * 
skip_tag(const uint8_t *start, const uint8_t *end) {
    if (start >= end) {
        return end;
    }
    if ((end - start) <= 4) {  // end_ptr is NULL byte
        if ((end - start) <= 2) {
            PyErr_SetString(PyExc_ValueError, "truncated tag");
        } else {
            PyErr_Format(PyExc_ValueError, "truncated tag %c%c", 
                start[0], start[1]);
        }
        return NULL; 
    }
    uint8_t type = start[2];
    const uint8_t *value_start = start + 3;
    const uint8_t *value_end;
    int value_length;
    size_t max_length;
    switch(type) {
        case 'H':
        case 'Z':
            max_length = end - value_start;
            // strnlen could also have been appropriate but that does not allow
            // for easy error checking when no terminating NULL is found.
            uint8_t *string_stop = memchr(value_start, 0, max_length);
            if (string_stop == NULL) {
                PyErr_Format(
                    PyExc_ValueError, 
                    "truncated tag: %c%c has a string value not terminated by NULL",
                    start[0], start[1]);
                return NULL;
            }
            return string_stop + 1;
        case 'B':
            type = start[3];
            value_length = value_type_size(type);
            if (value_length == 0) {
                return NULL;
            }
            value_start = start + 8;
            if (value_start >= end) {
                PyErr_Format(PyExc_ValueError, "truncated tag %c%c", 
                    start[0], start[1]);
                    return NULL;
            }
            uint32_t array_length = *((uint32_t *)(start + 4));
            value_end = value_start + (array_length * value_length);
            if (value_end > end) {
                PyErr_Format(PyExc_ValueError, "truncated tag %c%c", 
                    start[0], start[1]);
                    return NULL;          
            }
            return value_end;
        default:
            value_length = value_type_size(type);
            if (value_length == 0) {
                return NULL;
            }
            value_end = value_start + value_length;
            if (value_end > end) {
                PyErr_Format(PyExc_ValueError, "truncated tag %c%c", 
                    start[0], start[1]);
                    return NULL;          
            }
            return value_end;
    }
}

static int
find_tag(const uint8_t *tags, size_t tags_length, const uint8_t *tag, 
         const uint8_t **found_tag)
{
    uint8_t tag_left_char = tag[0];
    uint8_t tag_right_char = tag[1];
    const uint8_t *end_ptr = tags + tags_length;
    const uint8_t *tag_ptr = tags;
    const uint8_t *next_tag_ptr;
    while (tag_ptr < end_ptr) {
        if (tag_ptr + 2 >= (end_ptr)) {
            PyErr_SetString(PyExc_ValueError, "truncated tag");
            return -1;
        }
        if ((tag_ptr[0] == tag_left_char) && (tag_ptr[1]) == tag_right_char) {
            *found_tag = tag_ptr;
            return 0;
        }
        next_tag_ptr = skip_tag(tag_ptr, end_ptr);
        if (next_tag_ptr == NULL) {
            return -1;
        }
        tag_ptr = next_tag_ptr;
    }
    *found_tag = NULL;
    return 0;
}

static PyObject * 
tag_ptr_to_pyobject(const uint8_t *start, const uint8_t *end, PyObject *tag_object){
    if (start >= end) {
        PyErr_SetString(PyExc_RuntimeError, 
                        "tag_ptr_to_pyobject called with end before start.");
        return NULL;
    }
    if ((end - start) < 4) {
        if ((end -start) < 2) {
            PyErr_SetString(PyExc_ValueError, "truncated tag");
        } else {
            PyErr_Format(PyExc_ValueError, "truncated tag %c%c", 
                start[0], start[1]);
        }
        return NULL; 
    }
    uint8_t type = start[2];
    const uint8_t *value_start = start + 3;
    uint8_t *value_end;
    Py_ssize_t max_length = end - value_start;
    switch(type) {
        case 'A':
            // Only ASCII chars are allowed.
            return PyUnicode_DecodeASCII((char *)value_start, 1, NULL);
        case 'Z':
            value_end = memchr(value_start, 0, max_length);
            if (value_end == NULL) {
                break;
            }
            return PyUnicode_DecodeASCII((char *)value_start, 
                                         value_end - value_start, NULL);
        case 'c':
            return PyLong_FromLong(*(int8_t *)value_start);
        case 'C':
            return PyLong_FromLong(*(uint8_t *)value_start);
        case 's':
            if (max_length < 2) break;
            return PyLong_FromLong(*(int16_t *)value_start);
        case 'S':
            if (max_length < 2) break;
            return PyLong_FromLong(*(uint16_t *)value_start);
        case 'i':
            if (max_length < 4) break;
            return PyLong_FromLong(*(int32_t *)value_start);
        case 'I':
            if (max_length < 4) break;
            return PyLong_FromLong(*(uint32_t *)value_start);
        case 'f':
            if (max_length < 4) break;
            return PyFloat_FromDouble(*(float *)value_start);
        case 'B':
            // Export a memoryview, as this is the closest thing to a builtin 
            // array in python. Iterating over it works.
            if (max_length < 5) 
                break;
            uint8_t array_type = value_start[0];
            uint32_t array_count = *(uint32_t *)(value_start + 1);
            const uint8_t *array_start = value_start + 5;
            size_t array_max_length = end - array_start;
            int itemsize = value_type_size(array_type);
            if (!itemsize) 
                return NULL;
            char * python_array_type = bam_array_type_to_python_type(array_type);
            if (!python_array_type) 
                return NULL;
            size_t array_size = itemsize * array_count;
            if (array_size > array_max_length) break;
            Py_INCREF(tag_object);
            Py_buffer array_view = {
                .buf = (void *)(value_start + 5),  // Discard const, readonly is set.
                .obj = tag_object,
                .len = array_size,
                .itemsize = itemsize,
                .readonly = 1,
                .format = python_array_type,
                .ndim = 1,
                // Fields below are automatically set by PyMemoryView_FromBuffer
                .shape = NULL,
                .strides = NULL,
                .suboffsets = NULL,
                .internal = NULL,
            };
            PyObject *memview = PyMemoryView_FromBuffer(&array_view);
            // PyMemoryView_FromBuffer sets this to NULL, but the underlying 
            // object should be properly decreased in reference count.
            ((PyMemoryViewObject *)memview)->mbuf->master.obj = tag_object;
            return memview;
        case 'H':
            PyErr_SetString(PyExc_NotImplementedError,
                            "Decoding 'H' type tags is not yet supported.");
            return NULL;
        default:
            PyErr_Format(PyExc_ValueError, 
                         "Unknown tag type: %c for tag %c%c", 
                         type, start[0], start[1]);
            return NULL;
    }
    PyErr_Format(PyExc_ValueError, "truncated tag %c%c", 
                 start[0], start[1]);
    return NULL; 
}

PyDoc_STRVAR(BamRecord_get_tag__doc__,
"get_tag($self, tag, /)\n"
"--\n"
"\n"
"Returns the value of a tag. Raises a LookupError if not found.\n"
"\n"
"  tag\n"
"    A two-letter ASCII string.\n"
"\n");

#define BAMRECORD_GET_TAG_METHODDEF    \
    {"get_tag", (PyCFunction)(void(*)(void))BamRecord_get_tag, METH_O, \
     BamRecord_get_tag__doc__}

static PyObject * 
BamRecord_get_tag(BamRecord *self, PyObject *tag) {
    if (!PyUnicode_CheckExact(tag)) {
        PyErr_Format(PyExc_TypeError, "tag must be of type str, got %s.", 
                     Py_TYPE(tag)->tp_name);
        return NULL;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(tag)) {
        PyErr_SetString(PyExc_ValueError, "tag contains non-ASCII characters");
        return NULL;
    }
    if (!(PyUnicode_GET_LENGTH(tag) == 2)) {
        PyErr_Format(PyExc_ValueError, "tag must have length 2, got %ld", 
                                        PyUnicode_GET_LENGTH(tag));
        return NULL;
    }
    uint8_t *search_tag = (uint8_t *)PyUnicode_DATA(tag);
    uint8_t *tags = (uint8_t *)PyBytes_AS_STRING(self->tags);
    Py_ssize_t tags_length = PyBytes_GET_SIZE(self->tags);
    const uint8_t *found_tag = NULL;
    int ret = find_tag(tags, tags_length, search_tag, &found_tag);
    if (ret != 0) {
        return NULL;
    }
    if (found_tag == NULL) {
        PyErr_Format(PyExc_LookupError, "Tag not present: %S", tag);
        return NULL;
    }
    return tag_ptr_to_pyobject(found_tag, tags + tags_length, self->tags);
}

static const char *tag_to_value_type(const uint8_t *tag) {
    // Need to use a macro here for use in case statements.
    #define TAG_KEY(left, right) ((uint16_t)left << 8 | (uint16_t)(right))
    
    // Only uppercase tags are in the SAMtags spec, so alternatively a 26 by 36
    // lookup table can be used, but the switch statement is easier to use and maintain.
    // Let the compiler optimize for now.
    uint16_t key = TAG_KEY(tag[0], tag[1]);
    switch (key) {
        case TAG_KEY('T', 'S'):
            return "A";
        case TAG_KEY('A', 'M'):
        case TAG_KEY('A', 'S'):
        case TAG_KEY('C', 'M'):
        case TAG_KEY('C', 'P'):
        case TAG_KEY('F', 'I'):
        case TAG_KEY('H', '0'):
        case TAG_KEY('H', '1'):
        case TAG_KEY('H', '2'):
        case TAG_KEY('H', 'I'):
        case TAG_KEY('I', 'H'):
        case TAG_KEY('M', 'Q'):
        case TAG_KEY('N', 'H'):
        case TAG_KEY('N', 'M'):
        case TAG_KEY('O', 'P'):
        case TAG_KEY('P', 'Q'):
        case TAG_KEY('S', 'M'):
        case TAG_KEY('T', 'C'):
        case TAG_KEY('U', 'Q'):
            return "i";
        case TAG_KEY('B', 'C'):
        case TAG_KEY('B', 'Q'):
        case TAG_KEY('B', 'Z'):
        case TAG_KEY('C', 'B'):
        case TAG_KEY('C', 'C'):
        case TAG_KEY('C', 'O'):
        case TAG_KEY('C', 'Q'):
        case TAG_KEY('C', 'R'):
        case TAG_KEY('C', 'S'):
        case TAG_KEY('C', 'T'):
        case TAG_KEY('C', 'Y'):
        case TAG_KEY('E', '2'):
        case TAG_KEY('F', 'S'):
        case TAG_KEY('L', 'B'):
        case TAG_KEY('M', 'C'):
        case TAG_KEY('M', 'D'):
        case TAG_KEY('M', 'I'):
        case TAG_KEY('M', 'M'):
        case TAG_KEY('O', 'A'):
        case TAG_KEY('O', 'C'):
        case TAG_KEY('O', 'Q'):
        case TAG_KEY('O', 'X'):
        case TAG_KEY('P', 'G'):
        case TAG_KEY('P', 'T'):
        case TAG_KEY('P', 'U'):
        case TAG_KEY('Q', '2'):
        case TAG_KEY('Q', 'T'):
        case TAG_KEY('Q', 'X'):
        case TAG_KEY('R', '2'):
        case TAG_KEY('R', 'G'):
        case TAG_KEY('R', 'X'):
        case TAG_KEY('S', 'A'):
        case TAG_KEY('U', '2'):
            return "Z";
        case TAG_KEY('M', 'L'):
            return "BC";
        case TAG_KEY('F', 'Z'):
            return "BS";
        case TAG_KEY('C', 'G'):
            return "BI";
        default:
            return NULL;
    }
}

static const char *PyObject_to_value_type(PyObject *value) {
    if (PyUnicode_CheckExact(value)) {
        return "Z";
    }
    if (PyLong_CheckExact(value)) {
        return "I";
    }
    if (PyFloat_CheckExact(value)) {
        return "f";
    }
    if (PyObject_CheckBuffer(value)) {
        return "B";
    }
    PyErr_Format(PyExc_ValueError, 
                 "Could not determine appropriate tag type for %R", value);
    return NULL;
}

/**
 * @brief Add, replace or delete a tag on a BamRecord
 * 
 * This function checks if the tag already exists, if so it replaces it.
 * If the tag does not exist, it adds the new tag. To delete a tag, use
 * tag_marker_length=0 and tag_value_length=0. 
 *  
 * @param self BamRecord object
 * @param tag A pointer to a two-character tag (uint8_t[2]) or equivalent string.
 * @param tag_marker The marker for the replacing tag. Pointer to array with 
 *                   tag, value_type and optional array type + count.
 * @param tag_marker_length 3 for singular values 8 for array types. 
 * @param tag_value A pointer to the value of the new tag
 * @param tag_value_length The length of the value.
 * @return int 0 on success -1 on failure.
 */
static int _BamRecord_replace_tag(BamRecord *self, 
                                  const uint8_t *tag,
                                  const uint8_t *tag_marker,
                                  size_t tag_marker_length,
                                  uint8_t *tag_value,
                                  size_t tag_value_length)
{
    const uint8_t *tags = (uint8_t *)PyBytes_AS_STRING(self->tags);
    Py_ssize_t tags_length = PyBytes_GET_SIZE(self->tags);
    const uint8_t *tags_end = tags + tags_length;
    const uint8_t *old_tag_start = NULL;
    if (find_tag(tags, tags_length, tag, &old_tag_start) != 0) {
        return -1;
    }
    size_t before_tag_length = tags_length;
    const uint8_t *after_tag_start = NULL;
    size_t after_tag_length = 0;
    if (old_tag_start != NULL) {
        before_tag_length = old_tag_start - tags; 
        after_tag_start = skip_tag(old_tag_start, tags_end);
        if (after_tag_start == NULL) {
            return -1;
        }
        after_tag_length = tags_end - after_tag_start;
    }
    Py_ssize_t new_size = before_tag_length + 
                          after_tag_length + 
                          tag_marker_length + 
                          tag_value_length;
    size_t old_block_size = self->block_size;
    size_t new_block_size = old_block_size + new_size - tags_length;
    if (new_block_size > UINT32_MAX) {
        PyErr_SetString(PyExc_OverflowError, "Value too big to store in BamRecord");
        return -1;
    }
    PyObject *tmp = PyBytes_FromStringAndSize(NULL, new_size);
    if (tmp == NULL) {
        PyErr_NoMemory();
        return -1;
    }
    void *new_tags = PyBytes_AsString(tmp);
    memcpy(new_tags, tags, before_tag_length);
    new_tags += before_tag_length;
    memcpy(new_tags, after_tag_start, after_tag_length);
    new_tags += after_tag_length;
    memcpy(new_tags, tag_marker, tag_marker_length);
    new_tags += tag_marker_length;
    memcpy(new_tags, tag_value, tag_value_length);
    PyObject *old_tag_object = self->tags;
    self->tags = tmp;
    self->block_size = new_block_size;
    Py_DECREF(old_tag_object);
    return 0;
}

static int _BamRecord_set_array_tag(BamRecord *self,
                                    const uint8_t *tag,
                                    uint8_t array_type, 
                                    PyObject *value)
{
    Py_buffer buffer;
    uint8_t tp = array_type;
    if (tp) {
        if (PyObject_GetBuffer(value, &buffer, PyBUF_SIMPLE) != 0) {
            return -1;
        }
    }
    else {
        if (PyObject_GetBuffer(value, &buffer, PyBUF_FORMAT) != 0) {
            return -1;
        }
        tp = python_array_type_to_bam_type(buffer.format);
    }
    int itemsize = value_type_size(tp);
    if (buffer.len % itemsize) {
        PyErr_Format(
            PyExc_ValueError, 
            "Cannot set tag '%c%c' with type 'B%c' using %R. " 
            "Buffer size not a multiple of %d.",
            tag[0], tag[1], tp, value, itemsize
        );
        goto error;
    }
    size_t array_size = buffer.len / itemsize;
    if (array_size > UINT32_MAX) {
        PyErr_Format(
            PyExc_OverflowError,
            "Array size of %ld, is larger than %d", array_size, UINT32_MAX
        );
        goto error;
    }
    uint8_t tag_marker[8] = {tag[0], tag[1], 'B', tp, 0, 0, 0, 0};
    ((uint32_t *)tag_marker)[1] = (uint32_t)array_size;
    return _BamRecord_replace_tag(self, tag, tag_marker, 8, buffer.buf, buffer.len);
error:
    PyBuffer_Release(&buffer);
    return -1;
}

static int _BamRecord_set_tag(BamRecord *self, 
                              const uint8_t *tag, 
                              const uint8_t *value_type, 
                              PyObject *value) 
{
    uint8_t tag_value_store[8];
    void *tag_value = tag_value_store;
    size_t tag_value_size = 0;
    double dbl;

    switch(value_type[0]) {
        default:
            PyErr_Format(PyExc_ValueError, "Unkown format: %c", value_type[0]);
            return -1;
        case 'B':
            // Array tags are quite different and thus require quite a different
            // approach than single value tags and strings.
            return _BamRecord_set_array_tag(self, tag, value_type[1], value);
        case 'A': 
        case 'Z':
            if (!PyUnicode_CheckExact(value)) {
                PyErr_Format(
                    PyExc_TypeError, 
                    "Tag '%c%c' with value_type '%c' only accepts type str "
                    "inputs, got %s", 
                    tag[0], tag[1], value_type[0], Py_TYPE(value)->tp_name);
                return -1;
            }
            if (!PyUnicode_IS_COMPACT_ASCII(value)) {
                PyErr_Format(
                    PyExc_ValueError, 
                    "Tag '%c%c' with value_type '%c' only accepts valid ASCII "
                    "characters.",
                    tag[0], tag[1], value_type[0]);
                return -1;
            }
            tag_value = PyUnicode_DATA(value);
            if (value_type[0] == 'A') {
                if (PyUnicode_GET_LENGTH(value) != 1) {
                    PyErr_Format(
                        PyExc_ValueError, 
                        "Tag '%c%c' with value_type 'A', needs a string with "
                        "exactly one character, got %d",
                        tag[0], tag[1], PyUnicode_GET_LENGTH(value));
                    return -1;
                }
                tag_value_size = 1;
            } else {
                // + 1 to the length, as CPython ensures a NULL byte at the end of 
            // + 1 to the length, as CPython ensures a NULL byte at the end of 
                // + 1 to the length, as CPython ensures a NULL byte at the end of 
                // the string.
                tag_value_size = PyUnicode_GET_LENGTH(value) + 1;
            }
            break;
        case 'c':
            tag_value_size = StorePyObjectValue_c(value, tag_value_store, tag);
            break;
        case 'C':
            tag_value_size = StorePyObjectValue_C(value, tag_value_store, tag);
            break;
        case 's':
            tag_value_size = StorePyObjectValue_s(value, tag_value_store, tag);
            break;
        case 'S':
            tag_value_size = StorePyObjectValue_S(value, tag_value_store, tag);
            break;
        case 'i':
            tag_value_size = StorePyObjectValue_i(value, tag_value_store, tag);
            break;
        case 'I':
            tag_value_size = StorePyObjectValue_I(value, tag_value_store, tag);
            break;
        case 'f':
            dbl = PyFloat_AsDouble(value);
            if ((dbl == -1.0L) && PyErr_Occurred()) {
                return -1;
            }
            ((float *)tag_value_store)[0] = (float)dbl;
            tag_value_size = 4;
            break;
    }
    if (tag_value_size == 0) {
        return -1;
    }
    uint8_t tag_marker[3] = {tag[0], tag[1], value_type[0]};
    return _BamRecord_replace_tag(self, tag, tag_marker, 3, tag_value, tag_value_size);
}

PyDoc_STRVAR(BamRecord_set_tag__doc__,
"set_tag($self, tag, value, /, value_type=None)\n"
"--\n"
"\n"
"Set the value of a tag.\n"
"\n"
"  tag\n"
"    A two-letter ASCII string.\n"
"  value\n"
"    The value to store in the tag.\n"
"  value_type\n"
"    The value type of the tag. \n"
"    By default this automatically determined by the tag if tag is in the \n" 
"    SAMtags specification or else the value.\n"
"\n");

#define BAMRECORD_SET_TAG_METHODDEF    \
    {"set_tag", (PyCFunction)(void(*)(void))BamRecord_set_tag, \
     METH_VARARGS | METH_KEYWORDS, BamRecord_set_tag__doc__}

static PyObject *BamRecord_set_tag(BamRecord *self, PyObject *args, PyObject *kwargs)
{
    PyObject *tag_obj = NULL;
    uint8_t *tag = NULL; 
    PyObject *value = NULL;
    PyObject *value_type_obj = NULL;
    const uint8_t *value_type = NULL;
    static char *format = "O!O|O!:BamRecord.set_tag()";
    static char *keywords[] = {"", "", "value_type", NULL};
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, format, keywords, 
            &PyUnicode_Type, &tag_obj, 
            &value, 
            &PyUnicode_Type, &value_type_obj)) { 
                return NULL;
    }
    if (!PyUnicode_IS_COMPACT_ASCII(tag_obj)) {
        PyErr_SetString(PyExc_ValueError, 
                        "tag should only consist of ASCII characters");
        return NULL;
    }
    if (PyUnicode_GET_LENGTH(tag_obj) != 2) {
        PyErr_Format(
            PyExc_ValueError, 
            "tag must have a length of exactly 2, got %d",
            PyUnicode_GET_LENGTH(tag_obj)
        );
        return NULL;
    }
    tag = PyUnicode_DATA(tag_obj);
    if (value_type_obj) {
        if (!PyUnicode_IS_COMPACT_ASCII(value_type_obj)) {
            PyErr_SetString(
                PyExc_ValueError, 
                "value_type should only consist of ASCII characters");
            return NULL;
        }
        Py_ssize_t value_type_length = PyUnicode_GET_LENGTH(value_type_obj);
        if (value_type_length != 1 && value_type_length != 2) {
            PyErr_Format(
                PyExc_ValueError, 
                "value_type must have a length of 1 or 2, got %d",
                value_type_length);
            return NULL;
        }
        value_type = PyUnicode_DATA(value_type_obj);
    }
    else {
        value_type = (uint8_t *)tag_to_value_type(tag);
        if (value_type == NULL) {
            value_type = (uint8_t *)PyObject_to_value_type(value);
        }
        if (value_type == NULL) {
            return NULL;
        }
    }
    if (_BamRecord_set_tag(self, tag, value_type, value) != 0) {
        return NULL;
    }
    Py_RETURN_NONE;
}


static void
BamRecord_to_ptr(BamRecord *self, char * dest) {
    memcpy(dest, (char *)self + BAM_PROPERTIES_STRUCT_START,
         BAM_PROPERTIES_STRUCT_SIZE);
    Py_ssize_t cursor = BAM_PROPERTIES_STRUCT_SIZE;
    
    Py_ssize_t read_name_size = PyBytes_GET_SIZE(self->read_name);
    memcpy(dest + cursor, PyBytes_AS_STRING(self->read_name), read_name_size);
    cursor += read_name_size;

    // Terminate read_name with NULL byte
    dest[cursor] = 0; cursor += 1;

    Py_ssize_t cigar_char_size = Py_SIZE(self->bamcigar) * sizeof(uint32_t);
    memcpy(dest + cursor, BamCigar_GET_CIGAR(self->bamcigar), cigar_char_size);
    cursor += cigar_char_size;

    Py_ssize_t seq_size = PyBytes_GET_SIZE(self->seq);
    memcpy(dest + cursor, PyBytes_AS_STRING(self->seq), seq_size);
    cursor += seq_size;

    Py_ssize_t qual_size = PyBytes_GET_SIZE(self->qual);
    memcpy(dest + cursor, PyBytes_AS_STRING(self->qual), qual_size);
    cursor += qual_size;

    Py_ssize_t tag_size = PyBytes_GET_SIZE(self->tags);
    memcpy(dest + cursor, PyBytes_AS_STRING(self->tags), tag_size);
}


PyDoc_STRVAR(BamRecord_to_bytes__doc__,
"Return the BAM record as a bytesobject that can be written into a file.");

#define BAMRECORD_TO_BYTES_METHODDEF    \
    {"to_bytes", (PyCFunction)(void(*)(void))BamRecord_to_bytes, METH_NOARGS, \
     BamRecord_to_bytes__doc__}

static PyObject *
BamRecord_to_bytes(BamRecord *self, PyObject *NoArgs)
{
    PyObject * ret_val = PyBytes_FromStringAndSize(
        NULL, self->block_size + sizeof(self->block_size));
    if (ret_val == NULL){
        return PyErr_NoMemory();
    }
    char * bam_bytes = PyBytes_AS_STRING(ret_val);
    BamRecord_to_ptr(self, bam_bytes);
    return ret_val;
}


static PyMethodDef BamRecord_methods[] = {
    BAMRECORD_TO_BYTES_METHODDEF,
    BAMRECORD_GET_SEQUENCE_METHODDEF,
    BAMRECORD_SET_SEQUENCE_METHODDEF,
    BAMRECORD_GET_TAG_METHODDEF,
    BAMRECORD_SET_TAG_METHODDEF,
    {NULL}
};

PyDoc_STRVAR(BamRecord__doc__, 
"An object that represents an alignment record from a BAM file.");

static PyTypeObject BamRecord_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_bam.BamRecord",
    .tp_basicsize = sizeof(BamRecord),
    .tp_dealloc = (destructor)BamRecord_dealloc,      
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = BamRecord__doc__,
    .tp_methods = BamRecord_methods,
    .tp_members = BamRecord_members,   
    .tp_getset = BamRecord_properties,
    .tp_init = (initproc)BamRecord_init,
    .tp_new = PyType_GenericNew,
};

PyDoc_STRVAR(BamBlockBuffer__doc__,
"A structure to create a BGZF block from BamRecord objects.\n");

typedef struct {
    PyObject_HEAD
    Py_ssize_t buffersize;
    Py_ssize_t pos;
    char * buffer;
} BamBlockBuffer;

static void
BamBlockBuffer_dealloc(BamBlockBuffer *self) {
    PyMem_Free(self->buffer);
    Py_TYPE(self)->tp_free(self);
}

static int
BamBlockBuffer__init__(BamBlockBuffer * self, PyObject *args, PyObject *kwargs) {
    Py_ssize_t buffersize = BGZF_BLOCK_SIZE;
    self->buffer = NULL;
    char * tmp;
    char * keywords[] = {"", NULL};
    if (!PyArg_ParseTupleAndKeywords(
        args, kwargs, "|n:BamBlockBuffer", keywords, &buffersize)) {
            return -1;
    }
    if (buffersize < 0) {
        PyErr_Format(PyExc_ValueError,
                     "buffer size must be larger than 0. Got: %ld", buffersize);
        return -1;
    }
    tmp = (char *)PyMem_Malloc(buffersize);
    if (tmp == NULL) {
        PyErr_NoMemory();
        return -1;
    }
    self->buffer = tmp;
    self->buffersize = buffersize;
    self->pos = 0;
    return 0;
}

static PyMemberDef BamBlockBuffer_members[] = {
    {"buffersize", T_PYSSIZET, offsetof(BamBlockBuffer, buffersize), READONLY,
     "The size of the internal buffer."},
    {"bytes_written", T_PYSSIZET, offsetof(BamBlockBuffer, pos), READONLY,
     "The number of bytes written in the internal buffer."},
    {NULL}
};


PyDoc_STRVAR(BamBlockBuffer_write_doc,
"Write a BamRecord object into the BamBlockBuffer.\n"
"\n"
"Returns the amount of bytes written. Returns 0 if the BamRecord does not\n"
"fit in the buffer anymore");

#define BAMBLOCKBUFFER_WRITE_METHODDEF    \
    {"write", (PyCFunction)(void(*)(void))BamBlockBuffer_write, METH_O, \
     BamBlockBuffer_write_doc}

static PyObject * 
BamBlockBuffer_write(BamBlockBuffer * self, BamRecord * bam_record) {
    if (Py_TYPE(bam_record) != &BamRecord_Type) {
        PyErr_Format(PyExc_TypeError, "Type must be BamRecord, got: %s",
                     Py_TYPE(bam_record)->tp_name);
        return NULL;
    }
    Py_ssize_t record_size = bam_record->block_size + sizeof(bam_record->block_size);
    Py_ssize_t final_pos = self->pos + record_size;
    if (final_pos > self->buffersize) {
        return PyLong_FromSsize_t(0);
    }
    BamRecord_to_ptr(bam_record, self->buffer + self->pos);
    self->pos = final_pos;
    return PyLong_FromSsize_t(record_size);
}
PyDoc_STRVAR(BamBlockBuffer_reset_doc,
"Remove all records from the buffer.");

#define BAMBLOCKBUFFER_RESET_METHODDEF    \
    {"reset", (PyCFunction)(void(*)(void))BamBlockBuffer_reset, METH_NOARGS, \
     BamBlockBuffer_reset_doc}

static PyObject *
BamBlockBuffer_reset(BamBlockBuffer *self, PyObject *Py_UNUSED(ignore)) {
    self->pos = 0;
    Py_RETURN_NONE;
}

PyDoc_STRVAR(BamBlockBuffer_get_block_view_doc,
"Return a memoryview over all the memory blocks written so far.");

#define BAMBLOCKBUFFER_GET_BLOCK_VIEW_METHODDEF    \
    {"get_block_view", (PyCFunction)(void(*)(void))BamBlockBuffer_get_block_view, \
     METH_NOARGS, BamBlockBuffer_get_block_view_doc}

static PyObject *
BamBlockBuffer_get_block_view(BamBlockBuffer *self, PyObject *Py_UNUSED(ignore)) {
    return PyMemoryView_FromMemory(self->buffer, self->pos, PyBUF_READ);
}

static PyMethodDef BamBlockBuffer_methods[] = {
    BAMBLOCKBUFFER_WRITE_METHODDEF,
    BAMBLOCKBUFFER_RESET_METHODDEF,
    BAMBLOCKBUFFER_GET_BLOCK_VIEW_METHODDEF,
    {NULL},
};

static PyTypeObject BamBlockBuffer_Type = {
    .tp_name = "_bam.BamBlockBuffer",
    .tp_basicsize = sizeof(BamBlockBuffer),
    .tp_dealloc = (destructor)BamBlockBuffer_dealloc,
    .tp_init = (initproc)BamBlockBuffer__init__,
    .tp_new = PyType_GenericNew,
    .tp_doc = BamBlockBuffer__doc__,
    .tp_members = BamBlockBuffer_members,
    .tp_methods = BamBlockBuffer_methods,
};

typedef struct {
    PyObject_HEAD 
    Py_buffer view; 
    char * buf;
    Py_ssize_t pos;
    Py_ssize_t len; 
} BamIterator;

static void
BamIterator_dealloc(BamIterator *self) {
    PyBuffer_Release(&(self->view));
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static BamIterator *
BamIterator_iter(BamIterator *self){
    Py_INCREF(self);
    return self;
}

static PyObject *
BamIterator_iternext(BamIterator *self){
    if (self->pos >= self->len){
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    Py_ssize_t start_pos = self->pos;

    BamRecord * bam_record = PyObject_New(BamRecord, &BamRecord_Type);
    bam_record->read_name = NULL;
    bam_record->seq = NULL;
    bam_record->bamcigar = NULL;
    bam_record->qual = NULL;
    bam_record->tags = NULL;

    if ((self->len - self->pos) < BAM_PROPERTIES_STRUCT_SIZE) {
        PyErr_SetString(PyExc_EOFError, "Truncated BAM record");
        Py_DECREF(bam_record);
        return NULL;
    }
    // Copy the bam file data directly into the struct.
    memcpy((char *)bam_record + BAM_PROPERTIES_STRUCT_START,
            self->buf + self->pos,
            BAM_PROPERTIES_STRUCT_SIZE);

    // Block_size is excluding the block_size field it self.
    Py_ssize_t record_length = bam_record->block_size + 4;
    if (self->pos + record_length > self->len) {
        PyErr_SetString(PyExc_EOFError, "Truncated BAM record");
        Py_DECREF(bam_record);
        return NULL;
    }
    self->pos += BAM_PROPERTIES_STRUCT_SIZE;
    bam_record->read_name = PyBytes_FromStringAndSize(
        self->buf + self->pos, bam_record->l_read_name -1);
    self->pos += bam_record->l_read_name;

    Py_ssize_t cigar_length = bam_record->n_cigar_op * sizeof(uint32_t);
    bam_record->bamcigar = BamCigar_FromPointerAndSize(
        (uint32_t *)(self->buf + self->pos), bam_record->n_cigar_op);
    self->pos += cigar_length;

    Py_ssize_t seq_length = (bam_record->l_seq + 1) / 2;
    bam_record->seq = PyBytes_FromStringAndSize(
        self->buf + self->pos, seq_length);
    self->pos += seq_length;

    bam_record->qual = PyBytes_FromStringAndSize(
        self->buf + self->pos, bam_record->l_seq);
    self->pos += bam_record->l_seq;

    // Tags are in the remaining block of data.
    Py_ssize_t tags_length = start_pos + record_length - self->pos;
    bam_record->tags = PyBytes_FromStringAndSize(
        self->buf + self-> pos, tags_length);
    
    // Should be equal to record_length
    self->pos += tags_length;

    // Check if any of the bytes objects was NULL. This means there was 
    // no memory available.
    if ((bam_record->read_name == NULL) | (bam_record->tags == NULL) |
        (bam_record->seq == NULL) | (bam_record->qual == NULL) |
        (bam_record->tags == NULL)) {
        Py_DECREF(bam_record);
        return PyErr_NoMemory();
    }

    return (PyObject *)bam_record;
}

static PyTypeObject BamIterator_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_bam.BamIterator", 
    .tp_basicsize = sizeof(BamIterator),
    .tp_dealloc =(destructor)BamIterator_dealloc,  
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_iter = (getiterfunc)BamIterator_iter,
    .tp_iternext = (iternextfunc)BamIterator_iternext
};

PyDoc_STRVAR(bam_iterator_doc,
"bam_iterator($module, data, /)\n"
"--\n"
"\n"
"Return an iterator that yields BamRecord objects.\n"
"\n"
"  data\n"
"    A block of raw BAM Record data. May be any object\n"
"    that supports the buffer protocol: bytes, bytearray, memoryview.\n"
);
static PyObject * 
bam_iterator(PyObject *module, PyObject * obj) {
    BamIterator *self = PyObject_New(BamIterator, &BamIterator_Type);
    if (!PyObject_GetBuffer(obj, &(self->view), PyBUF_SIMPLE) == 0) {
        Py_DECREF(self);
        return NULL;
    }
    self->buf = self->view.buf;
    self->pos = 0;
    self->len = self->view.len;
    return (PyObject *)self;
}

static PyMethodDef _bam_methods[] = {
    {"bam_iterator", (PyCFunction)(void(*)(void))bam_iterator,
     METH_O, bam_iterator_doc},
    {NULL}
};

static struct PyModuleDef _bam_module = {
    PyModuleDef_HEAD_INIT,
    "_bam",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    _bam_methods  /* module methods */
};


PyMODINIT_FUNC
PyInit__bam(void)
{
    PyObject *m;

    m = PyModule_Create(&_bam_module);
    if (m == NULL)
        return NULL;

    if (PyType_Ready(&BamIterator_Type) < 0)
        return NULL;
    PyObject * BamiteratorType = (PyObject *)&BamIterator_Type;
    Py_INCREF(BamiteratorType);
    if (PyModule_AddObject(m, "BamIterator", BamiteratorType) < 0)
        return NULL;

    if (PyType_Ready(&BamRecord_Type) < 0)
        return NULL;
    PyObject * BamRecordType = (PyObject *)&BamRecord_Type;
    Py_INCREF(BamRecordType);
    if (PyModule_AddObject(m, "BamRecord", BamRecordType) < 0) {
        return NULL;
    }

    if (PyType_Ready(&BamBlockBuffer_Type) < 0)
        return NULL;
    PyObject * BamBlockBufferType = (PyObject *)&BamBlockBuffer_Type;
    Py_INCREF(BamBlockBufferType);
    if (PyModule_AddObject(m, "BamBlockBuffer", BamBlockBufferType) < 0) {
        return NULL;
    }

    if (PyType_Ready(&BamCigar_Type) < 0)
        return NULL;
    PyObject * BamCigarType = (PyObject *)&BamCigar_Type;
    Py_INCREF(BamCigarType);
    if (PyModule_AddObject(m, "Cigar", BamCigarType) < 0)
        return NULL;

    // Ready BamCigarIterType but do not expose it.
    if (PyType_Ready(&BamCigarIter_Type) < 0) {
        return NULL;
    }

    PyModule_AddIntMacro(m, BAM_CMATCH);
    PyModule_AddIntMacro(m, BAM_CINS);
    PyModule_AddIntMacro(m, BAM_CDEL);
    PyModule_AddIntMacro(m, BAM_CREF_SKIP);
    PyModule_AddIntMacro(m, BAM_CSOFT_CLIP);
    PyModule_AddIntMacro(m, BAM_CHARD_CLIP);
    PyModule_AddIntMacro(m, BAM_CPAD);
    PyModule_AddIntMacro(m, BAM_CEQUAL);
    PyModule_AddIntMacro(m, BAM_CDIFF);
    PyModule_AddIntMacro(m, BAM_CBACK);
    PyModule_AddIntMacro(m, BAM_CIGAR_SHIFT);

    PyModule_AddIntMacro(m, BAM_FPAIRED);
    PyModule_AddIntMacro(m, BAM_FPROPER_PAIR);
    PyModule_AddIntMacro(m, BAM_FUNMAP);
    PyModule_AddIntMacro(m, BAM_FREVERSE);
    PyModule_AddIntMacro(m, BAM_FMREVERSE);
    PyModule_AddIntMacro(m, BAM_FREAD1);
    PyModule_AddIntMacro(m, BAM_FREAD2);
    PyModule_AddIntMacro(m, BAM_FSECONDARY);
    PyModule_AddIntMacro(m, BAM_FQCFAIL);
    PyModule_AddIntMacro(m, BAM_FDUP);
    PyModule_AddIntMacro(m, BAM_FSUPPLEMENTARY);
    return m;
}
