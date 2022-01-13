#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"         // PyMemberDef
#include <stdint.h>

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
    PyObject * read_name;
    PyObject * cigar;
    PyObject * seq;
    PyObject * qual;
    PyObject * tags;
} BamRecord;

static void
BamRecord_dealloc(BamRecord *self) {
    Py_CLEAR(self->read_name);
    Py_CLEAR(self->cigar);
    Py_CLEAR(self->seq);
    Py_CLEAR(self->qual);
    Py_CLEAR(self->tags);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyMemberDef BamRecord_members[] = {
    {"read_name", T_OBJECT_EX, offsetof(BamRecord, read_name), 0},
    {"flag", T_USHORT, offsetof(BamRecord, flag), READONLY},
    {"pos", T_INT, offsetof(BamRecord, pos), READONLY},
    {"mapq", T_UBYTE, offsetof(BamRecord, mapq), READONLY},
    {"cigar", T_OBJECT_EX, offsetof(BamRecord, cigar), READONLY},
    {"next_pos", T_INT, offsetof(BamRecord, next_pos), READONLY},
    {"pnext", T_INT, offsetof(BamRecord, next_pos), READONLY}, // SAM name for next pos.
    {"tlen", T_INT, offsetof(BamRecord, tlen), READONLY},
    {"seq", T_OBJECT_EX, offsetof(BamRecord, seq), READONLY},
    {"qual", T_OBJECT_EX, offsetof(BamRecord, qual), READONLY},
    {"tags", T_OBJECT_EX, offsetof(BamRecord, tags), 0},
    {NULL}
};

static PyObject * BamRecord_get_qname(BamRecord * self, void* closure) {
    return PyUnicode_FromEncodedObject(self->read_name, "ascii", "strict");
}

static int BamRecord_set_qname(BamRecord * self, PyObject * new_qname, void* closure) {
    PyObject * new_read_name = PyUnicode_AsASCIIString(new_qname);
    if (new_read_name == NULL)
        return -1;
    self->read_name = new_read_name;
    return 0;
}

PyDoc_STRVAR(BamRecord_qname_doc,
"The name of the aligned read as a string.\n"
"WARNING: this attribute is a property that converts 'read_name' \n"
"To ASCII For faster access use the 'read_name' attribute which \n"
"is an ASCII-encoded bytes object.");

static PyGetSetDef BamRecord_properties[] = {
    {"qname", (getter)BamRecord_get_qname, (setter)BamRecord_set_qname, BamRecord_qname_doc, NULL},
    {NULL}
};

static PyMethodDef BamRecord_methods[] = {
    {NULL}
};

PyDoc_STRVAR(BamRecord__doc__, 
"An object that represents an alignment record from a BAM file.");

static PyTypeObject BamRecord_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_bamrecord.BamRecord",             /* tp_name */
    sizeof(BamRecord),                  /* tp_basicsize */
    0,                                  /* tp_itemsize */
    (destructor)BamRecord_dealloc,      /* tp_dealloc */
    0,                                  /* tp_vectorcall_offset */
    0,                                  /* tp_getattr */
    0,                                  /* tp_setattr */
    0,                                  /* tp_as_async */
    0,                                  /* tp_repr */
    0,                                  /* tp_as_number */
    0,                                  /* tp_as_sequence */
    0,                                  /* tp_as_mapping */
    0,                                  /* tp_hash  */
    0,                                  /* tp_call */
    0,                                  /* tp_str */
    0,                                  /* tp_getattro */
    0,                                  /* tp_setattro */
    0,                                  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                 /* tp_flags */
    BamRecord__doc__,                   /* tp_doc */
    0,                                  /* tp_traverse */
    0,                                  /* tp_clear */
    0,                                  /* tp_richcompare */
    0,                                  /* tp_weaklistoffset */
    0,                                  /* tp_iter */
    0,                                  /* tp_iternext */
    BamRecord_methods,                  /* tp_methods */
    BamRecord_members,                  /* tp_members */
    BamRecord_properties,               /* tp_getset */
    0,                                  /* tp_base */
    0,                                  /* tp_dict */
    0,                                  /* tp_descr_get */
    0,                                  /* tp_descr_set */
    0,                                  /* tp_dictoffset */
    0,                                  /* tp_init */
    0,                                  /* tp_alloc */
    0,                                  /* tp_new */
};


static struct PyModuleDef _sequence_module = {
    PyModuleDef_HEAD_INIT,
    "_sequence",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,
    NULL  /* module methods */
};


PyMODINIT_FUNC
PyInit__sequence(void)
{
    PyObject *m;

    m = PyModule_Create(&_sequence_module);
    if (m == NULL)
        return NULL;

    PyObject * BamRecordType = (PyObject *)&BamRecord_Type;
    Py_INCREF(BamRecordType);
    if (PyModule_AddObject(m, "BamRecord", BamRecordType) < 0) {
        return NULL;
    }
    return m;
}
