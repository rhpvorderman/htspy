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

typedef struct {
    PyObject_HEAD
    PyObject * raw;
    uint32_t * cigar;
    Py_ssize_t n_cigar_op;

} BamCigar;

static PyTypeObject BamCigar_Type;  // Forward declaration

static void 
BamCigar_dealloc(BamCigar *self) {
    Py_CLEAR(self->raw);
}

/**
 * @brief Creates a new BamCigar without checks. For internal calls
 * 
 * @param bytes a PyBytesObject (not checked)
 * @param n_cigar_op (number of cigar units. Not checked with length of Bytes object.)
 * @return PyObject* 
 */
static PyObject *
BamCigar_FromBytesAndSize(PyObject * bytes, Py_ssize_t n_cigar_op) {
    BamCigar *cigar = PyObject_New(BamCigar, &BamCigar_Type);
    Py_INCREF(bytes);
    cigar->raw = bytes;
    cigar->cigar = (uint32_t *)PyBytes_AS_STRING(bytes);
    cigar->n_cigar_op;
}

static PyTypeObject BamCigar_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "_pybam.BamCigar",
    .tp_basicsize = sizeof(BamCigar),
    .tp_dealloc = (destructor)BamCigar_dealloc
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
    PyObject * cigar;
    PyObject * seq;
    PyObject * qual;
    PyObject * tags;
} BamRecord;

# define BAM_PROPERTIES_STRUCT_START offsetof(BamRecord, block_size)
# define BAM_PROPERTIES_STRUCT_SIZE 36

static void
BamRecord_dealloc(BamRecord *self) {
    Py_CLEAR(self->read_name);
    Py_CLEAR(self->cigar);
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
    self->cigar = PyBytes_FromStringAndSize("", 0);
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
    {"_read_name", T_OBJECT_EX, offsetof(BamRecord, read_name), READONLY},
    {"_cigar", T_OBJECT_EX, offsetof(BamRecord, cigar), READONLY},
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

PyDoc_STRVAR(BamRecord_query_name_doc,
"The name of the aligned read as a string.\n"
"WARNING: this attribute is a property that converts 'read_name' \n"
"to ASCII For faster access use the 'read_name' attribute which \n"
"is an ASCII-encoded bytes object.");

static PyObject * 
BamRecord_get_query_name(BamRecord * self, void* closure) 
{
    return PyUnicode_FromEncodedObject(self->read_name, "ascii", "strict");
}

static int 
BamRecord_set_query_name(BamRecord * self, PyObject * new_qname, void* closure) 
{
    PyObject * new_read_name = PyUnicode_AsASCIIString(new_qname);
    if (new_read_name == NULL)
        return -1;
    Py_ssize_t read_name_size = PyBytes_GET_SIZE(new_read_name);
    if (read_name_size > 254) {
        PyErr_SetString(PyExc_ValueError, 
            "read_name may not be larger than 254 characters.");
        Py_DecRef(new_read_name);
        return -1;
    }
    PyObject * old_read_name = self->read_name;
    self->read_name = new_read_name;
    Py_DECREF(old_read_name);
    uint8_t old_l_read_name = self->l_read_name;
    self->l_read_name = (uint8_t)read_name_size + 1;
    self->block_size = self->block_size + self->l_read_name - old_l_read_name;
    return 0;
}

PyDoc_STRVAR(BamRecord_read_name_doc,
"The name of the aligned read as an ASCII encoded bytes object.\n");

static PyObject * 
BamRecord_get_read_name(BamRecord * self, void* closure) {
    Py_INCREF(self->read_name);
    return self->read_name;
}

static int 
BamRecord_set_read_name(BamRecord * self, PyObject * new_read_name, void* closure) 
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
    {"query_name", (getter)BamRecord_get_query_name, (setter)BamRecord_set_query_name, 
     BamRecord_query_name_doc, NULL},
    {"read_name", (getter)BamRecord_get_read_name, (setter)BamRecord_set_read_name,
     BamRecord_read_name_doc, NULL},
    {"tags", (getter)BamRecord_get_tags, (setter)BamRecord_set_tags,
     BamRecord_tags_doc, NULL},
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
    
    memcpy(bam_bytes, (char *)self + BAM_PROPERTIES_STRUCT_START,
         BAM_PROPERTIES_STRUCT_SIZE);
    Py_ssize_t cursor = BAM_PROPERTIES_STRUCT_SIZE;
    
    Py_ssize_t read_name_size = PyBytes_GET_SIZE(self->read_name);
    memcpy(bam_bytes + cursor, PyBytes_AS_STRING(self->read_name), read_name_size);
    cursor += read_name_size;

    // Terminate read_name with NULL byte
    bam_bytes[cursor] = 0; cursor += 1;

    Py_ssize_t cigar_size = PyBytes_GET_SIZE(self->cigar);
    memcpy(bam_bytes + cursor, PyBytes_AS_STRING(self->cigar), cigar_size);
    cursor += cigar_size;

    Py_ssize_t seq_size = PyBytes_GET_SIZE(self->seq);
    memcpy(bam_bytes + cursor, PyBytes_AS_STRING(self->seq), seq_size);
    cursor += seq_size;

    Py_ssize_t qual_size = PyBytes_GET_SIZE(self->qual);
    memcpy(bam_bytes + cursor, PyBytes_AS_STRING(self->qual), qual_size);
    cursor += qual_size;

    Py_ssize_t tag_size = PyBytes_GET_SIZE(self->tags);
    memcpy(bam_bytes + cursor, PyBytes_AS_STRING(self->tags), tag_size);

    return ret_val;
}


static PyMethodDef BamRecord_methods[] = {
    BAMRECORD_TO_BYTES_METHODDEF,
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
    bam_record->cigar = NULL;
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
    bam_record->cigar = PyBytes_FromStringAndSize(
        self->buf + self->pos, cigar_length);
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
    .tp_iter = (PyObject * (*)(PyObject *))BamIterator_iter,
    .tp_iternext = (PyObject * (*)(PyObject *))BamIterator_iternext
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

    if (PyType_Ready(&BamCigar_Type) < 0)
        return NULL;
    PyObject * BamCigarType = (PyObject *)&BamCigar_Type;
    Py_INCREF(BamCigarType);
    if (PyModule_AddObject(m, "BamCigar", BamCigarType) < 0)
        return NULL;

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
