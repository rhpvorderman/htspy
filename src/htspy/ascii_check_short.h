// Copyright (c) 2022 Leiden University Medical Center

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to 
// deal in the Software without restriction, including without limitation the 
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// This file is maintained and tested at
// https://github.com/rhpvorderman/ascii-check
// Please report bugs and feature requests there.

#include <stddef.h>
#include <stdint.h>

#define ASCII_MASK_8BYTE 0x8080808080808080ULL
#define ASCII_MASK_4BYTE 0x80808080U
#define ASCII_MASK_1BYTE 0x80U

/**
 * @brief Check if a string of given length only contains ASCII characters. 
 *        Optimized for short strings. (No memory alignment.)
 *
 * @param string A char pointer to the start of the string.
 * @param length The length of the string. This funtion does not check for 
 *               terminating NULL bytes.
 * @returns 1 if the string is ASCII-only
 */
static int
string_is_ascii(const char * string, size_t length) {
    size_t n = length;
    const char * char_ptr = string;
    const uint64_t *longword_ptr = (uint64_t *)char_ptr;
    while (n >= sizeof(uint64_t)) {
        if (*longword_ptr & ASCII_MASK_8BYTE){
            return 0;
        }
        longword_ptr += 1;
        n -= sizeof(uint64_t);
    }
    char_ptr = (char *)longword_ptr;
    while (n != 0) {
        if (*char_ptr & ASCII_MASK_1BYTE) {
            return 0;
        }
        char_ptr += 1;
        n -= 1;
    }
    return 1;
}
