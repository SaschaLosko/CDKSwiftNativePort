#ifndef IUPAC_INCHI_API_H
#define IUPAC_INCHI_API_H

/* The public InChI and IXA API of the vendored IUPAC InChI library, for packages
   that drive the library directly instead of through cdk_inchi_bridge.h. */
#include "../official/INCHI-1-SRC/INCHI_BASE/src/inchi_api.h"

/* IXA_ATOMID_IMPLICIT_H is a macro cast that Swift cannot import. */
static inline IXA_ATOMID iupac_inchi_implicit_h(void) {
    return IXA_ATOMID_IMPLICIT_H;
}

#endif
