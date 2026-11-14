#include "cdk_inchi_bridge.h"
#include "inchi_api.h"
#include <stdlib.h>
#include <string.h>

void cdk_inchi_free_stereo(CDKInChIStereoResult *result)
{
    if (!result) return;
    free(result->records);
    free(result->atom_symbols);
    free(result->hydrogen_parents);
    free(result->atom_hydrogens);
    free(result->atom_charges);
    free(result->bonds);
    memset(result, 0, sizeof(*result));
}

int cdk_inchi_read_stereo(const char *inchi, CDKInChIStereoResult *result)
{
    if (!inchi || !result) return 1;
    memset(result, 0, sizeof(*result));
    inchi_InputINCHI input = {0};
    inchi_OutputStruct output = {0};
    input.szInChI = (char *)inchi;
    input.szOptions = "";
    int status = GetStructFromINCHI(&input, &output);
    if ((status != inchi_Ret_OKAY && status != inchi_Ret_WARNING) ||
        output.num_atoms < 0 || output.num_stereo0D < 0)
    {
        FreeStructFromINCHI(&output);
        return 1;
    }
    result->atom_count = output.num_atoms;
    result->count = output.num_stereo0D;
    if (result->atom_count) result->atom_symbols = calloc((size_t)result->atom_count, 8);
    if (result->atom_count) result->hydrogen_parents = malloc((size_t)result->atom_count * sizeof(int));
    if (result->atom_count) result->atom_hydrogens = calloc((size_t)result->atom_count, sizeof(int));
    if (result->atom_count) result->atom_charges = calloc((size_t)result->atom_count, sizeof(int));
    for (int i = 0; i < output.num_atoms; ++i)
        for (int j = 0; j < output.atom[i].num_bonds; ++j)
            if (output.atom[i].neighbor[j] > i) ++result->bond_count;
    if (result->bond_count) result->bonds = calloc((size_t)result->bond_count * 3, sizeof(int));
    if (result->count) result->records = calloc((size_t)result->count * 7, sizeof(int));
    if ((result->atom_count && (!result->atom_symbols || !result->hydrogen_parents ||
                               !result->atom_hydrogens || !result->atom_charges)) ||
        (result->count && !result->records) || (result->bond_count && !result->bonds))
    {
        FreeStructFromINCHI(&output);
        cdk_inchi_free_stereo(result);
        return 2;
    }
    int bond_index = 0;
    for (int i = 0; i < result->atom_count; ++i)
    {
        strncpy(result->atom_symbols + i * 8, output.atom[i].elname, 7);
        inchi_Atom *atom = &output.atom[i];
        result->atom_hydrogens[i] = atom->num_iso_H[0];
        result->atom_charges[i] = atom->charge;
        for (int j = 0; j < atom->num_bonds; ++j)
            if (atom->neighbor[j] > i) {
                result->bonds[bond_index * 3] = i;
                result->bonds[bond_index * 3 + 1] = atom->neighbor[j];
                result->bonds[bond_index * 3 + 2] = atom->bond_type[j];
                ++bond_index;
            }
        result->hydrogen_parents[i] =
            !strcmp(atom->elname, "H") && atom->num_bonds == 1 && atom->isotopic_mass == 0
            ? atom->neighbor[0] : -1;
    }
    for (int i = 0; i < result->count; ++i)
    {
        int *record = result->records + i * 7;
        record[0] = output.stereo0D[i].type;
        record[1] = output.stereo0D[i].central_atom;
        for (int j = 0; j < 4; ++j) record[j + 2] = output.stereo0D[i].neighbor[j];
        record[6] = output.stereo0D[i].parity;
    }
    FreeStructFromINCHI(&output);
    return 0;
}
