#include "cdk_inchi_bridge.h"

#include <stdlib.h>
#include <string.h>

#include "inchi_api.h"
#include "mode.h"

static char *cdk_inchi_strdup(const char *source)
{
    if (!source)
    {
        return NULL;
    }

    size_t length = strlen(source);
    char *copy = (char *) malloc(length + 1);
    if (!copy)
    {
        return NULL;
    }

    memcpy(copy, source, length + 1);
    return copy;
}

static void cdk_inchi_reset_result(CDKInChIResult *result)
{
    if (!result)
    {
        return;
    }

    result->inchi = NULL;
    result->inchi_key = NULL;
    result->aux_info = NULL;
    result->message = NULL;
    result->log = NULL;
    result->return_code = 0;
    result->key_return_code = 0;
}

static int cdk_inchi_finish_output(int return_code,
                                   inchi_Output *output,
                                   int is_standard,
                                   CDKInChIResult *result)
{
    result->return_code = return_code;
    result->inchi = cdk_inchi_strdup(output->szInChI);
    result->aux_info = cdk_inchi_strdup(output->szAuxInfo);
    result->message = cdk_inchi_strdup(output->szMessage);
    result->log = cdk_inchi_strdup(output->szLog);

    if (output->szInChI && output->szInChI[0] != '\0')
    {
        char key[28];
        memset(key, 0, sizeof(key));
        result->key_return_code = is_standard
            ? GetStdINCHIKeyFromStdINCHI(output->szInChI, key)
            : GetINCHIKeyFromINCHI(output->szInChI, 0, 0, key, NULL, NULL);
        if (result->key_return_code == INCHIKEY_OK)
        {
            result->inchi_key = cdk_inchi_strdup(key);
        }
    }

    if ((return_code == inchi_Ret_OKAY || return_code == inchi_Ret_WARNING ||
         return_code == mol2inchi_Ret_OKAY || return_code == mol2inchi_Ret_WARNING) &&
        result->inchi && result->inchi[0] != '\0' &&
        result->inchi_key && result->inchi_key[0] != '\0')
    {
        return 0;
    }

    return return_code == 0 ? 2 : return_code;
}

int cdk_inchi_generate_standard_molfile(const char *molfile_text, CDKInChIResult *result)
{
    inchi_Output output;
    char options[] = "";

    if (!molfile_text || !result)
    {
        return 2;
    }

    cdk_inchi_reset_result(result);
    memset(&output, 0, sizeof(output));

    int return_code = MakeINCHIFromMolfileText(molfile_text, options, &output);
    int bridge_code = cdk_inchi_finish_output(return_code, &output, 1, result);

    FreeINCHI(&output);
    return bridge_code;
}

int cdk_inchi_generate_fixed_h_molfile(const char *molfile_text, CDKInChIResult *result)
{
    inchi_Output output;
    char options[] = "-FixedH";

    if (!molfile_text || !result)
    {
        return 2;
    }

    cdk_inchi_reset_result(result);
    memset(&output, 0, sizeof(output));

    int return_code = MakeINCHIFromMolfileText(molfile_text, options, &output);
    int bridge_code = cdk_inchi_finish_output(return_code, &output, 0, result);

    FreeINCHI(&output);
    return bridge_code;
}

static int cdk_inchi_add_bond(inchi_Atom *atoms, int atom_count, int from, int to, int order)
{
    if (!atoms || from < 0 || to < 0 || from >= atom_count || to >= atom_count || from == to)
    {
        return 2;
    }

    if (atoms[from].num_bonds >= MAXVAL || atoms[to].num_bonds >= MAXVAL)
    {
        return 2;
    }

    int normalized_order;
    switch (order)
    {
    case INCHI_BOND_TYPE_SINGLE:
    case INCHI_BOND_TYPE_DOUBLE:
    case INCHI_BOND_TYPE_TRIPLE:
    case INCHI_BOND_TYPE_ALTERN:
        normalized_order = order;
        break;
    default:
        normalized_order = INCHI_BOND_TYPE_SINGLE;
        break;
    }

    int from_slot = atoms[from].num_bonds++;
    atoms[from].neighbor[from_slot] = (AT_NUM) to;
    atoms[from].bond_type[from_slot] = (S_CHAR) normalized_order;
    atoms[from].bond_stereo[from_slot] = 0;

    int to_slot = atoms[to].num_bonds++;
    atoms[to].neighbor[to_slot] = (AT_NUM) from;
    atoms[to].bond_type[to_slot] = (S_CHAR) normalized_order;
    atoms[to].bond_stereo[to_slot] = 0;

    return 0;
}

static int cdk_inchi_generate_atoms(int atom_count,
                                    const char *atom_symbols,
                                    int atom_symbol_stride,
                                    const int *atom_charges,
                                    const int *atom_isotopes,
                                    const int *atom_implicit_hydrogens,
                                    int bond_count,
                                    const int *bond_from,
                                    const int *bond_to,
                                    const int *bond_order,
                                    int stereo_count,
                                    const int *stereo_centers,
                                    const int *stereo_neighbors,
                                    const int *stereo_parities,
                                    int fixed_h,
                                    CDKInChIResult *result)
{
    inchi_Output output;
    inchi_Input input;
    char standard_options[] = "";
    char fixed_h_options[] = "-FixedH";

    if (!result || atom_count <= 0 || !atom_symbols || atom_symbol_stride <= 0 ||
        (bond_count > 0 && (!bond_from || !bond_to || !bond_order)) ||
        (stereo_count > 0 && (!stereo_centers || !stereo_neighbors || !stereo_parities)))
    {
        return 2;
    }

    cdk_inchi_reset_result(result);
    memset(&output, 0, sizeof(output));
    memset(&input, 0, sizeof(input));

    inchi_Atom *atoms = (inchi_Atom *) calloc((size_t) atom_count, sizeof(inchi_Atom));
    inchi_Stereo0D *stereo0D = NULL;
    if (!atoms)
    {
        return 3;
    }

    if (stereo_count > 0)
    {
        stereo0D = (inchi_Stereo0D *) calloc((size_t) stereo_count, sizeof(inchi_Stereo0D));
        if (!stereo0D)
        {
            free(atoms);
            return 3;
        }
    }

    for (int i = 0; i < atom_count; i++)
    {
        const char *symbol = atom_symbols + ((size_t) i * (size_t) atom_symbol_stride);
        strncpy(atoms[i].elname, symbol, ATOM_EL_LEN - 1);
        atoms[i].elname[ATOM_EL_LEN - 1] = '\0';
        atoms[i].charge = atom_charges ? (S_CHAR) atom_charges[i] : 0;
        atoms[i].isotopic_mass = atom_isotopes ? (AT_NUM) atom_isotopes[i] : 0;
        atoms[i].radical = 0;
        for (int h = 0; h < NUM_H_ISOTOPES + 1; h++)
        {
            atoms[i].num_iso_H[h] = 0;
        }
        atoms[i].num_iso_H[0] = atom_implicit_hydrogens ? (S_CHAR) atom_implicit_hydrogens[i] : -1;
    }

    for (int i = 0; i < bond_count; i++)
    {
        int status = cdk_inchi_add_bond(atoms, atom_count, bond_from[i], bond_to[i], bond_order[i]);
        if (status != 0)
        {
            free(stereo0D);
            free(atoms);
            return status;
        }
    }

    for (int i = 0; i < stereo_count; i++)
    {
        int center = stereo_centers[i];
        if (center < 0 || center >= atom_count)
        {
            free(stereo0D);
            free(atoms);
            return 2;
        }
        stereo0D[i].central_atom = (AT_NUM) center;
        stereo0D[i].type = INCHI_StereoType_Tetrahedral;
        stereo0D[i].parity = (S_CHAR) stereo_parities[i];
        for (int j = 0; j < 4; j++)
        {
            int neighbor = stereo_neighbors[i * 4 + j];
            if (neighbor < 0 || neighbor >= atom_count)
            {
                free(stereo0D);
                free(atoms);
                return 2;
            }
            stereo0D[i].neighbor[j] = (AT_NUM) neighbor;
        }
    }

    input.atom = atoms;
    input.stereo0D = stereo0D;
    input.szOptions = fixed_h ? fixed_h_options : standard_options;
    input.num_atoms = (AT_NUM) atom_count;
    input.num_stereo0D = (AT_NUM) stereo_count;

    int return_code = fixed_h ? GetINCHI(&input, &output) : GetStdINCHI(&input, &output);
    int bridge_code = cdk_inchi_finish_output(return_code, &output, fixed_h ? 0 : 1, result);

    if (fixed_h)
    {
        FreeINCHI(&output);
    }
    else
    {
        FreeStdINCHI(&output);
    }
    free(stereo0D);
    free(atoms);
    return bridge_code;
}

int cdk_inchi_generate_standard_atoms(int atom_count,
                                      const char *atom_symbols,
                                      int atom_symbol_stride,
                                      const int *atom_charges,
                                      const int *atom_isotopes,
                                      const int *atom_implicit_hydrogens,
                                      int bond_count,
                                      const int *bond_from,
                                      const int *bond_to,
                                      const int *bond_order,
                                      int stereo_count,
                                      const int *stereo_centers,
                                      const int *stereo_neighbors,
                                      const int *stereo_parities,
                                      CDKInChIResult *result)
{
    return cdk_inchi_generate_atoms(atom_count,
                                    atom_symbols,
                                    atom_symbol_stride,
                                    atom_charges,
                                    atom_isotopes,
                                    atom_implicit_hydrogens,
                                    bond_count,
                                    bond_from,
                                    bond_to,
                                    bond_order,
                                    stereo_count,
                                    stereo_centers,
                                    stereo_neighbors,
                                    stereo_parities,
                                    0,
                                    result);
}

int cdk_inchi_generate_fixed_h_atoms(int atom_count,
                                     const char *atom_symbols,
                                     int atom_symbol_stride,
                                     const int *atom_charges,
                                     const int *atom_isotopes,
                                     const int *atom_implicit_hydrogens,
                                     int bond_count,
                                     const int *bond_from,
                                     const int *bond_to,
                                     const int *bond_order,
                                     int stereo_count,
                                     const int *stereo_centers,
                                     const int *stereo_neighbors,
                                     const int *stereo_parities,
                                     CDKInChIResult *result)
{
    return cdk_inchi_generate_atoms(atom_count,
                                    atom_symbols,
                                    atom_symbol_stride,
                                    atom_charges,
                                    atom_isotopes,
                                    atom_implicit_hydrogens,
                                    bond_count,
                                    bond_from,
                                    bond_to,
                                    bond_order,
                                    stereo_count,
                                    stereo_centers,
                                    stereo_neighbors,
                                    stereo_parities,
                                    1,
                                    result);
}

int cdk_inchi_key_from_standard_inchi(const char *inchi, char *key_buffer, int key_buffer_length)
{
    if (!inchi || !key_buffer || key_buffer_length < 28)
    {
        return INCHIKEY_UNKNOWN_ERROR;
    }

    memset(key_buffer, 0, (size_t) key_buffer_length);
    return GetStdINCHIKeyFromStdINCHI(inchi, key_buffer);
}

void cdk_inchi_free_result(CDKInChIResult *result)
{
    if (!result)
    {
        return;
    }

    free(result->inchi);
    free(result->inchi_key);
    free(result->aux_info);
    free(result->message);
    free(result->log);
    cdk_inchi_reset_result(result);
}

const char *cdk_inchi_version(void)
{
    return INCHI_VERSION;
}
