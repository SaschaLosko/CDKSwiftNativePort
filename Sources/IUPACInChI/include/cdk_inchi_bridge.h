#ifndef CDK_INCHI_BRIDGE_H
#define CDK_INCHI_BRIDGE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct CDKInChIResult {
    char *inchi;
    char *inchi_key;
    char *aux_info;
    char *message;
    char *log;
    int return_code;
    int key_return_code;
} CDKInChIResult;

int cdk_inchi_generate_standard_molfile(const char *molfile_text, CDKInChIResult *result);
int cdk_inchi_generate_fixed_h_molfile(const char *molfile_text, CDKInChIResult *result);
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
                                      CDKInChIResult *result);
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
                                     CDKInChIResult *result);
int cdk_inchi_key_from_standard_inchi(const char *inchi, char *key_buffer, int key_buffer_length);
void cdk_inchi_free_result(CDKInChIResult *result);
const char *cdk_inchi_version(void);

#ifdef __cplusplus
}
#endif

#endif
