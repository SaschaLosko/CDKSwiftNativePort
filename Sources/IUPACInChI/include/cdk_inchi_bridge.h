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
int cdk_inchi_key_from_standard_inchi(const char *inchi, char *key_buffer, int key_buffer_length);
void cdk_inchi_free_result(CDKInChIResult *result);
const char *cdk_inchi_version(void);

#ifdef __cplusplus
}
#endif

#endif
