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
    result->return_code = return_code;

    result->inchi = cdk_inchi_strdup(output.szInChI);
    result->aux_info = cdk_inchi_strdup(output.szAuxInfo);
    result->message = cdk_inchi_strdup(output.szMessage);
    result->log = cdk_inchi_strdup(output.szLog);

    if (output.szInChI && output.szInChI[0] != '\0')
    {
        char key[28];
        memset(key, 0, sizeof(key));
        result->key_return_code = GetStdINCHIKeyFromStdINCHI(output.szInChI, key);
        if (result->key_return_code == INCHIKEY_OK)
        {
            result->inchi_key = cdk_inchi_strdup(key);
        }
    }

    FreeINCHI(&output);

    if ((return_code == mol2inchi_Ret_OKAY || return_code == mol2inchi_Ret_WARNING) &&
        result->inchi && result->inchi[0] != '\0' &&
        result->inchi_key && result->inchi_key[0] != '\0')
    {
        return 0;
    }

    return return_code == 0 ? 2 : return_code;
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
