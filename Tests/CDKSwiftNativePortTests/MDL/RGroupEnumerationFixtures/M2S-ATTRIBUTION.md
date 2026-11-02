# M2S regression backbones

103 unmodified CXSMILES backbone strings and source identifiers extracted from
`docling-project/MarkushGrapher-2-Datasets`, M2S split, revision
`e90868cdd3b62f1fe175d9ede9e0c65485c4a96c`.
Source: https://huggingface.co/datasets/docling-project/MarkushGrapher-2-Datasets
Authors: Tim Strohmeyer, Lucas Morin, Gerhard Ingmar Meijer, Valery Weber,
Ahmed Nassar, Peter W. J. Staar. MarkushGrapher (CVPR 2025) and
MarkushGrapher-2 (CVPR 2026), https://github.com/DS4SD/MarkushGrapher.
Licensed CC BY 4.0: https://creativecommons.org/licenses/by/4.0/
Changes: selected the id, image_name and cxsmiles fields from the original rows.
No structures corrected. Textual substituent definitions and images are not
included. These backbones test unresolved-variable detection, not membership.

## Reviewed definitions

`m2s-definitions.json` additionally retains all 103 source annotations and OCR
cell texts, with annotation hashes and source image identifiers. Five finite
annotation interpretations (5, 57, 72, 75, 101) add explicit fragment choices and
independently written membership checks. These are interpretations of the supplied
annotation, not a validation of complete patent claims or omitted source context.
The remaining 98 records are explicitly unreviewed and cannot be compiled.
No images or Java implementation are included. Nested variables, explicit
exclusions, dependencies, ambiguous position placeholders and enhanced stereo
are also covered by synthetic tests with declared semantics.
