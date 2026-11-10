#!/usr/bin/env python3
"""Compatibility entry point for the former converter script name.

New commands should use scripts/convert_mambacpg_to_merge.py. The public helper
functions and CLI remain available here for existing workflows.
"""

if __package__:
    from .convert_mambacpg_to_merge import (
        canonical_chrom,
        chrom_sort_key,
        clean_dna,
        infer_pos_base,
        validate_methylation_sample,
        main,
    )
else:
    from convert_mambacpg_to_merge import (
        canonical_chrom,
        chrom_sort_key,
        clean_dna,
        infer_pos_base,
        validate_methylation_sample,
        main,
    )


if __name__ == '__main__':
    main()
