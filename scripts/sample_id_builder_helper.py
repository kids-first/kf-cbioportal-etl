#!/usr/bin/env python3
import pdb


def build_samp_id(style, t_header, info):
    if style == 'cbttc_dna_std':
        # sample ID will consist of kf external_sample_id with .SUFFIX removed
        sid_idx = t_header.index('external_sample_id')
        samp_id = info[sid_idx].split('.')[0]
        return samp_id

    elif style == 'cbttc_norm_std':
        # sample ID will consist of kf external_sample_id with .SUFFIX removed
        sid_idx = t_header.index('external_sample_id')
        samp_id = info[sid_idx].split('.')[:-1][0]
        return samp_id

    elif style == 'cbttc_rna_std':
        # sample ID will consist of kf external_sample_id with .SUFFIX removed
        sid_idx = t_header.index('external_sample_id')
        ali_idx = t_header.index('external_aliquot_id')
        id_parts = info[sid_idx].split('.')
        ext_id_parts = id_parts[0].split('-')
        samp_id = ext_id_parts[0] + '-' + ext_id_parts[1] + '_' + info[ali_idx]
        return samp_id

    elif style == 'pnoc_dna_wgs':
        # sample ID will consist of kf external_sample_id with .SUFFIX removed and external_aliquot_id added
        sid_idx = t_header.index('external_sample_id')
        ali_idx = t_header.index('external_aliquot_id')
        samp_id = info[sid_idx].split('.')[0] + '-' + info[ali_idx]
        return samp_id

    elif style == 'pnoc_norm':
        # sample ID will consist of kf external_sample_id with .SUFFIX removed and external_aliquot_id added
        sid_idx = t_header.index('external_sample_id')
        ali_idx = t_header.index('external_aliquot_id')
        # Different styles for PNOC003 vs 008
        # 008 style has naked 7316 ID, tack on aliquot to differentiate from tumor
        samp_id = info[sid_idx] + '-' + info[ali_idx]

        if info[sid_idx][-1] == "S":
            # 003 style ending in .WGS or WXS
            samp_id = info[sid_idx].split('.')[:-1][0]
        return samp_id

    elif style == 'pnoc_rna_wes':
        # sample ID will consist of kf external_sample_id with .SUFFIX removed
        sid_idx = t_header.index('external_sample_id')
        ali_idx = t_header.index('external_aliquot_id')
        samp_id = info[sid_idx] + '_' + info[ali_idx]
        return samp_id

    elif style == 'maris_dna_wgs':
        # sample ID will consist of kf external_id, splitting on - and using only the last part
        sid_idx = t_header.index('external_id')
        samp_id = info[sid_idx].split("-")[2]
        return samp_id
    elif style == 'maris_norm_std':
        # sample ID will consist of kf external_id, splitting on - and using only the last part
        sid_idx = t_header.index('external_id')
        samp_id = info[sid_idx].split("-")[2]
        return samp_id

