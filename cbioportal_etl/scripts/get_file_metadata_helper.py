def get_file_metadata(table: str, ftype:str) -> dict[str, dict[str, dict[str, str]]]:
    """Convert table into dict for downsteam ETL use.

    Subsets table by ftypr (etl_file_type) column, then returns dict.
    cbio_project is primary key,  cbio_sample_name is secondary key.
    Tertiary keys are other attributes with string as values.
    """
    tbl_fh = open(table)
    head = next(tbl_fh)
    header = head.rstrip("\n").split("\t")
    cp_idx = header.index("cbio_project")
    ft_idx = header.index("etl_file_type")
    fn_idx = header.index("file_name")
    ct_idx = header.index("cbio_sample_name")
    cn_idx = header.index("cbio_matched_normal_name")
    kt_idx = header.index("affected_bs_id")
    kn_idx = header.index("reference_bs_id")
    meta_dict = {}
    for line in tbl_fh:
        info = line.rstrip("\n").split("\t")
        if info[ft_idx] == ftype:
            project = info[cp_idx]
            fname = info[fn_idx]
            cbio_tum_samp = info[ct_idx]
            cbio_norm_samp = info[cn_idx]
            kf_tum_samp = info[kt_idx]
            kf_norm_samp = info[kn_idx]
            if project not in meta_dict:
                meta_dict[project] = {}
            meta_dict[project][cbio_tum_samp] = {}
            meta_dict[project][cbio_tum_samp]["fname"] = fname
            meta_dict[project][cbio_tum_samp]["cbio_norm_id"] = cbio_norm_samp
            meta_dict[project][cbio_tum_samp]["kf_tum_id"] = kf_tum_samp
            meta_dict[project][cbio_tum_samp]["kf_norm_id"] = kf_norm_samp
    tbl_fh.close()
    return meta_dict
