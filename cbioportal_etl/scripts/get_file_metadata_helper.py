def get_file_metadata(table: str, ftype:str) -> dict[str, dict[str, dict[str, str]]]:
    """Convert table into dict for downsteam ETL use.

    Subsets table by ftype (etl_file_type) column, then returns dict.
    cbio_project is primary key,  cbio_sample_name is secondary key.
    Tertiary keys are other attributes with string as values.
    """
    tbl_fh = open(table)
    head = next(tbl_fh)
    header = head.rstrip("\n").split("\t")
    project_idx = header.index("cbio_project")
    etl_ft_idx = header.index("etl_file_type")
    manifest_ft_idx = header.index("file_type")
    fname_idx = header.index("file_name")
    cbio_sample_idx = header.index("cbio_sample_name")
    cbio_normal_idx = header.index("cbio_matched_normal_name")
    kf_affected_idx = header.index("affected_bs_id")
    kf_reference_idx = header.index("reference_bs_id")
    meta_dict = {}
    for line in tbl_fh:
        info = line.rstrip("\n").split("\t")
        if info[etl_ft_idx] == ftype:
            project = info[project_idx]
            fname = info[fname_idx]
            cbio_tum_samp = info[cbio_sample_idx]
            cbio_norm_samp = info[cbio_normal_idx]
            kf_tum_samp = info[kf_affected_idx]
            kf_norm_samp = info[kf_reference_idx]
            manifest_ftype = info[manifest_ft_idx]
            if project not in meta_dict:
                meta_dict[project] = {}
            meta_dict[project][cbio_tum_samp] = {}
            meta_dict[project][cbio_tum_samp]["fname"] = fname
            meta_dict[project][cbio_tum_samp]["cbio_norm_id"] = cbio_norm_samp
            meta_dict[project][cbio_tum_samp]["kf_tum_id"] = kf_tum_samp
            meta_dict[project][cbio_tum_samp]["kf_norm_id"] = kf_norm_samp
            meta_dict[project][cbio_tum_samp]["manifest_ftype"] = manifest_ftype
    tbl_fh.close()
    return meta_dict
