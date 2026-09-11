with subset as (
	select distinct task_id, se.aliquot, tc.case_id , c2.submitter_case_id 
	FROM radiant_jdbc.public.sequencing_experiment se
JOIN radiant_jdbc.public.task_context tc on se.id = tc.sequencing_experiment_id 
JOIN radiant_jdbc.public.cases c2 on c2.id = tc.case_id 
where se.aliquot in('BS_0BZKABDK','BS_0WDJRMAB','BS_1WRWVYTD','BS_46H4TYCX','BS_4AM9C504','BS_6GEDTXKM','BS_6N7KRTKM','BS_91CV203N','BS_99MW9ECE','BS_CX6SQF8T','BS_DE3GMTKY','BS_EK6AB5ZC','BS_GR07XC1S','BS_KVBVPYNF','BS_M0TPSWZG','BS_M4KBSYZ2','BS_NBM9CPE4','BS_PT3TDRH5','BS_QB0QEZ1E','BS_SKWPQHER','BS_T5ZFAC8V','BS_TJ03GC0C','BS_TJW2Z7AB','BS_VR74JPGE','BS_WBK2D5XC','BS_X4VBHJSZ','BS_XMMTGNX2')
),
valid_locus as (
	select s.*, occ.`filter`, occ.locus_id, occ.dp , occ.ad_ref , occ.ad_alt, occ.info_mq, ggv.af, ggv.ac, ggv.an
	from subset s
	JOIN radiant.germline__snv__occurrence occ 
    ON s.task_id = occ.task_id 
	left outer join radiant.gnomad_genomes_v3 ggv on ggv.locus_id = occ.locus_id
	where ggv.af is NULL or ggv.af < 0.01 and occ.`filter` = "PASS"
)

SELECT DISTINCT 
	c.symbol as SYMBOL,  
	c.symbol as Hugo_Symbol,
	vl.aliquot  AS Matched_Norm_Sample_Barcode,
	'.' as Center,
	'GRCh38' as NCBI_Build,
    CONCAT('chr',v.chromosome) as Chromosome,
    v.`start` as Start_Position ,
    v.`start` as vcf_pos ,
    v.`end` as End_Position,
    '+' as Strand,
    c.consequences[1] as Consequence,
    v.variant_class as VARIANT_CLASS,
    v.reference as Reference_Allele,
    v.alternate as Match_Norm_Seq_Allele1,
    'Germline' as Mutation_Status,
    case 
   		when v.is_canonical then "YES"
   		else "NO"
    end as 'CANONICAL',
    c.dna_change as HGVSc,
    c.aa_change as HGVSp,
    c.transcript_id as Transcript_ID,
    c.transcript_id as Feature,
    vl.dp as n_depth,
    vl.ad_ref as n_ref_count,
    vl.ad_alt as n_alt_count,
    c.biotype as BIOTYPE,
    concat(c.exon_rank, "/", c.exon_total) as Exon_Number,
    concat(c.exon_rank, "/", c.exon_total) as 'EXON',
    c.is_picked as PICK,
    vl.`filter` as `FILTER`,
    c.vep_impact AS IMPACT,
    v.hgvsg as HGVSg,
    vl.ac as gnomad_3_1_1_AC,
    vl.an as gnomad_3_1_1_AN,
    vl.af as gnomad_3_1_1_AF,
    vl.info_mq as MQ,
    'DRAGEN' as CAL

FROM valid_locus vl
JOIN radiant.snv__consequence c 
    ON vl.locus_id = c.locus_id
JOIN radiant.snv__variant v 
    ON vl.locus_id = v.locus_id
where  c.vep_impact not in ('MODIFIER', 'LOW') and c.is_picked
order by vl.aliquot, v.hgvsg ASC
