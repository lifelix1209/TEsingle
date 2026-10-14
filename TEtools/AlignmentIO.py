import shlex

import pysam


def validate_alignment(filename, allow_truncated=False):
    """Catch the STAR setting that makes a multimapper look uniquely aligned."""
    with pysam.AlignmentFile(filename, check_sq=False) as bam:
        for program in bam.header.to_dict().get('PG', []):
            if program.get('PN') != 'STAR':
                continue
            tokens = shlex.split(program.get('CL', ''))
            def option(name, default):
                return int(tokens[tokens.index(name) + 1]) if name in tokens else default
            reported = option('--outSAMmultNmax', -1)
            permitted = option('--outFilterMultimapNmax', 10)
            if 0 < reported < permitted and not allow_truncated:
                raise ValueError(
                    'STAR BAM has truncated multimapping alignments '
                    '(--outSAMmultNmax %s, --outFilterMultimapNmax %s). '
                    'Realign FASTQs with --outSAMmultNmax -1. '
                    'The missing alignments cannot be recovered from this BAM. '
                    '--allow-truncated-bam is for diagnostic comparisons only.' % (reported, permitted))


def annotation_direction(stranded, read_is_reverse):
    if stranded == 'no':
        return '.'
    else:
        direction = [1, -1][stranded == 'reverse']
        algn_dir = [1, -1][read_is_reverse]
        return ['.', '+', '-'][algn_dir * direction]


def read_in_alignment(filename, geneIdx, teIdx, stranded):
    # given .sam or .bam file, returns dictionary of {barcodes:UMIs} and {UMIs:{reads:{alignments}}}
    # alignments are in format [chromosome, [0-indexed blocks in tuples], direction]

    psave = pysam.set_verbosity(0)  # disable htslib verbosity to avoid "no index file" warning
    samfile = pysam.AlignmentFile(filename, check_sq=False)
    pysam.set_verbosity(psave)  # re-enable htslib verbosity

    transc_dict = {}
    err_reads = 0

    try:
        while 1:
            aligned_read = next(samfile)
            if aligned_read.is_unmapped or aligned_read.is_duplicate or aligned_read.is_qcfail:
                continue

            rbc = False
            rumi = False

            try:
                rbc = aligned_read.get_tag('CB')
                if rbc == '-':
                    rbc = False
            except:
                err_reads += 1

            if not rbc:
                continue

            try:
                rumi = aligned_read.get_tag('UR')
                if rumi == '-':
                    rumi = False

            except:
                err_reads += 1

            if not (rbc and rumi):
                continue

            if rbc not in transc_dict:
                transc_dict[rbc] = {}

            if rumi not in transc_dict[rbc]:
                transc_dict[rbc][rumi] = {}

            cur_read_name = aligned_read.query_name

            if cur_read_name not in transc_dict[rbc][rumi]:
                transc_dict[rbc][rumi][cur_read_name] = []

            algn = aligned_read.get_blocks()  # 0-indexed
            chrom = aligned_read.reference_name
            direction = annotation_direction(stranded, aligned_read.is_reverse)

            itv_list = []

            for block in algn:
                itv_list.append([chrom, block[0]+1, block[1], direction])

            tmp_genes = geneIdx.Gene_annotation(itv_list)
            tmp_tes = teIdx.TE_annotation(itv_list)
            transc_dict[rbc][rumi][cur_read_name].append(tuple([tmp_genes, tmp_tes]))

    except StopIteration:
        pass

    finally:
        samfile.close()

    return transc_dict
