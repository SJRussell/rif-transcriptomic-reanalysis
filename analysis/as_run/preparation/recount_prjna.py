from pathlib import Path
import subprocess, json, datetime, hashlib

root = Path('<server_workdir>')
out = root / 'PRJNA576073'
out.mkdir(parents=True, exist_ok=True)
bams = sorted(Path('<server_data>/PRJNA576073_rnaseq/STAR_PRJNA576073').glob('*_Aligned.sortedByCoord.out.bam'))
assert len(bams) == 6
gtf = Path('/data/human_genome_T2T/ncbi_dataset/T2T.gtf')
assert gtf.is_file()
manifest = {'started_utc': datetime.datetime.utcnow().isoformat(), 'bams': [{'path': str(p), 'bytes': p.stat().st_size} for p in bams], 'annotation': str(gtf), 'commands': []}
for strand in [0, 2]:
    target = out / ('gene_counts_s%d.txt' % strand)
    if target.exists():
        raise RuntimeError('Refusing to overwrite count output: '+str(target))
    cmd = ['nice', '-n', '10', '/usr/bin/featureCounts', '-T', '4', '-t', 'exon', '-g', 'gene_id', '-s', str(strand), '-p', '--countReadPairs', '-a', str(gtf), '-o', str(target)] + [str(p) for p in bams]
    manifest['commands'].append(cmd)
    (out/'run_manifest.json').write_text(json.dumps(manifest, indent=2))
    with (out/('featureCounts_s%d.log' % strand)).open('w') as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    print('Completed strand', strand, flush=True)
manifest['completed_utc'] = datetime.datetime.utcnow().isoformat()
manifest['annotation_sha256'] = hashlib.sha256(gtf.read_bytes()).hexdigest()
(out/'run_manifest.json').write_text(json.dumps(manifest, indent=2))
print('RECOUNT_COMPLETE', flush=True)
