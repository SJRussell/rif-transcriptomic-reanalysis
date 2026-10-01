"""Run inside the existing transcriptomics container under SLURM; source references read-only."""
from pathlib import Path
import os,sys,json,urllib.request,hashlib,gzip,subprocess,datetime,time
ROOT=Path('<cluster_scratch>/subtypes')
SOURCE=Path('<cluster_project>/ref/gencode_human')
mode=sys.argv[1];samples=sorted(json.loads((ROOT/'reference/Xu2026_sample_metadata.json').read_text()),key=lambda r:r['run_accession'])
assert len(samples)==22 and len({s['title'] for s in samples})==22
def digest(p):
 with p.open('rb') as f:
  h=hashlib.sha256()
  for b in iter(lambda:f.read(1024*1024),b''):
   h.update(b)
 return h.hexdigest()
def run(cmd,log):
 with log.open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
def save(p,data):p.write_text(json.dumps(data,indent=2))
if mode=='download':
 s=samples[int(os.environ['SLURM_ARRAY_TASK_ID'])];out=ROOT/'raw';out.mkdir(exist_ok=True);p=out/(s['run_accession']+'.fastq.gz')
 if not p.exists():
  u='https://'+s['fastq_ftp'];tmp=p.with_suffix('.gz.part')
  for attempt in range(3):
   try:
    with urllib.request.urlopen(u,timeout=90) as r,tmp.open('wb') as f:
     while True:
      b=r.read(1024*1024)
      if not b:break
      f.write(b)
    tmp.replace(p);break
   except Exception:
    if attempt==2:raise
    time.sleep(5)
 h=hashlib.md5()
 with p.open('rb') as f:
  while True:
   b=f.read(1024*1024)
   if not b:break
   h.update(b)
 assert h.hexdigest()==s['fastq_md5'] and p.stat().st_size==int(s['fastq_bytes'])
 save(out/(s['run_accession']+'.verified.json'),dict(sample=s['title'],condition='RIF' if s['title'].startswith('RIF') else 'Control',run=s['run_accession'],md5=h.hexdigest(),bytes=p.stat().st_size))
 print('DOWNLOAD_COMPLETE',s['title'],flush=True)
elif mode=='index':
 ref=ROOT/'reference';idx=ref/'STAR_GRCh38_gencode49';idx.mkdir(exist_ok=True)
 fasta=SOURCE/'GRCh38.primary_assembly.genome.fa';fai=Path(str(fasta)+'.fai');gtf=ref/'gencode.v49.primary_contigs.gtf'
 allowed={line.split('\t')[0] for line in fai.open()};kept=excluded=0
 if not gtf.exists():
  with gzip.open(SOURCE/'gencode.v49.annotation.gtf.gz','rt') as src,gtf.open('w') as out:
   for line in src:
    if line.startswith('#') or line.split('\t',1)[0] in allowed:out.write(line);kept+=1
    else:excluded+=1
 else:kept=None;excluded=None
 cmd=['STAR','--runMode','genomeGenerate','--runThreadN','8','--genomeDir',str(idx),'--genomeFastaFiles',str(fasta),'--sjdbGTFfile',str(gtf),'--sjdbOverhang','149','--limitGenomeGenerateRAM','60000000000','--outFileNamePrefix',str(ROOT/'logs/index_')]
 if not (ref/'index_complete.json').exists():run(cmd,ROOT/'logs/STAR_index.log')
 assert (idx/'Genome').stat().st_size>1000000000 and (idx/'SA').stat().st_size>1000000000
 save(ref/'index_complete.json',dict(command=cmd,genome=str(fasta),genome_sha256=digest(fasta),annotation=str(gtf),annotation_sha256=digest(gtf),annotation_lines_retained=kept,annotation_lines_excluded=excluded,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 print('INDEX_COMPLETE',flush=True)
elif mode=='sample':
 s=samples[int(os.environ['SLURM_ARRAY_TASK_ID'])];runid=s['run_accession'];out=ROOT/'samples'/runid;out.mkdir(parents=True,exist_ok=True)
 assert (ROOT/'raw'/(runid+'.verified.json')).exists() and (ROOT/'reference/index_complete.json').exists()
 if (out/'complete.json').exists():print('ALREADY_COMPLETE',runid);sys.exit(0)
 raw=ROOT/'raw'/(runid+'.fastq.gz');trim=out/'trimmed.fastq.gz'
 # Fixed, outcome-blind trimming. Retain pre/post FastQC for inspection of adapters and read quality.
 fastqc=out/'fastqc';fastqc.mkdir(exist_ok=True)
 run(['fastqc','--threads','4','--outdir',str(fastqc),str(raw)],out/'fastqc_raw.log')
 trimcmd=['cutadapt','-j','4','-a','AGATCGGAAGAGCACACGTCTGAACTCCAGTCA','-q','20','-m','30','--json',str(out/'cutadapt.json'),'-o',str(trim),str(raw)]
 run(trimcmd,out/'cutadapt.log')
 run(['fastqc','--threads','4','--outdir',str(fastqc),str(trim)],out/'fastqc_trimmed.log')
 cmd=['STAR','--runThreadN','4','--genomeDir',str(ROOT/'reference/STAR_GRCh38_gencode49'),'--readFilesIn',str(trim),'--readFilesCommand','zcat','--outFileNamePrefix',str(out/'STAR_'),'--outSAMtype','BAM','Unsorted','--quantMode','GeneCounts','--outSAMattributes','NH','HI','AS','nM']
 run(cmd,out/'STAR_command.log')
 counts=out/'STAR_ReadsPerGene.out.tab';assert counts.exists() and counts.stat().st_size>100000
 save(out/'complete.json',dict(run=runid,sample=s['title'],condition='RIF' if s['title'].startswith('RIF') else 'Control',trim_command=trimcmd,alignment_command=cmd,counts_sha256=digest(counts),completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 print('SAMPLE_COMPLETE',runid,s['title'],flush=True)
elif mode=='collect':
 import csv
 out=ROOT/'results';out.mkdir(exist_ok=True);tables=[];qc=[]
 for s in samples:
  p=ROOT/'samples'/s['run_accession'];assert (p/'complete.json').exists()
  d=[line.rstrip().split('\t') for line in (p/'STAR_ReadsPerGene.out.tab').open()];tables.append(d)
  q=json.loads((p/'cutadapt.json').read_text());alignment={}
  for line in (p/'STAR_Log.final.out').read_text().splitlines():
   if '|' in line:
    k,v=line.split('|',1);alignment[k.strip()]=v.strip()
  qc.append(dict(run=s['run_accession'],sample=s['title'],condition='RIF' if s['title'].startswith('RIF') else 'Control',raw_reads=q['read_counts']['input'],raw_mean_length=q['basepair_counts']['input']/q['read_counts']['input'],retained_reads=q['read_counts']['output'],trimmed_mean_length=q['basepair_counts']['output']/q['read_counts']['output'],reads_with_adapter=q['read_counts']['read1_with_adapter'],cutadapt_version=q['cutadapt_version'],unique_alignment=alignment.get('Uniquely mapped reads %'),assigned_unstranded=sum(int(r[1]) for r in d[4:]),assigned_forward=sum(int(r[2]) for r in d[4:]),assigned_reverse=sum(int(r[3]) for r in d[4:])))
 genes=[r[0] for r in tables[0][4:]];assert all([r[0] for r in t[4:]]==genes for t in tables)
 for col,label in [(1,'unstranded'),(2,'forward'),(3,'reverse')]:
  with gzip.open(out/('Xu2026_'+label+'_counts.csv.gz'),'wt',newline='') as f:
   w=csv.writer(f);w.writerow(['EnsemblID']+[s['title'] for s in samples]);w.writerows([g]+[t[i+4][col] for t in tables] for i,g in enumerate(genes))
 with (out/'Xu2026_sample_qc.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=qc[0]);w.writeheader();w.writerows(qc)
 save(out/'complete.json',dict(samples=22,genes=len(genes),pipeline_sha256=digest(Path(__file__)),container_sha256=digest(Path('<cluster_project>/NXF_SINGULARITY_CACHEDIR/transcriptomics-core_0.2.1.sif')),completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 import tarfile
 bundle=ROOT/'Xu2026_results_and_qc.tar.gz'
 with tarfile.open(bundle,'w:gz') as tf:
  for folder in ['results','logs','scripts']:tf.add(ROOT/folder,arcname=folder)
  for name in ['index_complete.json','Xu2026_sample_metadata.json']:tf.add(ROOT/'reference'/name,arcname='reference/'+name)
  for s in samples:
   runid=s['run_accession'];p=ROOT/'samples'/runid
   tf.add(ROOT/'raw'/(runid+'.verified.json'),arcname='raw/'+runid+'.verified.json')
   for name in ['complete.json','cutadapt.json','cutadapt.log','STAR_Log.final.out','STAR_Log.out','STAR_ReadsPerGene.out.tab']:
    tf.add(p/name,arcname='samples/'+runid+'/'+name)
   for f in (p/'fastqc').glob('*.zip'):tf.add(f,arcname='samples/'+runid+'/fastqc/'+f.name)
 save(ROOT/'Xu2026_results_and_qc.manifest.json',dict(file=bundle.name,bytes=bundle.stat().st_size,sha256=digest(bundle)))
 print('XU_COUNTS_COMPLETE',len(genes),flush=True)
else:raise ValueError(mode)
