"""Audit archived exact-distance evidence. Does not launch a distance search."""
import argparse,hashlib,json,sys,tempfile,sqlite3,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--case',default='standard',choices=['standard','p192','p384','cpm10-P73','cpm10-P37','all']);a=p.parse_args();E=a.evidence.resolve();sys.path.insert(0,str(E/'vendor/QEDist'));sys.path.insert(0,str(E))
from cssdist.code import CSSCode,digest
from cssdist.algebra import basis,remainder,vector,permute,same_span
from cssdist.audit import audit_run,check_receipt,configuration
from cssdist.lp_blocks import search_problems
from audit_adaptive import header

def read(p):return json.loads(p.read_text())
def same_physical(public,search):
 for field in ('hx','hz'):
  aa=public[field];bb=search[field]
  if aa!=bb:
   P=public['n']//8;T=search['cyclic_period'];g=P//T
   perm=[(q%8*g+(q//8)%g)*T+(q//8)//g for q in range(8*P)]
   assert sorted(perm)==list(range(public['n']))
   assert [sorted(perm[q] for q in r) for r in aa]==bb

def standard():
 for r in read(E/'standard-campaigns.json'):
  public=read(E/'public'/r['public_matrix']);lo=1;hi=public['n'];cnt=0
  for path in r['campaigns']:
   root=E/path;doc=read(root/'code.json');same_physical(public,doc)
   z=audit_run(root);lo=max(lo,z['distance_lower_bound']);hi=min(hi,z['distance_upper_bound']);cnt+=sum(x['attempts'] for x in z['jobs'])
  assert lo==hi==r['distance'],(r['id'],lo,hi)
  print(json.dumps(dict(id=r['id'],d=lo,receipts=cnt,status='PASS')),flush=True)

def adaptive(name):
 same_physical(read(E/f'public/data/codes/apm-P{name[1:]}.json'),read(E/name/'code.json'))
 with tempfile.TemporaryDirectory() as tmp:
  out=Path(tmp)/'audit.json'
  subprocess.run([sys.executable,str(E/name/'audit.py'),'--output',str(out)],check=True)
  got=read(out);old=read(E/name/'verified.json')
  assert got['exact']==got['lower']==got['upper']==old['exact']
  assert got['attempts']==old['attempts'] and got['coverage']==old['coverage']
 print(json.dumps(dict(id=name,d=got['exact'],status='PASS')),flush=True)

def cpm73():
 root=E/'cpm10-P73';code=CSSCode.load(root/'c00/code.json');base,_=header(root/'c00',code)
 public=read(E/'public/cpm10/P73.json');assert all(public[k]==code.document[k] for k in ('n','hx','hz'))
 assert code.exchange is not None and all(code.kernel_even_weight(s) for s in ['X','Z'])
 lower=1;upper=code.n;count=0
 for w in [22,24]:
  cert=read(root/f'certificate_w{w}.json');cfg=configuration(code,w,2000,8,True,True,'css',True,True);problem=search_problems(code,cfg)['Z'];good={};terms={};witnesses=[]
  with tempfile.TemporaryDirectory() as t:
   t=Path(t)
   for f in cert['record_files']:
    raw=(root/f['file']).read_bytes();assert digest(raw)==f['sha256']
    for line in raw.splitlines():
     r=json.loads(line);mod=r['modulus'];expected=dict(cfg,shards=mod);assert r['candidate']==0 and r['config']==expected and mod in [2000,20000,200000]
     rec=r['receipt'];s=rec['shard'];assert 0<=s<mod
     for n,data in [('receipt.json',json.dumps(rec)),('stdout.jsonl',r['stdout']),('stderr.txt',r['stderr'])]:(t/n).write_text(data)
     campaign=dict(base,config=expected);checked=check_receipt(t/'receipt.json',t,campaign,code,'Z',s,problem);count+=1
     if checked['status']=='excluded':
      good.setdefault(mod,set()).add(s);terms[(mod,s)]=checked['terminal']
     witnesses+=checked['witnesses']
  cover={m:set(good.get(m,set())) for m in [2000,20000,200000]}
  for m in [200000,20000]:
   parent=m//10
   cover[parent].update(s for s in {x%parent for x in cover[m]} if all(s+parent*j in cover[m] for j in range(10)))
  assert sorted(cover[2000])==cert['excluded']
  if len(cover[2000])==2000:
   assert not witnesses
   signatures={(t['frontier_seen'],t['frontier_digest'],t['common_nodes']) for t in terms.values()};assert len(signatures)==1
   def owned(mod,shard):
    if (mod,shard) in terms:return terms[(mod,shard)]['frontier_owned']
    assert mod<200000
    return sum(owned(mod*10,shard+mod*j) for j in range(10))
   assert sum(owned(2000,s) for s in range(2000))==next(iter(signatures))[0]
   lower=max(lower,w+2)
  for word in witnesses:assert code.is_logical(word,'Z');upper=min(upper,word.bit_count())
 assert lower==upper==24
 print(json.dumps(dict(id='cpm10-P73',d=24,receipts=count,status='PASS')),flush=True)

def cpm37():
 root=E/'cpm10-P37';meta=read(root/'source_distance_audit.json');eq=read(root/'source_equivalence.json')
 assert digest((root/'source_distance_audit.json').read_bytes())==eq['source_distance_audit_sha256']
 for sec in ['lower_bound','upper_bound','duality']:assert digest((root/f'{sec}.json').read_bytes())==meta[sec]['certificate_sha256']
 assert digest((root/'source_P37_d22.cpm').read_bytes())==meta['candidate_sha256']
 cov=read(root/'lower_bound.json');expected={(r['root'],r['shard_index']):r['sha256'] for r in cov['certificates']};assert len(expected)==len(cov['certificates'])==5005
 observed={};counts={}
 for f in (root/'units').glob('*.result.json'):
  r=read(f);key=(r['root'],r['shard_index']);assert key not in observed
  assert r['complete'] and not r['found'] and r['n']==370 and r['P']==37 and r['side']=='X' and r['max_weight']==20
  assert digest(f.read_bytes())==expected[key];observed[key]=r
  counts.setdefault(key[0],r['shard_count']);assert counts[key[0]]==r['shard_count']
 assert set(observed)==set(expected)=={(rt,s) for rt,c in counts.items() for s in range(c)} and set(counts)==set(range(10))
 lines=(root/'source_P37_d22.cpm').read_text().splitlines();ex=[list(map(int,l.split())) for l in lines[2:10]];P=37
 mats=[[[l*P+(v-ex[4*side+i][l])%P for l in range(10)] for i in range(4) for v in range(P)] for side in range(2)]
 X,Z=[[vector(r,370) for r in H] for H in mats];BX,BZ=basis(X),basis(Z)
 assert 370-len(BX)-len(BZ)==80 and all(not((x&z).bit_count()%2) for x in X for z in Z)
 assert all(not remainder((1<<370)-1,b) for b in [BX,BZ])
 dm=read(root/'duality.json')['side_swapping_maps'][0]
 perm=[dm['columns'][l]*P+(dm['multiplier']*u+dm['qubit_offsets'][l])%P for l in range(10) for u in range(P)]
 assert sorted(perm)==list(range(370));assert same_span([permute(x,perm) for x in X],Z) and same_span([permute(z,perm) for z in Z],X)
 transformed=read(root/'input_P37.json');p=eq['physical_permutation'];assert sorted(p)==list(range(370))
 for H,key in [(X,'hx'),(Z,'hz')]:assert same_span([permute(h,p) for h in H],[vector(r,370) for r in transformed[key]])
 code=CSSCode(dict(transformed,schema='cssdist.code.v1'));v=vector(eq['witness_support'],370);assert v.bit_count()==22 and code.is_logical(v,'X')
 print(json.dumps(dict(id='cpm10-P37',d=22,receipts=5005,status='PASS',scope='Archived PP-verifier certificates, coverage, physical parity, side exchange and seed equivalence; no per-node proof trace')),flush=True)
for name,fn in [('standard',standard),('p192',lambda:adaptive('p192')),('p384',lambda:adaptive('p384')),('cpm10-P73',cpm73),('cpm10-P37',cpm37)]:
 if a.case in [name,'all']:fn()
