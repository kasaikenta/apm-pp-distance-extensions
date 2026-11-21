"""Verify released matrix/witness bindings and optional archive checksums."""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text())
def sha_file(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
 return h.hexdigest()
def basis(rows):
 b={}
 for v in rows:
  while v:
   i=v.bit_length()-1
   if i not in b:b[i]=v;break
   v^=b[i]
 return b
def rem(v,b):
 while v:
  i=v.bit_length()-1
  if i not in b:return v
  v^=b[i]
 return 0
def bits(r,n):
 assert len(r)==len(set(r)) and all(type(x)==int and 0<=x<n for x in r)
 return sum(1<<x for x in r)
p=argparse.ArgumentParser();p.add_argument('--assets',type=Path);a=p.parse_args()
for entry in read(ROOT/'exact-distance/index.json'):
 path=ROOT/entry['matrix'];assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['matrix_sha256'];c=read(path);n=c['n'];X,Z=[[bits(r,n) for r in c[s]] for s in ['hx','hz']];BX,BZ=basis(X),basis(Z)
 assert n==entry['n'] and n-len(BX)-len(BZ)==entry['k']
 assert all((x&z).bit_count()%2==0 for x in X for z in Z)
 weights=[]
 for side in ['X','Z']:
  ws=entry.get('physical_witnesses',c['witnesses'])[side];ws=[ws['support']] if isinstance(ws,dict) else ws
  for support in ws:
   v=bits(support,n);assert v and all((v&r).bit_count()%2==0 for r in (Z if side=='X' else X));assert rem(v,BX if side=='X' else BZ);weights.append(v.bit_count())
 assert min(weights)==entry['d'];print(entry['id'],'physical PASS')
if a.assets:
 for item in read(ROOT/'exact-distance/assets.json'):
  f=a.assets/item['name'];assert f.stat().st_size==item['bytes'];assert sha_file(f)==item['sha256'];print(item['name'],'SHA256 PASS')
print('Matrix identities and physical upper witnesses verified; lower exclusions require audit_release.py.')
