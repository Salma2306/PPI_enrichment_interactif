#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Single-script PPI + enrichment workflow v1.0.0.
Primary: five fixed shortlisted targets vs study-derived target universe.
Secondary: high-confidence STRING first-shell network context.
No re-ranking, no permutations, no forced significance.
"""
from __future__ import annotations
import argparse,csv,io,json,math,re,time
from pathlib import Path
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt

VERSION='1.0.0'; DEFAULT_SEEDS=['REN','SELP','ANGPT2','CXCL10','IFNG']
STRING='https://version-12-0.string-db.org'; GPROF='https://biit.cs.ut.ee/gprofiler/api/gost/profile/'; SPECIES=9606
SOURCES=['GO:BP','GO:MF','GO:CC','REAC','KEGG','WP']; DISPLAY=['GO:BP','REAC','KEGG','WP']

def ask(p,d,cast=str):
 s=input(f'{p} [{d}]: ').strip()
 if not s:return d
 try:return cast(s)
 except:print('Invalid value; default retained.');return d
def genes(x):return list(dict.fromkeys(z.upper() for z in re.split(r'[;,\s]+',str(x)) if z.strip()))
def req(ss,m,u,**kw):
 err=None
 for i in range(5):
  try:
   r=ss.request(m,u,timeout=120,**kw)
   if r.status_code==429 or r.status_code>=500:raise RuntimeError(f'HTTP {r.status_code}')
   r.raise_for_status();return r
  except Exception as e:err=e;time.sleep(min(2**i,16))
 raise RuntimeError(err)
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')
def locate_background(root,explicit=''):
 r=Path(root);cand=[]
 if explicit:cand.append(Path(explicit).expanduser())
 cand += [r/'06_cross_database/target_evidence_matrix.csv',r/'target_evidence_matrix.csv']
 if r.exists():cand += list(r.rglob('target_evidence_matrix.csv'))
 for p in cand:
  if p.is_file():
   d=pd.read_csv(p,dtype=str,keep_default_na=False,encoding_errors='replace');col=next((c for c in ['Target','HGNC_Symbol','Gene','Symbol'] if c in d.columns),None)
   if col:
    bg=list(dict.fromkeys(str(v).strip().upper() for v in d[col] if str(v).strip()))
    return p,bg
 raise FileNotFoundError('target_evidence_matrix.csv not found')
def smap(ss,ids):
 txt=req(ss,'POST',STRING+'/api/tsv/get_string_ids',data={'identifiers':'\r'.join(ids),'species':SPECIES,'limit':1,'caller_identity':'PPIEnrichmentSingleFinal'}).text
 return pd.DataFrame(list(csv.DictReader(txt.splitlines(),delimiter='\t')))
def snet(ss,ids,score):
 if not ids:return pd.DataFrame()
 txt=req(ss,'POST',STRING+'/api/tsv/network',data={'identifiers':'\r'.join(ids),'species':SPECIES,'required_score':score,'add_nodes':0,'caller_identity':'PPIEnrichmentSingleFinal'}).text.strip()
 return pd.read_csv(io.StringIO(txt),sep='\t') if txt else pd.DataFrame()
def spart(ss,sid,k,score):
 txt=req(ss,'POST',STRING+'/api/tsv/interaction_partners',data={'identifiers':sid,'species':SPECIES,'limit':k,'required_score':score,'caller_identity':'PPIEnrichmentSingleFinal'}).text.strip()
 return pd.read_csv(io.StringIO(txt),sep='\t') if txt else pd.DataFrame()
def enrich(ss,q,bg=None):
 body={'organism':'hsapiens','query':q,'sources':SOURCES,'user_threshold':0.05,'all_results':True,'ordered':False,'significance_threshold_method':'fdr','domain_scope':'custom_annotated' if bg else 'annotated'}
 if bg:body['background']=bg
 d=pd.DataFrame(req(ss,'POST',GPROF,json=body).json().get('result',[]))
 if not d.empty:
  keep=[c for c in ['source','native','name','p_value','significant','term_size','query_size','intersection_size','effective_domain_size','precision','recall','intersection'] if c in d.columns];d=d[keep].copy()
  if 'intersection' not in d:d['intersection']=''
  else:d['intersection']=d.intersection.map(lambda x:'; '.join(x) if isinstance(x,list) else str(x))
  d=d.sort_values('p_value')
 return d,body
def graph(d,nodes):
 G=nx.Graph();G.add_nodes_from(nodes)
 if not d.empty:
  for _,r in d.iterrows():
   a,b=str(r.get('preferredName_A','')),str(r.get('preferredName_B',''))
   if a and b:G.add_edge(a,b,score=float(r.get('score',0) or 0))
 return G
def draw_network(ax,G,seeds,title):
 pos=nx.spring_layout(G,seed=42,weight='score',k=max(.35,1/math.sqrt(max(len(G),1))))
 shared={n:sum(G.has_edge(n,s) for s in seeds if s in G) for n in G if n not in seeds}
 nx.draw_networkx_edges(G,pos,ax=ax,alpha=.28,edge_color='#8FA0AC',width=[.6+2*G[u][v].get('score',0) for u,v in G.edges])
 other=[n for n in G if n not in seeds];colors=['#159A96' if shared.get(n,0)>=2 else '#C6D2DB' for n in other]
 if other:nx.draw_networkx_nodes(G,pos,nodelist=other,node_color=colors,node_size=250,edgecolors='white',linewidths=.5,ax=ax)
 sd=[n for n in G if n in seeds];nx.draw_networkx_nodes(G,pos,nodelist=sd,node_color='#D95D58',node_size=720,edgecolors='#7A2B28',linewidths=1,ax=ax)
 labels={n:n for n in sd};labels.update({n:n for n in other if shared.get(n,0)>=2});nx.draw_networkx_labels(G,pos,labels=labels,font_size=7,font_weight='bold',ax=ax);ax.set_title(title,loc='left',fontweight='bold');ax.axis('off')
def display_terms(d,n=10):
 if d.empty:return d
 x=d[d.source.isin(DISPLAY)&d.significant.astype(bool)].copy()
 if x.empty:return x
 return pd.concat([g.nsmallest(3,'p_value') for _,g in x.groupby('source')]).sort_values('p_value').head(n)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--targets');ap.add_argument('--target-output-root');ap.add_argument('--background');ap.add_argument('--out');ap.add_argument('--string-score',type=int,default=700);ap.add_argument('--partners-per-seed',type=int,default=5);x=ap.parse_args()
 seeds=genes(x.targets or ask('Shortlisted targets','; '.join(DEFAULT_SEEDS)));root=x.target_output_root or ask('Target-discovery output root',r'D:\Coag_gen_Literature_2020_2026\target_output');out=Path(x.out or str(Path(root)/'18_PPI_enrichment_final'));out.mkdir(parents=True,exist_ok=True)
 bgfile,bg=locate_background(root,x.background or '');print(f'Background: {bgfile} | {len(bg)} targets');
 if len(bg)!=197:print(f'WARNING: expected 197 study-derived targets; loaded {len(bg)}. No automatic correction.')
 ss=requests.Session();mp=smap(ss,seeds);mapped=list(mp.preferredName.astype(str));ids=list(mp.stringId.astype(str));mp.to_csv(out/'01_STRING_mapping.csv',index=False)
 sedges=snet(ss,ids,x.string_score);sedges.to_csv(out/'02_seed_PPI_edges.csv',index=False);Gs=graph(sedges,mapped)
 seed_enr,sreq=enrich(ss,mapped,bg);seed_enr.to_csv(out/'03_seed_enrichment_vs_197.csv',index=False);sig=seed_enr[seed_enr.significant.astype(bool)].copy() if not seed_enr.empty else seed_enr;sig.to_csv(out/'04_seed_enrichment_FDR_significant.csv',index=False)
 first=[];expanded=set(mapped)
 for _,r in mp.iterrows():
  d=spart(ss,str(r.stringId),x.partners_per_seed,x.string_score);seed=str(r.preferredName)
  for _,z in d.iterrows():
   a,b=str(z.get('preferredName_A','')),str(z.get('preferredName_B',''));p=b if a==seed else a
   if p:expanded.add(p);first.append({'Seed':seed,'Partner':p,'STRING_score':z.get('score','')})
 pd.DataFrame(first).to_csv(out/'05_first_shell_context.csv',index=False)
 em=smap(ss,sorted(expanded));eedges=snet(ss,list(em.stringId.astype(str)),x.string_score);eedges.to_csv(out/'06_expanded_PPI_edges.csv',index=False);Ge=graph(eedges,list(em.preferredName.astype(str)))
 ext,ereq=enrich(ss,list(em.preferredName.astype(str)),None);ext.to_csv(out/'07_expanded_enrichment_context_only.csv',index=False);pd.DataFrame({'Target':bg}).to_csv(out/'08_background_universe.csv',index=False)
 # Publication figure: direct PPI, primary ORA, expanded network context.
 fig,axs=plt.subplots(1,3,figsize=(18,6.8),gridspec_kw={'width_ratios':[1,1,1.4]});draw_network(axs[0],Gs,set(mapped),'A  Direct PPI of five shortlisted targets')
 t=display_terms(seed_enr)
 if t.empty:
  axs[1].axis('off');axs[1].set_title('B  Primary functional enrichment',loc='left',fontweight='bold');axs[1].text(.5,.64,f'{len(mapped)} shortlisted targets',ha='center',fontsize=14,fontweight='bold');axs[1].text(.5,.52,'versus',ha='center',fontsize=9,color='#687783');axs[1].text(.5,.42,f'{len(bg)}-target study universe',ha='center',fontsize=12,fontweight='bold',color='#173B5E');axs[1].text(.5,.26,f'FDR-significant terms: {len(sig)}',ha='center',fontsize=14,fontweight='bold',color='#D95D58' if len(sig)==0 else '#159A96')
 else:
  t=t.iloc[::-1];v=[-math.log10(max(float(p),1e-300)) for p in t.p_value];colors={'GO:BP':'#397CB7','REAC':'#D95D58','KEGG':'#E5A338','WP':'#8A6BBE'};axs[1].barh(range(len(t)),v,color=[colors.get(s,'#78909C') for s in t.source]);axs[1].set_yticks(range(len(t)));axs[1].set_yticklabels([str(z)[:50] for z in t.name],fontsize=7);axs[1].set_xlabel('-log10(FDR)');axs[1].set_title('B  Primary functional enrichment',loc='left',fontweight='bold');axs[1].spines[['top','right']].set_visible(False)
 draw_network(axs[2],Ge,set(mapped),'C  High-confidence first-shell PPI context');fig.suptitle('Network context of shortlisted pulmonary immunothrombosis targets',fontsize=16,fontweight='bold');fig.text(.02,.015,'Primary ORA compares the fixed five-target shortlist with the fixed study-derived target universe. STRING expansion provides descriptive context only and does not alter target selection.',fontsize=8,color='#586875');fig.tight_layout(rect=[0,.04,1,.95])
 for e in ['png','pdf','svg']:fig.savefig(out/f'Figure_3_PPI_enrichment.{e}',dpi=600 if e=='png' else None,bbox_inches='tight',facecolor='white')
 plt.close(fig)
 qc={'version':VERSION,'seeds_requested':seeds,'seeds_mapped':mapped,'background_file':str(bgfile),'background_n':len(bg),'STRING_score':x.string_score,'partners_per_seed':x.partners_per_seed,'seed_FDR_significant_terms':int(len(sig)),'seed_PPI_edges':Gs.number_of_edges(),'expanded_PPI_nodes':Ge.number_of_nodes(),'expanded_PPI_edges':Ge.number_of_edges(),'guards':['Five-target shortlist fixed before analysis.','Study-derived universe fixed as primary enrichment background.','Negative seed ORA retained.','STRING expansion contextual only.','No permutations, leave-one-out, target re-ranking, adaptive thresholding, or forced significance.']};dump(out/'09_QC_provenance.json',qc);dump(out/'10_gProfiler_requests.json',{'seed':sreq,'expanded_context':ereq});print('\nCOMPLETE:',out);print(json.dumps(qc,indent=2))
if __name__=='__main__':main()
