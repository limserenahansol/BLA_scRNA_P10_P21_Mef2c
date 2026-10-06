// Uses the supplied Codex Artifact Tool runtime; the R/Python analysis is independent.
import fs from 'node:fs/promises';
import path from 'node:path';
import {pathToFileURL, fileURLToPath} from 'node:url';
import {Presentation, PresentationFile, FileBlob} from '@oai/artifact-tool';
const ROOT = path.resolve(process.env.BLA_V2_ROOT || 'analysis_development_v2');
const SKILL_DIR = process.env.BLA_PRESENTATION_SKILL || 'C:/Users/hsollim/.codex/plugins/cache/openai-primary-runtime/presentations/26.904.11930/skills/presentations';
const RUNTIME_PYTHON = process.env.BLA_RUNTIME_PYTHON || 'C:/Users/hsollim/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe';
process.env.RUNTIME_NODE_MODULES ||= path.resolve(path.dirname(fileURLToPath(import.meta.resolve('@oai/artifact-tool'))),'../../..');
process.env.RUNTIME_NODE ||= process.execPath;
const {resolvePresentationFont, applyPresentationChartFont, finalizePresentation} = await import(pathToFileURL(path.join(SKILL_DIR,'container_tools/artifact_tool_utils.mjs')).href);
const family = resolvePresentationFont({fontFamily:'Arial'});
const data = JSON.parse(await fs.readFile(path.join(ROOT,'results/tables/presentation_data.json'),'utf8'));
const p = Presentation.create({slideSize:{width:1280,height:720}});
const ink='#173449', teal='#087E8B', red='#C65050', muted='#536471';
const tableOwners=[], chartOwners=[];
const nfmt = x => Number(x).toLocaleString('en-US');
function text(s, value, x,y,w,h,size=24,color=ink,bold=false) {
  const sh=s.shapes.add({geometry:'textbox',position:{left:x,top:y,width:w,height:h},fill:'none',line:{fill:'none',width:0}});
  sh.text=String(value); sh.text.style={typeface:family,fontSize:size,color,bold,autoFit:'none'}; return sh;
}
function slide(title,subtitle='',notes='') {
  const s=p.slides.add(); s.background.fill='#FFFFFF';
  text(s,title,60,32,1160,95,42,ink,true);
  if(subtitle) text(s,subtitle,62,122,1156,56,23,muted);
  text(s,'Hansol Lim • Developmental amygdala scRNA-seq • 06 Oct 2026',62,684,1000,24,16,muted);
  text(s,String(p.slides.items.length),1160,684,60,24,16,muted);
  s.speakerNotes.textFrame.setText(notes); return s;
}
function note(s,v,y=605) {text(s,v,64,y,1152,64,23,muted);}
function table(s,values,{x=64,y=188,w=1152,h=380,font=23,widths}={}) {
  const t=s.tables.add({rows:values.length,columns:values[0].length,left:x,top:y,width:w,height:h,values,
    ...(widths?{columnTracks:widths.map(value=>({mode:'fr',value}))}:{})});
  t.cells.block({row:0,column:0,rowCount:values.length,columnCount:values[0].length}).assign({fill:'#FFFFFF',textStyle:{typeface:family,fontSize:font,color:ink},margins:{left:10,right:10,top:2,bottom:2}});
  t.cells.block({row:0,column:0,rowCount:1,columnCount:values[0].length}).assign({fill:ink,textStyle:{typeface:family,fontSize:font,bold:true,color:'#FFFFFF'}});
  for(let r=1;r<values.length;r++) if(r%2===0) t.cells.block({row:r,column:0,rowCount:1,columnCount:values[0].length}).fill='#F1F5F7';
  tableOwners.push(p.slides.items.length); return t;
}
function chart(s,type,config) {
  // Display data round to five decimals for portable Excel workbooks; full
  // numerical precision is retained in the analysis CSV/JSON source tables.
  for (const series of config.series ?? []) for (const key of ['values','xValues'])
    if (series[key]) series[key]=series[key].map(v=>Number(v.toFixed(5)));
  const c=s.charts.add(type,{chartFill:'#FFFFFF',plotAreaFill:'#FFFFFF',...config});
  applyPresentationChartFont(c,{fontFamily:family}); chartOwners.push(p.slides.items.length); return c;
}
const axis=(title,more={})=>({title:{text:title,textStyle:{typeface:family,fontSize:20,fill:ink}},textStyle:{typeface:family,fontSize:18,fill:ink},...more});
async function figure(s,name,{x=62,y=184,w=1156,h=402}={}) {
  s.images.add({blob:new Uint8Array(await fs.readFile(path.join(ROOT,'results/figures',name+'.png'))),contentType:'image/png',alt:name,fit:'contain',position:{left:x,top:y,width:w,height:h}});
}

let s=slide('Developmental amygdala scRNA-seq','E18 · P0 · P10 · P21 — an exploratory update for Wei-Ming Kao',
  'Source: input_manifest.csv, sample_qc_summary.csv, doublet_summary.csv. No labelled KO/WT data, animal IDs, or connectivity measurements. All experimental units are libraries, not cells.');
text(s,'What changes when the early libraries are added?',64,210,1100,55,34,ink,true);
text(s,`${nfmt(data.putative_singlets)} putative singlets\n${data.libraries} libraries across four stages\n${nfmt(data.common_genes)} shared genes`,64,292,570,230,34,teal,true);
text(s,'Population inventory\nP10 transcription factors\nGuidance / remodelling candidates',690,292,530,220,30,ink);
note(s,'This dataset can generate testable candidates. It cannot explain a MEF2C-KO projection phenotype by itself.');

s=slide('What is in the combined dataset?','Seven new libraries + three previous P10/P21 libraries',
  'Sources: sample_qc_summary.csv, doublet_summary.csv. New samples: E18_1 P188 SI-GA-F3; E18_2 P189 SI-GA-F4; P0_1/P0_2 P165 SI-GA-E6/E9; P0_3 P188 SI-GA-F2; P10_1/P10_2 P173 SI-GA-E10/E11. Original new object 32285 genes x156727 barcodes; counts-only. P56 old samples are excluded. Animal independence and reference build are undocumented.');
table(s,[['Library','Stage','Source','Input','QC-pass','Putative singlets'],...data.qc.map(r=>[r.sample,r.age,r.source==='new_2026'?'New':'Previous',nfmt(r.input_barcodes),nfmt(r.qc_primary),nfmt(data.doublets.find(z=>z.sample===r.sample).singlets)])],{y:180,h:428,font:21,widths:[1.3,.7,1,1,1,1.4]});
note(s,'Shared-gene analysis avoids treating the 7,256 new-only genes as absent in old samples.',613);

s=slide('QC is a major source of uncertainty','Simple QC retention and classifier-labelled doublets — separate measures',
 'Source: QC summaries. Main QC >=500 genes, >=1000 UMIs, mito<=20%. scDblFinder per library; dbr=.08,dbr.sd=1. Doublets are putative. P0 low-count barcodes cannot be proven empty droplets from this supplied object alone. https://www.bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html');
chart(s,'bar',{position:{left:62,top:188,width:1152,height:380},categories:data.qc.map(r=>r.sample),
 series:[{name:'QC-pass / input (%)',values:data.qc.map(r=>100*r.qc_primary/r.input_barcodes),fill:teal},{name:'Putative doublets / QC (%)',values:data.qc.map(r=>data.doublets.find(z=>z.sample===r.sample).predicted_doublet_pct),fill:red}],
 barOptions:{direction:'column',grouping:'clustered'},hasLegend:true,legend:{position:'bottom',textStyle:{fontSize:20}},xAxis:axis('Library'),yAxis:axis('Percent (%)',{min:0,max:100,majorUnit:20}),dataLabels:{showValue:false}});
note(s,'P0_2: only 8.0% pass QC. P10_2: 32.1% putative doublets. Old/new mitochondrial capture differs sharply.');

s=slide('A broad map, not a cell-type verdict','Both uncorrected and source-corrected UMAPs are preserved',
 'Source: V2_F02 and cell_metadata.csv. Embedding only: shared-gene RNA log normalization, HVG3000 excluding candidates/ribo/mito, PCA30, Harmony source theta1/lambda1, UMAP30/min_dist.3/cosine/seed1062026. Neither map is ground truth. RNA integrity hashes unchanged. Source correction does not fully mix old/new P10.');
await figure(s,'V2_F02_UMAP_baseline_corrected');
note(s,'Labels use excluded markers and reference checks. Ambiguous cells remain visible; fine developmental identities are provisional.');

s=slide('Stage and source are partly confounded','P10 overlaps sources; P21 exists only in the previous source',
 'Source: V2_F03, embedding_validation.csv. P10 balanced-source neighbor mixing rises but remains far below balanced random mixing. Do not interpret a source-specific island as maturation or a new subtype.');
await figure(s,'V2_F03_UMAP_age_source');
note(s,'Harmony is for visualization only. All expression calculations use unchanged, uncorrected RNA counts.');

s=slide('Captured populations differ by library','Fractions of retained barcodes — not proportions in intact tissue',
 'Source: population_composition_per_library.csv, primary_singlets. Minor classes grouped only in this chart; full CSV retains all classes. ITC-like and Rspo2-positive subsets are nested and not added to broad denominators. No independent animal n.');
const major=['Glutamatergic','GABAergic','Immature_neuron','Astrocyte','Progenitor','OPC','Oligodendrocyte','Microglia','Ambiguous'];
const colors=['#D97706','#5367B4','#A259A6','#00A896','#81B29A','#70A0D5','#35638C','#A05A2C','#C3C7CC','#AA7B8B'];
chart(s,'bar',{position:{left:60,top:182,width:1156,height:405},categories:data.qc.map(r=>r.sample),
 series:[...major,'Other non-neuronal'].map((pop,i)=>({name:pop.replace('_',' '),fill:colors[i],values:data.qc.map(q=>data.composition.filter(z=>z.sample===q.sample&&(pop==='Other non-neuronal'?!major.includes(z.population):z.population===pop)).reduce((a,z)=>a+z.percent,0))})),
 barOptions:{direction:'column',grouping:'stacked',overlap:100},hasLegend:true,legend:{position:'bottom',textStyle:{fontSize:17}},xAxis:axis('Library'),yAxis:axis('Retained barcodes (%)',{min:0,max:100,majorUnit:25}),dataLabels:{showValue:false}});
note(s,'Library capture, depth, ambiguous labels and doublet selection prevent a definitive developmental cell census.');

s=slide('Annotation checks expose unresolved cells','The adult-atlas-only baseline did not validate well',
 'Source: annotation_validation_per_sample.csv computed before doublet exclusion. Training1457 and heldout1457 reference genes exclude marker/candidate genes. Previous-reference predictions inherit v1 uncertainty. Leave-one-old-library-out label agreement84.6/90.5/90.5% is consistency, not independent truth. Adult reference Yao2023: https://www.nature.com/articles/s41586-023-06812-z');
table(s,[['Library','Adult atlas / marker','Previous broad split','Ambiguous (QC %)'],...data.annotation.map(r=>[r.sample,(100*r.independent_marker_atlas_agreement).toFixed(1)+'%',(100*r.previous_broad_split_agreement).toFixed(1)+'%',r.ambiguous_pct.toFixed(1)+'%'])],{y:181,h:392,font:21,widths:[1,1.8,1.8,1.2]});
note(s,'Broad marker-supported assignments are useful candidates. Novel or low-depth identities require further validation.',612);

s=slide('Mef2c expression: library-level evidence','Each dot is one library; old and new P10 are shown separately',
 'Source: sample_population_expression.csv.gz, primary singlets, n>=30 per group. CPM uses total shared-gene pseudobulk UMIs. No per-cell p-values or KO inference. MEF2C activity/binding is not measured.');
const pops=['Glutamatergic','GABAergic','ITC-like GABA (provisional)'];
const all=data.mef2c.filter(r=>pops.includes(r.population));
const ymax=Math.ceil(Math.max(...all.map(r=>Math.log2(r.cpm+.5)))+1);
for(let i=0;i<pops.length;i++) {
 const rows=all.filter(r=>r.population===pops[i]);
 text(s,pops[i].replace(' (provisional)','*'),65+i*393,186,375,42,25,ink,true);
 chart(s,'scatter',{position:{left:62+i*393,top:232,width:375,height:326},hasLegend:i===2,legend:{position:'bottom',textStyle:{fontSize:18}},
  series:['new_2026','old_2023'].map(src=>{const z=rows.filter(r=>r.source===src);return{name:src==='new_2026'?'New source':'Previous source',xValues:z.map(r=>data.stages.indexOf(r.age)+1),values:z.map(r=>Math.log2(r.cpm+.5)),fill:src==='new_2026'?teal:red,line:{fill:'none',width:0},marker:{symbol:src==='new_2026'?'circle':'square',size:8}};}),
  scatterOptions:{style:'marker'},xAxis:axis('1 E18 · 2 P0 · 3 P10 · 4 P21',{min:1,max:4,majorUnit:1}),yAxis:axis('log2(CPM + 0.5)',{min:-1,max:ymax,majorUnit:2})});
}
note(s,'*ITC-like identity is provisional and not spatially verified. Similar expression does not imply similar MEF2C function.');

s=slide('P10 transcription-factor expression','A candidate screen, not TF activity or a causal target ranking',
 'Source: TF_expression_per_library.csv.gz, supplied mouse_TF_list_1321.csv. Figure shows prespecified developmental/lineage TFs. Color is log2(mean library CPM+.5), size is expected detection at1000UMIs. Marker TFs such as Foxp2/Etv1 are annotation anchors.');
await figure(s,'V2_F09b_P10_TFs_slide',{x:64,y:178,w:1152,h:426});
note(s,'Use cell identity + candidate availability to plan experiments; verify protein, location and genotype-specific regulation.');

s=slide('Guidance and adhesion candidates across stages','Sema3e · Plxnd1/Nrp1 · Kirrel3 · cadherins · Slit/Robo',
 'Source: V2_F06 and candidate_expression_per_library.csv. Equal library weights, groups>=30 cells, common genes only. Rspo2-positive subset selected by Rspo2 detection; ITC-like tentative. Transcript detection is not receptor stoichiometry or attraction/repulsion. Chauvet2007 primary: https://doi.org/10.1016/j.neuron.2007.10.019');
await figure(s,'V2_F06_guidance',{x:60,y:177,w:1160,h:426});
note(s,'Expression nominates components. Functional cue response and the projecting cortical cell must be tested independently.');

s=slide('Synapse-remodelling and glial candidates','MEF2-associated activity genes and engulfment/complement machinery',
 'Source: V2_F07, raw pseudobulk CPM. C1qa/C1qb/Cx3cr1/P2ry12 are also annotation markers: expression in microglia is not independent cell-type enrichment proof. Flavell2006 activity-dependent MEF2 regulation in hippocampal synapses https://pubmed.ncbi.nlm.nih.gov/16484497/. No axonal pruning measurement here.');
await figure(s,'V2_F07_pruning',{x:60,y:177,w:1160,h:426});
note(s,'These observations do not distinguish retained axon branches, additional synapses, altered tracer labeling or changed cortical cell identity.');

s=slide('Cortical context is external adult reference','The new experiment samples amygdala, not the projecting cortex',
 'Source: prior results/tables/C01_cortex_glut_receptor_pct_adult_Allen.csv. Adult Allen atlas, not P10 or KO cortex. Detection does not indicate co-localized receptor complexes. Martin2015 Kirrel3 hippocampal targeting https://elifesciences.org/articles/09395.');
// A small prespecified editable table, sourced directly from the previous CSV.
const cortexCSV=await fs.readFile(path.join(ROOT,'../results/tables/C01_cortex_glut_receptor_pct_adult_Allen.csv'),'utf8');
const lines=cortexCSV.trim().split(/\r?\n/).map(l=>l.split(','));
const cortexGenes=['Kirrel3','Plxnd1','Nrp1','Nrp2','Sdk2','Robo1'];
table(s,[['Gene','L2/3 IT\nDetection %','L4/5 IT\nDetection %','L5 ET\nDetection %'],...cortexGenes.map(g=>{const r=lines.find(z=>z[0]===g);return[g,r?r[2]:'not available',r?r[3]:'not available',r?r[6]:'not available'];})],{h:340,font:25,widths:[1,1.3,1.3,1.3]});
note(s,'Do not claim Plxnd1 is exclusive to one cortical class, or that Nrp1-negative RNA means a repulsive receptor state.',582);

s=slide('P10–P21 directions must survive source checks','True pseudobulk CPM ratios; mito and putative-doublet sensitivity',
 'Source: candidate_P10_P21_descriptive_effects.csv and candidate_source_QC_sensitivity.csv. Gray means mixed directions or incomplete support; teal requires agreement across old/new sources,20/15% mito,inclusion/exclusion of putative doublets,and every library pair. No cell-level significance. Prior mean-log differences are not fold changes.');
await figure(s,'V2_F10_source_sensitivity',{x:60,y:177,w:1160,h:422});
note(s,'A consistent direction is stronger descriptive evidence, but P21-only-old-source confounding and unknown animal n remain.');

s=slide('Core hypotheses: consistent and unstable patterns','Prespecified genes shown whether they pass or fail the sensitivity rule',
 'Source: candidate_P10_P21_descriptive_effects.csv and candidate_source_QC_sensitivity.csv. Log2(CPM ratio) P10/P21, separately old/new P10. Consistency means all source/mito/doublet/library-pair directions agree, not significance or causality. All candidates, including failures, remain in CSVs.');
const core=[['Mef2c','Glutamatergic'],['Mef2c','GABAergic'],['Pcdh10','Glutamatergic'],['Kirrel3','ITC-like GABA (provisional)'],['Sema3e','Glutamatergic'],['Plxnd1','ITC-like GABA (provisional)'],['Nrp1','Glutamatergic']];
table(s,[['Gene / population','Old P10 / P21\nlog2 ratio','New P10 / P21\nlog2 ratio','All direction checks'],...core.map(([gene,pop])=>{
 const get=src=>data.core_effects.find(r=>r.gene===gene&&r.population===pop&&r.P10_source===src);
 const rr=data.sensitivity.find(r=>r.gene===gene&&r.population===pop);
 return[`${gene} — ${pop.replace('Glutamatergic','Glut').replace('GABAergic','GABA').replace(' (provisional)','*')}`,get('old_2023')?.log2_CPM_ratio_P10_over_P21.toFixed(2)??'n/a',get('new_2026')?.log2_CPM_ratio_P10_over_P21.toFixed(2)??'n/a',rr?.all_P10_P21_library_pair_directions_agree?'Consistent direction':'Mixed / unstable'];
})],{h:372,font:22,widths:[2.2,1.2,1.2,1.6]});
note(s,'Pcdh10 is direction-consistent. Kirrel3 ITC-like change reverses by source. Neither observation is a KO result.');

s=slide('Three scenarios remain plausible — not established','More cortical labeling in BLA can have several explanations',
 'Hypotheses only. MEF2 family regulation is not equivalent to Mef2c-specific cortical axon pruning. Harrington2016 context-dependent cortical outcomes https://pmc.ncbi.nlm.nih.gov/articles/PMC5094851/; Flavell2006 https://pubmed.ncbi.nlm.nih.gov/16484497/; Chauvet2007 doi10.1016/j.neuron.2007.10.019; Martin2015 https://elifesciences.org/articles/09395.');
table(s,[['Scenario','Candidate components','Decisive independent test'],
 ['Retained branches / synapses','MEF2 activity pathway, Pcdh10;\nglial engulfment machinery','Sparse axons at several ages +\nboutons/synapses, normalized to cortical labeling'],
 ['Altered cue response / matching','Sema3e–Plxnd1/Nrp1,\nKirrel3, cadherins, Slit/Robo','KO/WT projecting cortical neurons;\nprotein/cue-response and rescue tests'],
 ['Cell identity / labeling difference','Cortical fate, survival, reporter\nexpression or tracer uptake','Count labeled cortical neurons;\nmatch layer/type, injection and survival']
 ],{h:358,font:23,widths:[1.4,2,2.7]});
note(s,'BLA RNA availability helps choose a test. It cannot decide which scenario caused the projection phenotype.');

s=slide('A practical follow-up for the wiring project','Prioritize measurements that distinguish mechanisms',
 'Experimental proposals, not observed results. No new perturbation is inferred from these transcriptomes. Biological replicates and confirmed genotype are required. ITC-like / Rspo2-positive hypotheses must be spatially verified.');
text(s,'1. Establish the anatomical phenotype',64,210,1100,45,30,teal,true);
text(s,'Same-age KO/WT; normalize axon signal to cortical labeling and injection.\nSeparate axon branches, boutons, synapses and reporter intensity.',64,265,1130,85,25);
text(s,'2. Identify both sides of the connection',64,374,1100,45,30,teal,true);
text(s,'Spatially validate P10 candidate BLA / ITC-like cells; profile or label\nthe actual cortical projection neurons rather than using adult atlas alone.',64,430,1130,86,25);
text(s,'3. Test one candidate mechanism with a rescue',64,542,1100,45,30,teal,true);
text(s,'Choose guidance/matching versus branch retention only after these checks.',64,594,1130,55,25);

s=slide('Reproducible pipeline and data handoff','Original inputs are preserved; all important decisions are documented',
 'Sources: ANALYSIS_PLAN.md, METHODS.md, README.md, scripts/run_pipeline.R. Package versions in results/validation/. Public/private sharing choice must be respected; large data are release assets only if authorized. Analysis works independently of presentation runtime.');
table(s,[['Step','What it does','Saved evidence'],['00 / 01','Input audit → shared-gene QC → putative doublets','Manifest, QC, barcode check, classifier calls'],['02 / 03','Marker/reference checks → two maps','Predictions, ambiguous labels, RNA integrity'],['04 / 05','Library pseudobulk → sensitivity → figures','CPM, detection, TF/candidate tables'],['PowerPoint','Editable charts/tables + scientific plots','Builder, source JSON, rendered/validated deck']],{h:330,font:23,widths:[.8,2.2,2.2]});
note(s,'Start: read analysis_development_v2/README.md. Load the processed object with readRDS(); use RNA counts for reanalysis.',570);

s=slide('Bottom line for collaboration','A useful developmental resource — with explicit uncertainty',
 'This deck reports exploratory analysis, not validated biological causality. Exact sample-level values, failed baselines and sensitivity are supplied. Missing animal metadata/reference/Cell Ranger background prevent definitive inference.');
text(s,'Available now',64,216,510,50,32,teal,true);
text(s,'Four-stage maps and provisional populations\nP10 TF and candidate-gene expression\nPer-library P10/P21 source/QC checks\nReproducible R/Python pipeline',64,281,590,230,27);
text(s,'Not established',720,216,480,50,32,red,true);
text(s,'MEF2C KO transcriptional changes\nDirect MEF2 targets or TF activity\nCortical axon pruning / cue switching\nDefinitive developmental tissue fractions',720,281,500,230,27);
note(s,'The strongest next step is matched KO/WT projection-neuron data plus anatomical validation, not a causal claim from UMAP.');

const staging=path.join(ROOT,'.codex-finalizer'); await fs.mkdir(staging,{recursive:true});
const out=path.resolve(process.env.BLA_PPTX_OUTPUT || path.join(ROOT,'results','BLA_development_E18_P0_P10_P21_2026-10-06_final.pptx'));
const candidate=path.join(staging,'candidate.pptx');
await (await PresentationFile.exportPptx(p)).save(candidate);
const finalized=await finalizePresentation({workspaceDir:ROOT,candidatePath:candidate,finalPath:out,
 pythonExecutable:RUNTIME_PYTHON,
 integrityValidatorPath:path.join(SKILL_DIR,'container_tools/inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(SKILL_DIR,'container_tools/inspect_presentation_layout_geometry.py'),
 layoutArgs:['--expected-slide-size-emu','12192000,6858000','--validate-bullet-geometry','--validate-heading-fit',...[...new Set(tableOwners)].flatMap(n=>['--require-native-table-slide',String(n)])],
 explicitTotalSlideCount:18,requiredNativeTableOwnerSlides:[...new Set(tableOwners)],requiredNativeChartOwnerSlides:[...new Set(chartOwners)],
 materializeLiteralChartWorkbooks:true,fontPolicy:{basis:'design',families:[family]},verifyArtifactToolImport:true,
 receiptPath:path.join(staging,path.basename(out)+'.validation.json')});
console.log(JSON.stringify(finalized));
// Render the actual finalized file, not only the in-memory draft.
const finalDeck=await PresentationFile.importPptx(await FileBlob.load(out));
const previews=path.join(ROOT,'.build','slides'); await fs.mkdir(previews,{recursive:true});
for(let i=0;i<finalDeck.slides.items.length;i++) {
 const png=await finalDeck.export({slide:finalDeck.slides.items[i],format:'png',scale:1});
 await fs.writeFile(path.join(previews,`slide-${String(i+1).padStart(2,'0')}.png`),new Uint8Array(await png.arrayBuffer()));
}
console.log(JSON.stringify({pptx:out,slides:finalDeck.slides.items.length,previews}));
