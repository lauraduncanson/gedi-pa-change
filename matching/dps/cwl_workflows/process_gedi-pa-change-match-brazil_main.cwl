cwlVersion: v1.2
$graph:
- class: Workflow
  label: gedi-pa-change-match-brazil
  doc: Veronika Leitold's Brazil V3 propensity matching (step 4) for one Brazil biome
    region, with GEDI L4B in the propensity model (with_l4b, her V3 setting) or without
    it (no_l4b). Reads her MATCHING_points_V3 inputs; writes MATCHING_results_V3/
    to the job output.
  id: gedi-pa-change-match-brazil
  inputs:
    region:
      doc: Brazil biome region code (BAe, BAw, BCa, BCe, BMa, BPp, BPt)
      label: Region
      type: string
    covariates:
      doc: with_l4b (propensity model as in her V3 run), no_l4b (gedi_l4b removed,
        all else unchanged) or l4b_classes (gedi_l4b removed; exact matching within
        50 Mg/ha L4B classes)
      label: Covariate set
      type: string
    n_cores:
      doc: Parallel workers over PAs, e.g. the number of vCPUs of the queue
      label: Number of cores
      type: string
    part:
      doc: Part i/n of the region's PAs to run (1/1 = all). Splitting does not change
        any match.
      label: Part
      type: string
    input_prefix:
      doc: s3:// folder with the step 1-3 inputs; her Brazil MATCHING_points_V3 folder,
        or the matching_inputs folder from matching/prep_global_inputs.py (then <prefix>/<region>/)
      label: Input folder
      type: string
  outputs:
    out:
      type: Directory
      outputSource: process/outputs_result
  steps:
    process:
      run: '#main'
      in:
        region: region
        covariates: covariates
        n_cores: n_cores
        part: part
        input_prefix: input_prefix
      out:
      - outputs_result
- class: CommandLineTool
  id: main
  requirements:
    DockerRequirement:
      dockerPull: ghcr.io/lauraduncanson/gedi-pa-change-match-brazil:main
    NetworkAccess:
      networkAccess: true
    ResourceRequirement:
      ramMin: 16384
      coresMin: 4
      outdirMax: 8192
  baseCommand: run_matching.sh
  inputs:
    region:
      type: string
      inputBinding:
        position: 1
        prefix: --region
    covariates:
      type: string
      inputBinding:
        position: 2
        prefix: --covariates
    n_cores:
      type: string
      inputBinding:
        position: 3
        prefix: --n_cores
    part:
      type: string
      inputBinding:
        position: 4
        prefix: --part
    input_prefix:
      type: string
      inputBinding:
        position: 5
        prefix: --input_prefix
  outputs:
    outputs_result:
      outputBinding:
        glob: ./output*
      type: Directory
s:author:
- class: s:Person
  s:name: lauraiduncanson
s:contributor:
- class: s:Person
  s:name: lauraiduncanson, leitoldv
s:citation: https://github.com/lauraduncanson/gedi-pa-change.git
s:codeRepository: https://github.com/lauraduncanson/gedi-pa-change.git
s:commitHash: 49a5a49fef01b519a65a34b8ee77d07c77653d4e
s:dateCreated: 2026-11-26
s:license: https://raw.githubusercontent.com/lauraduncanson/gedi-pa-change/refs/heads/main/LICENSE
s:softwareVersion: 1.0.0
s:version: v6
s:releaseNotes: None
s:keywords: ogc, gedi, protected areas, matching
$namespaces:
  s: https://schema.org/
$schemas:
- https://raw.githubusercontent.com/schemaorg/schemaorg/refs/heads/main/data/releases/9.0/schemaorg-current-http.rdf
