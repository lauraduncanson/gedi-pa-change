cwlVersion: v1.2
$graph:
- class: Workflow
  label: gedi-pa-change-extract-tile
  doc: Extracts quality-filtered, power-beam, growing-season GEDI shots for one 1-degree
    tile of gedi-tiled-v3 and aggregates them to 1 km EASE2 cells for T1 (2019-2021)
    and T2 (2024-2026).
  id: gedi-pa-change-extract-tile
  inputs:
    tile:
      doc: Tile id in gedi-tiled-v3, e.g. S13_E028
      label: Tile id
      type: string
  outputs:
    out:
      type: Directory
      outputSource: process/outputs_result
  steps:
    process:
      run: '#main'
      in:
        tile: tile
      out:
      - outputs_result
- class: CommandLineTool
  id: main
  requirements:
    DockerRequirement:
      dockerPull: ghcr.io/lauraduncanson/gedi-pa-change-extract-tile:main
    NetworkAccess:
      networkAccess: true
    ResourceRequirement:
      ramMin: 8192
      coresMin: 2
      outdirMax: 4096
  baseCommand: extract_tile.py
  inputs:
    tile:
      type: string
      inputBinding:
        position: 1
        prefix: --tile
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
  s:name: lauraiduncanson
s:citation: https://github.com/lauraiduncanson/gedi-pa-change.git
s:codeRepository: https://github.com/EDIT_ME/gedi-pa-change.git
s:commitHash: 088f4d1231ecc86b8088c4b9e615438969686718
s:dateCreated: 2026-11-27
s:license: https://raw.githubusercontent.com/EDIT_ME/gedi-pa-change/refs/heads/main/LICENSE
s:softwareVersion: 1.0.0
s:version: drift-v1
s:releaseNotes: None
s:keywords: ogc, gedi, protected areas
$namespaces:
  s: https://schema.org/
$schemas:
- https://raw.githubusercontent.com/schemaorg/schemaorg/refs/heads/main/data/releases/9.0/schemaorg-current-http.rdf
