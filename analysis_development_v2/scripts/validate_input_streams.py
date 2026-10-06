"""Verify RDS compression streams without changing inputs."""
import gzip
import json
import os
from pathlib import Path

paths = [Path(os.environ.get("BLA_OLD_RDS", r"K:\scRNA_BLA_phd\hansol\seu_harm_qc.rds")),
         Path(os.environ.get("BLA_PRIOR_RDS", r"K:\scRNA_BLA_phd\P10_P21_amygdala_2026-10\B06_P10_P21_annotated.rds")),
         Path(os.environ.get("BLA_NEW_RDS", r"K:\scRNA_BLA_phd\hansol2_combined_seurat.rds"))]
results = []
for path in paths:
    result = {"file": str(path), "bytes": path.stat().st_size}
    try:
        with gzip.open(path, "rb") as stream:
            size = 0
            while block := stream.read(4 * 1024 * 1024):
                size += len(block)
        result.update(gzip_stream_valid=True, uncompressed_bytes=size)
    except (OSError, EOFError) as error:
        result.update(gzip_stream_valid=False, error=str(error))
    print(json.dumps(result))
    results.append(result)
out=Path(os.environ.get("BLA_V2_ROOT","analysis_development_v2"))/"results/validation"
out.mkdir(parents=True,exist_ok=True)
(out/"input_compression_streams.json").write_text(json.dumps(results,indent=2),encoding="utf-8")
if not all(x["gzip_stream_valid"] for x in results): raise SystemExit("An input compression stream failed validation")
