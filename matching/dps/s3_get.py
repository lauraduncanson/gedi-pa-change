#!/usr/bin/env python3
"""Download one s3:// object into a folder with the job's AWS role: s3_get.py s3://bucket/key dest_dir/"""
import os
import sys

import boto3

url, dest = sys.argv[1], sys.argv[2]
bucket, key = url[len("s3://"):].split("/", 1)
out = os.path.join(dest, os.path.basename(key))
boto3.client("s3").download_file(bucket, key, out)
print(f"downloaded {url} ({os.path.getsize(out):,} bytes)")
