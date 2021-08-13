#!/usr/bin/env python

"""
Must have the following environment variables defined:
* BUCKET_NAME : AWS bucket name
* PREFIX : 'latest' or other version number

"""

import os
import pip
import tempfile
import subprocess
import thermopyl.version


BUCKET_NAME = 'thermopyl.org'
if not thermopyl.version.release:
    PREFIX = 'latest'
else:
    PREFIX = thermopyl.version.short_version

if not any(d.project_name == 's3cmd' for d in pip.get_installed_distributions()):
    raise ImportError('The s3cmd package is required. try $ pip install s3cmd')
# The secret key is available as a secure environment variable
# on travis-ci to push the build documentation to Amazon S3.
with tempfile.NamedTemporaryFile('w') as f:
    f.write('''[default]
access_key = {AWS_ACCESS_KEY_ID}
secret_key = {AWS_SECRET_ACCESS_KEY}
'''.format(**os.environ))
    f.flush()

    template = ('s3cmd --guess-mime-type --config {config} '
                'sync docs/_build/ s3://{bucket}/{prefix}/')
    cmd = template.format(
            config=f.name,
            bucket=BUCKET_NAME,
            prefix=PREFIX)
    return_val = subprocess.call(cmd.split())

    # Sync index file.
    template = ('s3cmd --guess-mime-type --config {config} '
                'sync devtools/ci/index.html s3://{bucket}/')
    cmd = template.format(
            config=f.name,
            bucket=BUCKET_NAME)
    return_val = subprocess.call(cmd.split())
