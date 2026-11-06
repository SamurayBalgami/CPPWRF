#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Print a downscoped GCS access token: read-only (roles/storage.objectViewer) for one bucket, valid at most 1 h.

    tools/colab/gcs_token.py BUCKET > TOKEN_FILE

The Colab runtime gets only this token (gcloud storage --access-token-file TOKEN_FILE), never the workstation's
gcloud credentials. The source token is `gcloud auth print-access-token` of the active account; Google's STS
exchanges it with a Credential Access Boundary (google.auth.downscoped) that allows reading this bucket's objects
and nothing else. A copy that takes longer than the token's lifetime asks for a new one.
"""
import subprocess
import sys

import google.auth.credentials
import google.auth.downscoped
import google.auth.transport.requests


class GcloudToken(google.auth.credentials.Credentials):
    """The active gcloud account's access token as source credentials."""

    def refresh(self, request):
        self.token = subprocess.run(["gcloud", "auth", "print-access-token"], check=True,
                                    capture_output=True, text=True).stdout.strip()


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    rule = google.auth.downscoped.AccessBoundaryRule(
        available_resource=f"//storage.googleapis.com/projects/_/buckets/{sys.argv[1]}",
        available_permissions=["inRole:roles/storage.objectViewer"])
    cred = google.auth.downscoped.Credentials(
        source_credentials=GcloudToken(),
        credential_access_boundary=google.auth.downscoped.CredentialAccessBoundary(rules=[rule]))
    cred.refresh(google.auth.transport.requests.Request())
    print(cred.token)


if __name__ == "__main__":
    main()
