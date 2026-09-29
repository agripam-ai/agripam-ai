# Version 1.0.0 release and DOI checklist

Do not create the final tag until the competition version has been frozen and
the checks below pass on the exact commit intended for publication.

## 1. Confirm release metadata

- Replace `AgriPAM-AI team` in `CITATION.cff` and `.zenodo.json` with the final
  creator list in publication order.
- Add affiliations and ORCID identifiers where available.
- Confirm the MIT licence and public contact address.
- Confirm that sensitive, private or restricted laboratory data are absent.
- Replace `release pending` in `CHANGELOG.md` with the release date.

## 2. Freeze and verify the competition version

```bash
git status --short
python -m unittest discover -s tests
docker build -t agripam-ai:1.0.0 .
git diff --check
```

Record the test count, application URL and release commit in the competition
report. Test the public example without authentication in a private browser
window and confirm that the exported files reproduce the documented example.

## 3. Connect GitHub to Zenodo

1. Sign in to Zenodo with the GitHub account that administers
   `agripam-ai/agripam-ai`.
2. Open the GitHub integration page in Zenodo.
3. Enable archiving for the `agripam-ai/agripam-ai` repository.
4. Confirm that `CITATION.cff`, `.zenodo.json`, `LICENSE`, `README.md` and the
   example inputs are present on the default branch.

Enabling the connection does not archive ordinary commits. Zenodo creates a
new deposit when a GitHub release is published.

## 4. Create the immutable release

After the final commit has been pushed and reviewed:

```bash
git tag -a v1.0.0 -m "AgriPAM-AI v1.0.0 competition and manuscript release"
git push origin main
git push origin v1.0.0
```

On GitHub, create a release from `v1.0.0`, use the corresponding section of
`CHANGELOG.md` as the release notes, and publish it. Do not move or reuse this
tag after Zenodo has archived it; corrections require a new version.

## 5. Verify the Zenodo archive

- Confirm that the deposit contains the intended commit and files.
- Check title, creators, version, licence, keywords and related URLs.
- Reserve or publish the record as appropriate for the competition schedule.
- Record both identifiers:
  - the **version DOI**, which identifies exactly `v1.0.0`;
  - the **concept DOI**, which resolves to the latest archived version.

Use the version DOI for the reproducible manuscript release. Use the concept
DOI on the project page when readers should reach the latest release.

## 6. Update citation surfaces

After Zenodo assigns the DOI:

- add `doi:` to `CITATION.cff` using the version DOI;
- add the DOI badge and citation text to `README.md`;
- replace the provisional repository statement in the manuscript;
- cite the software release in the manuscript reference list;
- confirm that the manuscript names the exact version and commit.

