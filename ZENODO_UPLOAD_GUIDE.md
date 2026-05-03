# 📦 Zenodo Upload Guide

This guide walks you through depositing the GitHub repository to Zenodo to obtain a citable DOI for the manuscript's Data Availability statement.

## Method 1 — GitHub-Zenodo Integration (Recommended, fully automated)

This is the easiest route: every GitHub release automatically generates a Zenodo DOI.

### Step 1.1 — Authorize Zenodo
1. Go to https://zenodo.org/account/settings/github/
2. Click **Connect** to GitHub (or **Login with GitHub** if not signed in)
3. Authorize Zenodo to access your public repositories

### Step 1.2 — Toggle the repository
1. After authorization, you'll see a list of your GitHub repositories
2. Find `wengmeishun123/ATL-osteolytic-transcriptomics`
3. Toggle the switch to **On**
4. Zenodo will now watch this repository for new releases

### Step 1.3 — Create the first release on GitHub
```powershell
cd "D:\ATL research\analysis\github_repo"

# Create v1.0.0 release
gh release create v1.0.0 `
   --title "v1.0.0 — Manuscript submission to Frontiers in Immunology" `
   --notes "Initial public release of the analysis code and gene panels accompanying the manuscript 'Single-cell and bulk transcriptomic dissection of osteolytic programs in ATL' (Weng et al., 2026)."
```

### Step 1.4 — Verify Zenodo deposit
1. Go to https://zenodo.org/account/settings/github/
2. The new release should appear with a DOI assigned
3. Copy the DOI (format: `10.5281/zenodo.XXXXXXX`)

### Step 1.5 — Update manuscript and README
Replace the placeholder `[ZENODO_DOI]` in your manuscript with the assigned DOI.

In `README.md`, replace:
```
[![Zenodo](https://img.shields.io/badge/Zenodo-DOI%20pending-blue)](https://zenodo.org/)
```
with:
```
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

Commit and push:
```powershell
git add README.md
git commit -m "Add Zenodo DOI badge"
git push
```

## Method 2 — Manual Zenodo Upload (alternative)

Use this if you prefer not to integrate with GitHub.

### Step 2.1 — Create Zenodo account
- Go to https://zenodo.org/signup
- Use your institutional email if possible

### Step 2.2 — Create a new upload
- Click **+ New Upload** at https://zenodo.org/deposit/new
- Upload type: **Software**
- License: **MIT**

### Step 2.3 — Fill metadata
Use the values from `.zenodo.json` (already in this repo):
- Title, Authors (with ORCIDs), Affiliation
- Keywords
- Publication date
- Description (paste from `.zenodo.json`)

### Step 2.4 — Upload the repository as a ZIP
```powershell
cd "D:\ATL research\analysis"
Compress-Archive -Path "github_repo\*" -DestinationPath "ATL-osteolytic-transcriptomics-v1.0.0.zip"
```

Upload `ATL-osteolytic-transcriptomics-v1.0.0.zip` via the Zenodo web UI.

### Step 2.5 — Publish
Click **Publish**. The DOI is assigned within seconds.

## Method 3 — Zenodo Sandbox (for testing first)

If you want to test the deposit without creating a permanent record:
- Use https://sandbox.zenodo.org/ instead of https://zenodo.org/
- Sandbox DOIs are not citable, but the workflow is identical
- Once tested, repeat with the production Zenodo

## Updating subsequent versions

Each new GitHub release (e.g., `v1.1.0` after manuscript revision) automatically creates a new Zenodo version with a new DOI; the **concept DOI** (master) always points to the latest version. Cite the **concept DOI** in your paper for stable referencing across versions.

## After deposit — what to update

1. **Manuscript** (`ATL_Manuscript_v2.0_FINAL.docx`): Section 5 Declarations → Data and code availability paragraph → replace `[ZENODO_DOI]` with the DOI
2. **GitHub README.md**: replace pending badge with DOI badge
3. **CITATION.cff**: optionally add `doi: 10.5281/zenodo.XXXXXXX`
4. **Cover letter**: optionally mention Zenodo DOI in the data availability sentence

## Troubleshooting

### "Zenodo: GitHub repository not visible"
- Ensure repository is **public**
- Re-authorize Zenodo at https://zenodo.org/account/settings/github/

### "Release created but no DOI assigned"
- Wait 5–10 minutes; Zenodo polls GitHub asynchronously
- Check https://zenodo.org/account/settings/github/ for status
- If still missing, try creating a fresh release with a different version tag (e.g., `v1.0.1`)

### "DOI badge shows wrong version"
- The DOI badge automatically updates with new releases
- For a specific version, use the version-specific DOI (e.g., `10.5281/zenodo.XXXXXXX`)
- For the always-latest version, use the concept DOI (e.g., `10.5281/zenodo.YYYYYYY`)

