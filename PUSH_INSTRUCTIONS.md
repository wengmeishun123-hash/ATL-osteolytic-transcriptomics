# 🔐 GitHub Push Instructions (run on your local machine, not in any AI sandbox)

## ⚠ Security note
**Do NOT paste your GitHub password into AI chats / scripts**. Use either:
- GitHub Personal Access Token (PAT) with `repo` scope, stored in your local Windows Credential Manager via `git credential-manager`, or
- GitHub CLI (`gh auth login`) with browser-based device-code flow (recommended)

Both methods avoid sending your password to any third party.

## Step 1 — Install GitHub CLI (one-time)

Download from https://cli.github.com/ (Windows MSI installer).

After install, open PowerShell:

```powershell
gh --version           # confirm install
gh auth login          # interactive: browser → Login with your GitHub account
```

When prompted:
- "Where do you use GitHub?" → `GitHub.com`
- "Preferred protocol for Git operations?" → `HTTPS`
- "Authenticate Git with your GitHub credentials?" → `Yes`
- "How would you like to authenticate?" → `Login with a web browser`
- Browser opens → enter your GitHub username (wengmeishun123@icloud.com) and password → authorize device

GitHub CLI will store the token securely in Windows Credential Manager. **Your password is never written to disk.**

## Step 2 — Initialize the repository

```powershell
cd "D:\ATL research\analysis\github_repo"
git init
git branch -M main
git config user.name "Meishun Weng"
git config user.email "wengmeishun2026@163.com"
```

## Step 3 — Stage and commit

```powershell
git add .
git status                       # review what will be committed (especially check .gitignore is excluding clinical files)
git commit -m "Initial commit: ATL osteolytic transcriptomics analysis pipeline (manuscript v2.0)"
```

## Step 4 — Create remote repository (one command)

```powershell
gh repo create wengmeishun123/ATL-osteolytic-transcriptomics --public `
   --description "Single-cell and bulk transcriptomic dissection of osteolytic programs in adult T-cell leukemia/lymphoma" `
   --homepage "https://www.frontiersin.org/journals/immunology" `
   --source=. --push
```

This will:
1. Create the repository on github.com
2. Add it as `origin` remote
3. Push the `main` branch
4. Configure default branch

After completion: visit https://github.com/wengmeishun123/ATL-osteolytic-transcriptomics

## Step 5 — Add Zenodo integration (for DOI on release)

1. Go to https://zenodo.org/account/settings/github/ and authorize Zenodo
2. Toggle "On" for the `wengmeishun123/ATL-osteolytic-transcriptomics` repository
3. On GitHub, create release: `gh release create v1.0 --title "v1.0 — manuscript submission" --notes "Initial public release accompanying Frontiers in Immunology submission"`
4. Zenodo will automatically create a DOI for v1.0
5. Replace `[ZENODO_DOI]` placeholder in manuscript with the assigned DOI (format: 10.5281/zenodo.XXXXXXX)

## Step 6 — Update README badge

After Zenodo assigns a DOI, edit `README.md` and replace:
```
[![Zenodo](https://img.shields.io/badge/Zenodo-DOI%20pending-blue)](https://zenodo.org/)
```
with:
```
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

Then commit + push:
```powershell
git add README.md
git commit -m "Add Zenodo DOI badge"
git push
```

## Troubleshooting

### `gh auth login` fails with "client error"
- Disable any VPN / proxy that intercepts HTTPS traffic
- Try `gh auth login --hostname github.com --git-protocol https`

### Repository name conflict
If `wengmeishun123/ATL-osteolytic-transcriptomics` already exists:
- Either delete the old one: `gh repo delete wengmeishun123/ATL-osteolytic-transcriptomics --yes`
- Or pick a different name in Step 4

### Push rejected: large file
GitHub rejects files >100 MB. Check `.gitignore` excludes all `.rds` / `.tar.gz` / `.gctx` files. If something slipped through:
```powershell
git rm --cached <large_file>
git commit -m "Remove large file"
git push
```

For very large data (e.g., GSE195674 .rds), use Git LFS or just keep them out of the repo (per .gitignore).
