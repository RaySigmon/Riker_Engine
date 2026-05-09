# Infrastructure — Remote and Mirror Map

## Git remotes (on Pi 5 `kai`)

| Remote | URL | Purpose |
|--------|-----|---------|
| `origin` | `https://github.com/RaySigmon/Riker_Engine.git` | GitHub source of truth |
| `web` | `riker-deploy@45.32.207.82:/srv/git/riker.git` | Vultr mirror → rikerengine.quickaffordablesites.com |

The `web` remote has a post-receive hook that auto-deploys on push:
- Updates `/var/www/riker/` working tree
- Regenerates `STATE.json` and `MANIFEST.json`
- Appends to audit log

## Push protocol

Both remotes should be pushed after any commit containing data or manifests:
```bash
git push origin main   # GitHub (requires PAT in credential store)
git push web main      # Mirror (SSH key auth)
```

The post-receive hook on `web` also pushes to `origin`, so `git push web main` alone keeps both in sync. Direct `git push origin main` is still preferred for explicit control.

## Hosts

| Hostname | Hardware | Role |
|----------|----------|------|
| `kai` (this Pi 5) | RPi 5, 8GB RAM | Primary compute for validation |
| `ghost` (other Pi 5) | RPi 5, 8GB RAM | Project Vulture (separate project) |
| `qas-production` (45.32.207.82) | Vultr VPS | Mirror site, STATE.json hosting |
| RunPod (on-demand) | Cloud GPU/CPU | AD blind, cross-hardware verification |

## Mirror site

- URL: https://rikerengine.quickaffordablesites.com
- Verification: compare `STATE.json → version.commit_short` against `commits/log.txt`
- Content: full repo checkout, browsable disease day artifacts
