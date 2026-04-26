# Session Receipt: <YYYY-MM-DD> — <task name>

**Started:** <ISO timestamp>
**Completed:** <ISO timestamp>
**Executor:** Kai
**Instructions source:** <reference to the block from Claude / Cody>

---

## 1. Commands executed

(Verbatim, in order, with timestamps for any command spanning >1 minute.
This is the raw shell record, not a narrative.)

## 2. Files created or modified

(Output of: find <relevant_paths> -newer <start_marker> -type f -exec sha256sum {} \; -exec stat -c '%y  %s  %n' {} \;
Group by path. Include sha256, byte size, and modification time for each file.)

## 3. Git state at end

- HEAD: <commit hash>
- Branch: <branch>
- `git status` output (verbatim):
- `git log --oneline -10` (verbatim):

## 4. Output verification

(For every numerical claim made in the report-back summary, show the raw command
and output that backs it. Format:

CLAIM: "<claim>"
VERIFIED BY:
$ <command>
<output>
INTERPRETATION: <how the output supports the claim>

If a claim cannot be verified by a deterministic command, mark as UNVERIFIED.)

## 5. Decisions or fixes outside original instructions

(If the original instructions did not anticipate a situation and Kai made a
choice, log it here. For each entry:
- Situation encountered
- Decision made
- Alternatives not chosen
- Rationale
- Was the change reversible? (yes/no)

If empty, write "None.")

## 6. Notes for next session

(Anything the next instance of Kai or Claude should know to pick up cleanly.
Empty is fine.)
