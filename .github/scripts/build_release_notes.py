#!/usr/bin/env python3
"""Build rich GitHub release notes from merged pull requests.

For each merged PR since the previous tag we extract the
``## Proposed Changes`` and ``## Usage or Results`` sections from the
PR body (matching the wording of `.github/PULL_REQUEST_TEMPLATE.md`)
and group them under category headings determined by labels.

Output is markdown printed to stdout, suitable to feed into

    gh release create "$TAG" --notes-file -
or
    gh release edit "$TAG" --notes-file release-notes.md --draft=false

Requires:
    - python 3.9+
    - `gh` CLI authenticated (set GH_TOKEN in CI)
    - run from inside a git checkout with `fetch-depth: 0`

Inputs (CLI):
    --tag TAG         the tag being published (e.g. v0.3.0). Required.
    --previous TAG    previous tag to diff against. Auto-detected if omitted.
    --repo OWNER/REPO defaults to env GITHUB_REPOSITORY.

Stability contract:
    The PR template's section headings (`## Proposed Changes` and
    `## Usage or Results`) are an implicit input to this script.
    If the template is renamed, update SECTIONS below.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from typing import Optional


# ---------------------------------------------------------------------
#  Configuration
# ---------------------------------------------------------------------

# Categories for grouping. Order matters — labels are matched in this
# order, and a PR is placed in the FIRST category whose label set
# contains one of the listed labels. PRs with no matching label fall
# into UNCATEGORISED.
CATEGORIES = [
    ("⚠️ Breaking Changes", ["breaking"]),
    ("🚀 Features",          ["enhancement", "feature"]),
    ("🐛 Bug Fixes",         ["bug", "fix"]),
    ("⚡ Performance",       ["performance"]),
    ("📖 Documentation",     ["documentation", "docs"]),
    ("🧰 Maintenance",       ["chore", "refactor", "ci"]),
]
UNCATEGORISED = "🔧 Other Changes"

# Labels that exclude a PR from the release notes entirely.
EXCLUDE_LABELS = {"skip-changelog", "wontfix", "duplicate", "invalid"}

# Section headers we extract from each PR body. Each tuple is:
#     (header in PR template, display label in release notes)
SECTIONS = [
    ("Proposed Changes", "What changed"),
    ("Usage or Results", "How to use"),
]


# ---------------------------------------------------------------------
#  Subprocess helpers
# ---------------------------------------------------------------------

def gh(*args: str) -> str:
    """Run `gh <args...>` and return stdout. Raises on non-zero."""
    return subprocess.check_output(["gh", *args], text=True)


def gh_json(*args: str):
    return json.loads(gh(*args))


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


# ---------------------------------------------------------------------
#  Tag / PR discovery
# ---------------------------------------------------------------------

def find_previous_tag(current_tag: str) -> Optional[str]:
    """Return the most recent tag before *current_tag* (by creator date)."""
    try:
        tags = git("tag", "--sort=-creatordate").splitlines()
    except subprocess.CalledProcessError:
        return None
    for t in tags:
        if t and t != current_tag:
            return t
    return None


def merged_prs_since(repo: str, previous_tag: Optional[str]):
    """Return a list of PR dicts merged since *previous_tag*'s commit date.

    If *previous_tag* is None this returns all merged PRs (used for the
    very first release).
    """
    search = "is:merged"
    if previous_tag:
        try:
            prev_date = git("log", "-1", "--format=%aI", previous_tag)
            search = f"merged:>{prev_date}"
        except subprocess.CalledProcessError:
            pass

    return gh_json(
        "pr", "list",
        "--repo", repo,
        "--state", "merged",
        "--base", "main",
        "--limit", "200",
        "--search", search,
        "--json", "number,title,author,body,labels,url,mergedAt",
    )


# ---------------------------------------------------------------------
#  Section extraction
# ---------------------------------------------------------------------

# Phrases that indicate a section was left at its template default.
PLACEHOLDER_PHRASES = [
    r"what did you change\??",
    r"-\s*write here",
    r"#\s*Example or verification script",
    r"_no description provided\._",
    r"placeholder",
]

PLACEHOLDER_RE = re.compile("|".join(PLACEHOLDER_PHRASES), re.IGNORECASE)
HTML_COMMENT_RE = re.compile(r"<!--.*?-->", re.DOTALL)


def strip_decoration(text: str) -> str:
    """Remove HTML comments and surrounding whitespace from `text`."""
    text = HTML_COMMENT_RE.sub("", text)
    return text.strip()


def is_meaningful(text: str) -> bool:
    """Return True if `text` looks like real content rather than placeholder.

    Strips HTML comments, code-fence markers and known placeholder
    phrases, then checks that something other than whitespace remains.
    """
    cleaned = HTML_COMMENT_RE.sub("", text)
    cleaned = re.sub(r"^```\w*$|^```$", "", cleaned, flags=re.MULTILINE)
    cleaned = PLACEHOLDER_RE.sub("", cleaned)
    return bool(cleaned.strip())


def extract_section(body: str, header: str) -> str:
    """Extract the body of `## <header>` until the next `## ` heading.

    Returns an empty string if the section is missing or contains only
    placeholder content. HTML comments inside the section are stripped
    so template hint comments do not leak into release notes.
    """
    if not body:
        return ""
    body = body.replace("\r\n", "\n").replace("\r", "\n")
    pattern = re.compile(
        r"^##\s+" + re.escape(header) + r"\s*\n(.*?)(?=^##\s|\Z)",
        re.DOTALL | re.MULTILINE,
    )
    m = pattern.search(body)
    if not m:
        return ""
    text = m.group(1)
    if not is_meaningful(text):
        return ""
    return strip_decoration(text)


# ---------------------------------------------------------------------
#  Categorisation + rendering
# ---------------------------------------------------------------------

def categorise(pr) -> str:
    label_names = {l["name"].lower() for l in pr.get("labels", [])}
    for title, labels in CATEGORIES:
        if any(lab in label_names for lab in labels):
            return title
    return UNCATEGORISED


def excluded(pr) -> bool:
    label_names = {l["name"].lower() for l in pr.get("labels", [])}
    return bool(label_names & EXCLUDE_LABELS)


def render_pr(pr) -> str:
    """Render one PR as a markdown sub-section."""
    title = pr["title"]
    number = pr["number"]
    author = pr.get("author") or {}
    login = author.get("login", "ghost")
    url = pr["url"]
    body = pr.get("body") or ""

    lines = [f"### #{number} — {title} (@{login})", ""]

    sections_present = False
    for header, display in SECTIONS:
        section = extract_section(body, header)
        if not section:
            continue
        sections_present = True
        # Use a collapsible details block so long Usage code blocks
        # don't blow out the release notes view.
        lines.append(f"<details><summary><b>{display}</b></summary>")
        lines.append("")
        lines.append(section)
        lines.append("")
        lines.append("</details>")
        lines.append("")

    if not sections_present:
        lines.append(f"_No detailed notes provided. See {url}._")
        lines.append("")

    return "\n".join(lines)


def build_notes(repo: str, tag: str, previous: Optional[str]) -> str:
    prs = [pr for pr in merged_prs_since(repo, previous) if not excluded(pr)]

    out: list[str] = []
    out.append(f"## What's Changed in {tag}")
    out.append("")

    if not prs:
        out.append("_No merged pull requests in this release._")
        out.append("")
    else:
        groups: dict[str, list] = {title: [] for title, _ in CATEGORIES}
        groups[UNCATEGORISED] = []
        for pr in prs:
            groups[categorise(pr)].append(pr)

        for title in [c[0] for c in CATEGORIES] + [UNCATEGORISED]:
            items = groups.get(title) or []
            if not items:
                continue
            out.append(f"## {title}")
            out.append("")
            for pr in items:
                out.append(render_pr(pr))
                out.append("")

    if previous:
        out.append(
            f"**Full Changelog**: "
            f"https://github.com/{repo}/compare/{previous}...{tag}"
        )
        out.append("")

    return "\n".join(out).rstrip() + "\n"


# ---------------------------------------------------------------------
#  Entrypoint
# ---------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tag", required=True, help="tag being published, e.g. v0.3.0")
    p.add_argument("--previous", default=None,
                   help="previous tag (auto-detected if omitted)")
    p.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY"),
                   help="owner/repo (defaults to GITHUB_REPOSITORY env var)")
    args = p.parse_args()

    if not args.repo:
        sys.exit("error: --repo or GITHUB_REPOSITORY required")

    previous = args.previous or find_previous_tag(args.tag)
    print(build_notes(args.repo, args.tag, previous))


if __name__ == "__main__":
    main()
